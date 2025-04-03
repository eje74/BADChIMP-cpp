import numpy as np
import pyvista as pv
from numba import njit
from numba import prange
from matplotlib import pyplot as plt
from vtklb import vtklb
from clean_geometry import bwlabel3D_set

# ======================================================================== Data format
# Solid geometry:
#   element > 0: solid
# Fluid geometry:
#   element > 0: wetting

# ######################################################################## FUNCTIONS
# ======================================================================== Remove non-perculating domains
def remove_non_perculating_domains(geo, marker=(1,)):
    """
    Remove the non percolating clusters, in the z-direction,
    from the geometry.

    Paramters
    ---------
    geo : nd.array(int)
        3d matrix defining the geometry
    
    marker: tuple
        Loist of lables in `geo` that make up the domains

    Return
    ------
    numpy.array(int)
        Cleand geometry fluid=1, solid=0 
    """
    # ------------------------------------------------------------------------ neighborhood
    cc = np.ones((19, 3), dtype=np.int64)
    cc[:10, :] = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 0],
        [1, 0, 1],
        [0, 1, 1],
        [-1, 1, 0],
        [-1, 0, 1],
        [0, -1, 1],
        [0, 0, 0]
        ]
    cc[:-10:-1, :] = -cc[:9, :]

    # ======================================================================== Domains in geo
    domains, max_label = bwlabel3D_set(cc, geo, (1,))

    # ========================================================================== Remove non-percolating pore space
    percolating_domains = set(domains[:, :, 0].flatten()) & set(domains[:, :, -1].flatten())
    non_percolating_domains = set(range(1,max_label+1)) - percolating_domains
    for n in non_percolating_domains:
        domains[domains==n] = 0
    domains[domains>0] = 1
    return domains

# ======================================================================== Remove non-perculating domains
def domain_decomposition_for_mpi(geo_mpi, nn=(2,)*3):
    nx, ny, nz = geo_mpi.shape
    xbe = tuple((x[0], x[-1]+1) for x in np.array_split(np.arange(nx), nn[0]))
    ybe = tuple((x[0], x[-1]+1) for x in np.array_split(np.arange(ny), nn[1]))
    zbe = tuple((x[0], x[-1]+1) for x in np.array_split(np.arange(nz), nn[2]))

    nproc = 0
    for x0, x1 in xbe:
        for y0, y1 in ybe:
            for z0, z1 in zbe:
                if np.any(geo_mpi[x0:x1, y0:y1, z0:z1]>0):
                    nproc += 1
                    geo_mpi[x0:x1, y0:y1, z0:z1] *= nproc
    return geo_mpi, nproc

def object_name(obj, namespace):
    return [n for n in namespace if namespace[n] is obj]    

# ======================================================================== Path to data
# ------------------------------------------------------------------------ root path 
pathinput = r"/home/AD.NORCERESEARCH.NO/esje/OneDrive/NORCE/CSSR/RelPerm LB LS/"
# ------------------------------------------------------------------------ system
pathinput += r"Castlegate_Tow20_LVC_Oil/Sat_control/PrimaryDrainage/"

# ======================================================================== Rock data
# ------------------------------------------------------------------------ load data
pore = np.load(pathinput + r"CG_PoreSolid_200x200x200_SDF_PD.npy")
#pore = pore[:50, :50, :50] # Test size
# ------------------------------------------------------------------------ generate geometry
geo = np.ones(pore.shape, dtype=np.int32)
geo[pore>0] = 0

# ======================================================================== Fluid phase data
# ------------------------------------------------------------------------ load data
fluid = np.load(pathinput + r"CG_NWP_Sw0746_200x200x200_SDF_PD.npy")
#fluid = fluid[:50, :50, :50]
# ------------------------------------------------------------------------ wetting generate geometry
wphase = np.ones_like(geo)
wphase[fluid<=0] = 0
wphase[geo==0] = 0
# ------------------------------------------------------------------------ non-wetting generate geometry
nwphase = np.ones_like(geo)
nwphase[fluid>0] = 0
nwphase[geo==0] = 0

if np.all( geo == (wphase + nwphase) ):
    print("All is well")

# ======================================================================== Remove non-percolating
geo = remove_non_perculating_domains(geo)
nwphase = remove_non_perculating_domains(nwphase)
wphase = remove_non_perculating_domains(wphase)

# ======================================================================== Setup mpi MPI
print(np.any(geo>0))
for m in (geo, nwphase, wphase):
    m, n = domain_decomposition_for_mpi(m)
    print(object_name(m, globals())[0], ": ", n)

# ======================================================================== Inputfile badchimp
# ------------------------------------------------------------------------ rock
vtk = vtklb(geo, "D3Q19", "xyz", path="/home/AD.NORCERESEARCH.NO/esje/GitHub/BADChIMP-cpp/input/mpi/")


# ------------------------------------------------------------------------ Plot rock
grid = pv.ImageData()
grid.dimensions = geo.shape
grid.origin = (0, 0, 0)
grid.spacing = (1, 1, 1)
grid.point_data["geo"] = geo.flatten(order="F")
#grid.plot(show_edges=False, opacity="linear")
grid.plot(volume=True, opacity=[0, 0.5], shade=False)

grid = pv.ImageData()
grid.dimensions = geo.shape
grid.origin = (0, 0, 0)
grid.spacing = (1, 1, 1)
grid.point_data["geo"] = nwphase.flatten(order="F")
#grid.plot(show_edges=False, opacity="linear")
grid.plot(volume=True, opacity=[0, 0.5], shade=False)

grid = pv.ImageData()
grid.dimensions = geo.shape
grid.origin = (0, 0, 0)
grid.spacing = (1, 1, 1)
grid.point_data["geo"] = wphase.flatten(order="F")
#grid.plot(show_edges=False, opacity="linear")
grid.plot(volume=True, opacity=[0, 0.5], shade=False)
