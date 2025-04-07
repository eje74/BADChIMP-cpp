import numpy as np
import pyvista as pv
import os
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
def domain_decomposition_for_mpi(geo_mpi, num_proc=(2,)*3):
    """
    Domain decomposition of geometries using a cubical grid
    defined by `num_proc` that gives the linear number of 
    bins for each direction.

    Paramaters
    ----------
    geo_mpi: np.ndarray[int]
        Binary geometry images where it is assumed that
        0=solid and 1=fluid

    num_proc: list-like, default: (2, 2, 2)
        Number of bins in each spatial direction so that
        `prod(num_proc) is the maximum number of domains

    Return
    ------
    np.ndarray[int]:
        Domain decomposed geometry where the domains are
        enuerated from 1 up to the total number of domains

    int:
        Total number of domains

    """
    nx, ny, nz = geo_mpi.shape
    xbe = tuple((x[0], x[-1]+1) for x in np.array_split(np.arange(nx), num_proc[0]))
    ybe = tuple((x[0], x[-1]+1) for x in np.array_split(np.arange(ny), num_proc[1]))
    zbe = tuple((x[0], x[-1]+1) for x in np.array_split(np.arange(nz), num_proc[2]))

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
#pathinput = r"/home/AD.NORCERESEARCH.NO/esje/OneDrive/NORCE/CSSR/RelPerm LB LS/"
#pathinput = r"/home/AD.NORCERESEARCH.NO/esje/OneDrive/NORCE/CSSR/RelPerm LB LS/"
# need to do a "tilde expansion", that is, insert the path for '~'
pathinput = os.path.expanduser(r"~/OneDrive/NORCE/CSSR/RelPerm LB LS/") 
# ------------------------------------------------------------------------ system
pathinput += r"Castlegate_Tow20_LVC_Oil/Sat_control/PrimaryDrainage/"

# ======================================================================== Rock data
# ------------------------------------------------------------------------ load data
pore = np.load(pathinput + r"CG_PoreSolid_200x200x200_SDF_PD.npy")
# ------------------------------------------------------------------------ generate geometry
geo = np.ones(pore.shape, dtype=np.int32)
geo[pore>0] = 0

# ======================================================================== Fluid phase data
# ------------------------------------------------------------------------ load data
fluid = np.load(pathinput + r"CG_NWP_Sw0746_200x200x200_SDF_PD.npy")
# ------------------------------------------------------------------------ wetting generate geometry
wphase = np.ones_like(geo)
wphase[fluid<=0] = 0
wphase[geo==0] = 0
# ------------------------------------------------------------------------ non-wetting generate geometry
nwphase = np.ones_like(geo)
nwphase[fluid>0] = 0
nwphase[geo==0] = 0
# ------------------------------------------------------------------------ test
if np.all( geo == (wphase + nwphase) ):
    print("All is well")

# ======================================================================== Remove non-percolating
geo_org = remove_non_perculating_domains(geo)
nwphase_org = remove_non_perculating_domains(nwphase)
wphase_org = remove_non_perculating_domains(wphase)
# ------------------------------------------------------------------------ Number of fluid nodes in geos 
print("sum geo = ", np.sum(geo_org))

# ======================================================================== Boundary conditions
# ------------------------------------------------------------------------ Extend the geometry
new_shape = tuple(x+2 for x in geo_org.shape)
geo = np.zeros(new_shape, dtype=np.int64)
geo[1:-1, 1:-1, 1:-1] = geo_org
del geo_org 

# ------------------------------------------------------------------------ Tag inlet and outlet 
geo_tag = np.zeros_like(geo)
geo_tag[1:-1, 1:-1, 0] = 1
geo_tag[1:-1, 1:-1, -1] = 1 

# ======================================================================== Setup mpi MPI
print(np.any(geo>0))
for m in (geo, nwphase, wphase):
    m, n = domain_decomposition_for_mpi(m)
    print(object_name(m, globals())[0], ": ", n)

# ======================================================================== Inputfile badchimp
# ------------------------------------------------------------------------ rock
#vtk = vtklb(geo, "D3Q19", "xyz", path="/home/AD.NORCERESEARCH.NO/esje/GitHub/BADChIMP-cpp/input/mpi/")
#vtk = vtklb(geo, "D3Q19", "", path="/home/AD.NORCERESEARCH.NO/esje/GitHub/BADChIMP-cpp/input/mpi/")
vtk = vtklb(geo, "D3Q19", "", path="../../../input/mpi/")
# ------------------------------------------------------------------------ boundary tag
vtk.append_data_set("geo_tag", geo_tag)

if False:

    # ------------------------------------------------------------------------ Plot rock
    grid = pv.ImageData()
    grid.dimensions = geo_tag.shape
    grid.origin = (0, 0, 0)
    grid.spacing = (1, 1, 1)
    grid.point_data["geo"] = geo_tag.flatten(order="F")
    #grid.plot(show_edges=False, opacity="linear")
    grid.plot(volume=True, opacity=[0.0, 0.9], shade=False)

if False:
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
