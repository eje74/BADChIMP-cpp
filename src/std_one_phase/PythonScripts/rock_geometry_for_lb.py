import numpy as np
import pyvista as pv
from matplotlib import pyplot as plt
from vtklb import vtklb

# ======================================================================== Data format
# Solid geometry:
#   element > 0: solid
# Fluid geometry:
#   element > 0: wetting
   
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
