import numpy as np
import os
import time
from generate_geometry_mpi import generate_geometry_mpi
from generate_geometry_mpi import write_input_file

lbpath = r"/home/AD.NORCERESEARCH.NO/esje/GitHub/BADChIMP-cpp/"

filename = r"/home/AD.NORCERESEARCH.NO/esje/OneDrive/NORCE/CSSR/RelPerm LB LS/Gildehauser/PD_Scontrol_Tow20/GH_PoreSolid_400x400x400_SDF_PD.npy"

pore = np.load(filename)

# ------------------------------------------------------------------------ generate geometry
geo = np.ones(pore.shape, dtype=np.int32)
geo[pore>0] = 0

num_proc = (3,)*3

nproc = generate_geometry_mpi(geo, num_proc, lbpath + r"input/mpi/")

write_input_file(
    lbpath,
    500,
    50,
    0.8,
    1e-6,
    "TestSolid"
)