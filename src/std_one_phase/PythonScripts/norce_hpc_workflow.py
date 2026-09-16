import numpy as np
import os
import argparse
import time
from generate_geometry_mpi import generate_geometry_mpi
from generate_geometry_mpi import write_input_file

p = argparse.ArgumentParser()
p.add_argument("--outputpath", type=str, required=True)
args = p.parse_args()

outputpath = args.outputpath

print("outputpath = ", outputpath)


path = "/cluster/home/esje/"
lbpath = path + r"github/BADChIMP-cpp/"
filename = lbpath + r"input/data/GH_PoreSolid_400x400x400_SDF_PD.npy"


# pore = np.load(filename)

# # ------------------------------------------------------------------------ generate geometry
# geo = np.ones(pore.shape, dtype=np.int32)
# geo[pore>0] = 0

# num_proc = (4,)*3

# nproc = generate_geometry_mpi(geo, num_proc, lbpath + r"input/mpi/")

# write_input_file(
#     lbpath,
#     5000,
#     500,
#     0.8,
#     1e-6,
#     "TestSolid2"
# )