import numpy as np
import os
import argparse
import time
from generate_geometry_mpi import generate_geometry_mpi
from generate_geometry_mpi import write_input_hpc

p = argparse.ArgumentParser()
p.add_argument("--outputpath", type=str, required=True)
p.add_argument("--solidfilename", type=str, required=True)
p.add_argument("--fluidfilename", type=str, required=True)
p.add_argument("--outfilename", type=str, required=True)
args = p.parse_args()

outputpath = args.outputpath
solidfilename = args.solidfilename
fluidfilename = args.fluidfilename
outfilename = args.outfilename

inputpath = r"./"
datapath = r"/cluster/home/esje/github/BADChIMP-cpp/input/data/"


# MPI setup
# Initial regular composition of geo-array
num_proc = (4,)*3


# ------------------------------------------------------------------------ generate geometry
pore = np.load(
    datapath + solidfilename +  r".npy")

fluid = np.load(
    datapath + fluidfilename +  r".npy")

geo = np.ones(
    pore.shape, 
    dtype=np.int32)

# -------------------------------------------------------------------- Wetting phase
geo[fluid<=0] = 0
geo[pore>0] = 0

basefilename = outfilename + "_W"
# nproc is the actually number of processors used
#  sbatch reads the from the ntasks_datafile
nproc = generate_geometry_mpi(
    geo, 
    num_proc, 
    inputpath + r"input/mpi/", 
    vtklbfilename=basefilename
    )

with open(outputpath + "ntasks" + basefilename + ".txt", "w") as f:
    f.write(str(nproc) + "\n")

hpc_write_args = (
    inputpath,
    20000,
    1000,
    0.8,
    1e-6,
    outputpath
)

write_input_hpc(
    *hpc_write_args,    
    basefilename
)

# -------------------------------------------------------------------- Non-wetting phase
geo[:] = 1
geo[fluid>0] = 0
geo[pore>0] = 0

basefilename = outfilename + "_NW"
# nproc is the actually number of processors used
#  sbatch reads the from the ntasks_datafile
nproc = generate_geometry_mpi(
    geo, 
    num_proc, 
    inputpath + r"input/mpi/", 
    vtklbfilename=basefilename
    )

with open(outputpath + "ntasks" + basefilename + ".txt", "w") as f:
    f.write(str(nproc) + "\n")

write_input_hpc(
    *hpc_write_args,    
    basefilename
)
