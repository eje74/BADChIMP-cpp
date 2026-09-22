import numpy as np
import os
import argparse
import time
from generate_geometry_mpi import generate_geometry_mpi
from generate_geometry_mpi import write_input_hpc

p = argparse.ArgumentParser()

p.add_argument("--outputpath", type=str, required=True)
p.add_argument("--datafilename", type=str, required=True)
p.add_argument("--outfilename", type=str, required=True)
args = p.parse_args()

outputpath = args.outputpath
datafilename = args.datafilename
outfilename = args.outfilename

#username = "olau" #"esje"
inputpath = r"./"
#datapath = r"/cluster/home/"+username+"/github/BADChIMP-cpp/input/data/"
datapath = inputpath+"LSdata/"

# filenamedata = r"GH_PoreSolid_400x400x400_SDF_PD"
# filename = datapath + filenamedata +  r".npy"
filename = datapath + datafilename +  r".npy"

# ------------------------------------------------------------------------ generate geometry
pore = np.load(filename)
geo = np.ones(pore.shape, dtype=np.int32)
geo[pore>0] = 0

num_proc = (4,)*3

nproc = generate_geometry_mpi(
    geo, 
    num_proc, 
    inputpath + r"input/mpi/", 
    vtklbfilename=outfilename
    )

with open(outputpath + "ntasks" + outfilename + ".txt", "w") as f:
    f.write(str(nproc) + "\n")

write_input_hpc(
    inputpath,
    20000,
    500,
    0.8,
    1e-6,
    outputpath,
    outfilename
)
