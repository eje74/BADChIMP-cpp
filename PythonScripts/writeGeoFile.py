import numpy as np
import matplotlib.pyplot as plt
from random import random as rand

def cirular_shift(geo, dr):
    M = np.roll(geo, dr[0], axis=0)
    for dim in range(1, len(dr)):
        M = np.roll(M, dr[dim], axis=dim)
    return M


def write_geo_file(filename, geo, res=1.0):
    f = open(filename, "w")
    ## Write dimensions x y
    f.write("dim")
    for dim in np.arange(geo.ndim, 0, -1):
        f.write(" " + str(geo.shape[dim-1]))
    f.write("\n")
    ## Write resolution float
    f.write("res " + str(res) + "\n")
    ## <geo char>
    f.write("<geo char>\n")
    ## write geometry

    for z in np.arange(0, geo.shape[0]):
        for y in np.arange(0, geo.shape[1]):
            for x in np.arange(0, geo.shape[2]):
                f.write(str(int(geo[z, y, x])))
            f.write("\n")
    ## <end>
    f.write("<end>\n")
    f.close()



# System size
nx = 5
ny = 7
nz = 1

# SETUP GEOMETRY with rank (0: SOLID, 1:RANK0, 2:RANK1, ...)
geo = np.ones([nz, ny, nx])

# Position boundaries
geo[:,0,:] = 0
geo[:,-1,:] = 0

# Partition the process
geo[:, :4, :] = 2*geo[:, :4, :]

# Write file
write_dir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/PythonScripts/"  # Home
write_geo_file(write_dir + "walls.dat", geo)
