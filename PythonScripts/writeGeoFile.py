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

nx = 30
ny = 12
nz = 3

geo = np.zeros([nz, ny, nx])

geo[:,0,:] = 1
geo[:,-1,:] = 1

write_dir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/PythonScripts/"  # Home

write_geo_file(write_dir + "walls.dat", geo)

print(sum(geo.flatten())/geo.size)

plt.figure(1)
plt.pcolormesh(geo[0,:,:])
plt.gca().set_aspect('equal')
plt.colorbar()

plt.figure(2)
plt.pcolormesh(geo[:,0,:])
plt.gca().set_aspect('equal')
plt.colorbar()

plt.figure(3)
plt.pcolormesh(geo[:,:,0])
plt.gca().set_aspect('equal')
plt.colorbar()
plt.show()