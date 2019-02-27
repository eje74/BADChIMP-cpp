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

nx = 120
ny = 140
nz = 120

geo = np.zeros([nz, ny, nx])

geo_one_sphere = np.zeros(geo.shape)
ZZ, YY, XX = np.meshgrid(np.linspace(0, nz, nz), np.linspace(0, ny, ny), np.linspace(0, nx, nx), indexing='ij')

xm = nx // 2
ym = ny // 2
zm = nz // 2
XX2 = (XX - xm)**2
YY2 = (YY - ym)**2
ZZ2 = (ZZ - zm)**2
r2 = 14**2

geo_one_sphere[XX2 + YY2 + ZZ2 < r2] = 1

#translation operation that shifts first sphere to origin (0,0)
#geo_shift=cirular_shift(geo_one_sphere, [-ym, -xm])

#Generates spheres and consequentially shifts them
for x in np.arange(0, nx-1, 40):
    row_num = 0
    for y in np.arange(0, ny-1, 35):
        level_num = 0
        for z in np.arange(0, nz-1, 35):
            vec_shift = [z-zm, y-ym, x-xm]
            vec_shift[1] += (20 * (level_num % 2))
            vec_shift[2] += (20 * (row_num % 2))#(20 * (row_num % 2))
            for d in np.arange(len(vec_shift)):
                vec_shift[d] += 0#15*rand()
                vec_shift[d] = int(np.round(vec_shift[d]))
            geo_shift = cirular_shift(geo_one_sphere, vec_shift)
            geo = np.maximum(geo, geo_shift)
            level_num += 1
        row_num += 1
            #plt.figure(10*x+y)
            #plt.pcolormesh(geo)
            #plt.gca().set_aspect('equal')
            #plt.show()


write_geo_file("test2.dat", geo)

print(sum(geo.flatten())/geo.size)

plt.figure(1)
plt.pcolormesh(geo[0,:,:])
plt.gca().set_aspect('equal')

plt.figure(2)
plt.pcolormesh(geo[35,:,:])
plt.gca().set_aspect('equal')

plt.figure(3)
plt.pcolormesh(geo[:,:,0])
plt.gca().set_aspect('equal')
plt.show()