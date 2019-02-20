import numpy as np
import matplotlib.pyplot as plt
from random import random as rand

def cirular_shift(geo, dr):
    M = np.roll(geo, dr[0], axis=0)
    M = np.roll(M, dr[1], axis=1)
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

    for y in np.arange(0, geo.shape[0]):
        for x in np.arange(0, geo.shape[1]):
            f.write(str(int(geo[y, x])))
        f.write("\n")
    ## <end>
    f.write("<end>\n")
    f.close()

nx = 400
ny = 200

geo = np.zeros([ny, nx])

geo_one_sphere = np.zeros(geo.shape)
YY, XX = np.meshgrid(np.linspace(0, ny, ny), np.linspace(0, nx, nx), indexing='ij')

xm = nx // 2
ym = ny // 2
XX2 = (XX - xm)**2
YY2 = (YY - ym)**2
r2 = 14**2

geo_one_sphere[XX2 + YY2 < r2] = 1
cirular_shift(geo_one_sphere, [-ym, -xm])

for x in np.arange(0, nx-1, 40):
    row_num = 0
    for y in np.arange(0, ny-1, 40):
        vec_shift = [y, x]
        vec_shift[1] += (20 * (row_num % 2))
        for d in np.arange(len(vec_shift)):
            vec_shift[d] += 15*rand()
            vec_shift[d] = int(np.round(vec_shift[d]))
        geo_shift = cirular_shift(geo_one_sphere, vec_shift)
        geo = np.maximum(geo, geo_shift)
        row_num += 1


write_geo_file("test.dat", geo)

print(sum(geo.flatten())/geo.size)

plt.figure(1)
plt.pcolormesh(geo)
plt.gca().set_aspect('equal')
plt.show()