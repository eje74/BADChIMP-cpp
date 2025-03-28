#!/usr/bin/env python3
from vtklb import vtklb
import numpy as np
import matplotlib.pyplot as plt

# system size: nx  ny
sytems_size = (40, 32)
# set the size of the geometri
geo = np.ones(sytems_size, dtype=int)
# partition the system in two
geo[20:, :] = 2
# Add solids to
# - top and bottom plate
geo[:,0] = 0
geo[:,-1] = 0
# - lower left corner
geo[:10, :11] = 0
# - upper left corner
geo[:10, 21:] = 0
# - lower right corner
geo[30:, :11] = 0
# - upper right corner
geo[30:, 21:] = 0

# path to your badchimp folder
path_badchimp = "/home/AD.NORCERESEARCH.NO/esje/Tmp/"
# generate geometry input file(s)
vtk = vtklb(geo, "D2Q9", "x", "tmp", path_badchimp + "input/mpi/") 

# Set initial density
# - get the Cartesian coordinates
X, Y = np.mgrid[0:40, 0:32]
# - set rho
rho = 1.0 + 0.1*np.sin(2*np.pi*X/40)

# Add rho to the geometry files
vtk.append_data_set("init_rho", rho)

plt.figure()
plt.title("rho")
rho[geo==0] = np.nan
ax = plt.pcolormesh(rho.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("geo")
ax = plt.pcolormesh(geo.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()


plt.show()

