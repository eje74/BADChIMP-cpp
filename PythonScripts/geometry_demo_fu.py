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
# # - top and bottom plate
# geo[:,0] = 0
# geo[:,-1] = 0
# # - lower left corner
# geo[:10, :11] = 0
# # - upper left corner
# geo[:10, 21:] = 0
# # - lower right corner
# geo[30:, :11] = 0
# # - upper right corner
# geo[30:, 21:] = 0

# path to your badchimp folder
path_badchimp = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/"
# generate geometry input file(s)
vtk = vtklb(geo, "D2Q9", "x", "tmp", path_badchimp + "input/mpi/") 

# Set initial density
# - get the Cartesian coordinates
X, Y = np.mgrid[0:40, 0:32]
# - set rho
rho = np.zeros(sytems_size)

ii = int(sytems_size[0]/3)
jj = int(2*sytems_size[0]/5)
rho[:ii,:] = 1
vtk.append_data_set("init_rho_0", rho)
rho[:] = 0
rho[ii:jj,:] = 1
vtk.append_data_set("init_rho_1", rho)
rho[:] = 0
rho[jj:,:] = 1
vtk.append_data_set("init_rho_2", rho)

plt.figure(1)
rho[geo==0] = np.nan
ax = plt.pcolormesh(rho.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()
#plt.savefig(path_badchimp + "std_init_rho.pdf")
plt.show()

