#!/usr/bin/env python3
import vtklb
import numpy as np
import matplotlib.pyplot as plt

# system size: nx  ny
sytems_size = (40, 32)
# set the size of the geometri
geo = np.ones(sytems_size)
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
path_badchimp = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/"
vtk = vtklb.vtklb(geo, "D2Q9", "x", "tmp", path_badchimp + "input/mpi/") 

plt.figure(1)
ax = plt.pcolormesh(geo.transpose())
plt.axis('equal')
plt.axis('off')
plt.savefig(path_badchimp + "std_geo.pdf")
plt.show()

