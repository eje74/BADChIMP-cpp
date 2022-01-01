#!/usr/bin/env python3
from vtklb import vtklb
import numpy as np
import matplotlib.pyplot as plt

# system size: nx  ny
sytems_size = (40, 32)
# set the size of the geometri
geo = np.ones(sytems_size, dtype=int)
# Add solids to
# - top and bottom plate
geo[:,0] = 0
geo[:,-1] = 0

# path to your badchimp folder
path_badchimp = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/"
# generate geometry input file(s)
vtk = vtklb(geo, "D2Q9", "", "tmp", path_badchimp + "input/mpi/") 

# Set pressure boundary
pressure_boundary = np.zeros(sytems_size, dtype=int)
pressure_boundary[0, 1:-1] = 1
pressure_boundary[-1, 1:-1] = 2
vtk.append_data_set("pressure_boundary", pressure_boundary)

# Set boundary intersection
boundary_x_pos = np.zeros(sytems_size, dtype=int)
boundary_x_pos[:,1] = 0
boundary_x_pos[:,-2] = 0
vtk.append_data_set("x_pos", boundary_x_pos)


# Set boundary intersection
boundary_y_pos = np.zeros(sytems_size, dtype=int)
boundary_y_pos[:,1] = -1
boundary_y_pos[:,-2] = 1
vtk.append_data_set("y_pos", boundary_y_pos)

# Set relative fraction
boundary_q = np.zeros(sytems_size)
boundary_q[:,1] = 0.5
boundary_q[:,-2] = 0.5
vtk.append_data_set("q", boundary_q)

plt.figure(1)
plt.title("geometry")
ax = plt.pcolormesh(geo.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure(2)
plt.title("boundary number")
ax = plt.pcolormesh(pressure_boundary.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure(3)
plt.title("x pos")
ax = plt.pcolormesh(boundary_x_pos.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure(4)
plt.title("y pos")
ax = plt.pcolormesh(boundary_y_pos.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure(5)
plt.title("q")
ax = plt.pcolormesh(boundary_q.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()


plt.show()

