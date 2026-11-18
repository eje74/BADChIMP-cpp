#!/usr/bin/env python3
from vtklb import vtklb
import numpy as np
import matplotlib.pyplot as plt

# system size: nx  ny
sytems_size = (220 + 2, 41 + 2)
# sytems_size = (2, 4)
# set the size of the geometri
geo = np.ones(sytems_size, dtype=int)
# set boundaries
geo[:, 0] = 0
geo[:, -1] = 0

geo[100:120, 0:10] = 0
geo[100:120, -10:] = 0



# From physical to lb units
Re = 20
dx = 0.01  # m

D = 0.1
nu = 1e-3 # m^2 / s

Dlb = D/dx
Ulb = 0.01

tau = 3*Ulb*Dlb/Re + 0.5 

Fd = 0.5*(0.0555 + 0.0538)
Fw = 0.5*(0.000127 + 0.000107)

Cd = 2*Fd / (Ulb**2 * Dlb )

Cw = 2*Fw / (Ulb**2 * Dlb )

print("Dlb = ", Dlb)
print("tau = ", tau)
print("Cd = ", Cd)
print("Cw (x100) = ", Cw * 100)

print(geo.size)

# partition the system
# path to your badchimp folder
path_badchimp = "/home/AD.NORCERESEARCH.NO/esje/GitHub/BADChIMP-cpp/"
# generate geometry input file(s)
vtk = vtklb(geo, "D2Q9", "xy", "tmp", path_badchimp + "input/mpi/") 

# Set initial density
rho = np.ones(sytems_size, dtype=np.float32)
# Add rho to the geometry files
vtk.append_data_set("init_rho", rho)

# Boundary indicator
boundary = np.zeros(sytems_size, dtype=np.int32)
boundary[0, :] = 1
boundary[-1, :] = 2
boundary[:, 0] = 3
boundary[:, -1] = 3
# Add the boundary to the geometry file
vtk.append_data_set("boundary", boundary)

# Set particle indicator field
X, Y = np.mgrid[0:sytems_size[0], 0:sytems_size[1]]

# Set boundary
# X = 0.5*X
# Y = 0.5*Y


# d = np.abs( Y - 0.5*(sytems_size[1]) ) - 3

# Set particle
print(np.max(X), np.min(X))
d = np.sqrt( (X - 20.25)**2 + (Y - 20.25)**2 ) - Dlb/2


eps = 1.0
eps *= 0.5

gamma = 0.5 * (1 - np.sin(0.5*3.14159 * d / eps))
gamma[d>eps] = 0
gamma[d<-eps] = 1

delta = (0.5/eps)*(1 + np.cos(3.14159*d/eps))
delta[np.abs(d)>eps]  = 0

# normal vector
Nx = np.zeros_like(d)
Nx[1:-1, :] = 0.5 * (d[2:, :] - d[:-2, :])
Ny = np.zeros_like(d)
Ny[:, 1:-1] = 0.5 * (d[:, 2:] - d[:, :-2])

tmp = np.sqrt(Nx**2 + Ny**2)
Nx[tmp == 0.0] = 1.
tmp[tmp == 0.0]  = 1.
Nx /= tmp
Ny /= tmp

#gamma = np.zeros(sytems_size, dtype=np.float32)
# Add gamma to the geometry files
vtk.append_data_set("gamma", gamma)
vtk.append_data_set("delta", delta)
vtk.append_data_set("Nx", Nx)
vtk.append_data_set("Ny", Ny)

# plt.figure()
# plt.title("Nx")
# Nx[geo==0] = np.nan
# ax = plt.pcolormesh(Nx.transpose())
# plt.axis('equal')
# plt.axis('off')
# plt.colorbar()

# plt.figure()
# plt.title("Ny")
# Ny[geo==0] = np.nan
# ax = plt.pcolormesh(Ny.transpose())
# plt.axis('equal')
# plt.axis('off')
# plt.colorbar()


# plt.figure()
# plt.title("rho")
# rho[geo==0] = np.nan
# ax = plt.pcolormesh(rho.transpose())
# plt.axis('equal')
# plt.axis('off')
# plt.colorbar()

plt.figure()
plt.title("gamma")
gamma[geo==0] = np.nan
ax = plt.pcolormesh(gamma.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

# plt.figure()
# plt.title("delta")
# delta[geo==0] = np.nan
# ax = plt.pcolormesh(delta.transpose())
# plt.axis('equal')
# plt.axis('off')
# plt.colorbar()


plt.figure()
plt.title("geo")
ax = plt.pcolormesh(geo.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()


plt.show()

