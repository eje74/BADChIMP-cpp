#!/usr/bin/env python3
from operator import ge
from this import d
from vtklb import vtklb
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
#                S Y S T E M   G E O M E T R Y
# -------------------------------------------------------------------
# ------------------------------------------------------------------- System size: 
system_size = (81, 72) # nx, ny
geo = np.ones(system_size, dtype=int)
geo[:,36:] = 2
# ------------------------------------------------------------------- Add solid nodes
geo[:,0]  = 0
geo[:,-1] = 0
geo[-1,:] = 0
geo[:50,31:41] = 0

# ------------------------------------------------------------------- Set pressure boundary
pressure_boundary = np.zeros(system_size, dtype=int)
pressure_boundary[0, 1:31] = 1
pressure_boundary[0, 41:-1] = 2
#pressure_boundary[-1, 1:15] = 3

# pressure_boundary[0, :15] = 1
# pressure_boundary[0, 15:] = 2

# ------------------------------------------------------------------- Badchimp folder path
path_badchimp = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/"
# ------------------------------------------------------------------- Generate geometry input file(s)
vtk = vtklb(geo, "D2Q9", "", "tmp", path_badchimp + "input/mpi/") 
# ------------------------------------------------------------------- Set initial rho
rho = np.ones(system_size)
vtk.append_data_set("init_rho", rho)

# -------------------------------------------------------------------
#                       B O U N D A R Y
# -------------------------------------------------------------------
geo_padded = -np.ones( tuple(x+2 for x in system_size), dtype=int )
geo_padded[1:-1,1:-1] = geo
geo_padded[geo_padded>0] = 1
solid_wall_ind = np.where(geo_padded == 0)
# ------------------------------------------------------------------- Find fluid boundary
basis_vec = np.zeros((8,2), dtype=int)
basis_vec[0,:] = np.array([1, 0])
basis_vec[1,:] = np.array([0, 1])
basis_vec[2,:] = np.array([-1, 0])
basis_vec[3,:] = np.array([0, -1])
basis_vec[4,:] = np.array([1, 1])
basis_vec[5,:] = np.array([-1, 1])
basis_vec[6,:] = np.array([-1,-1])
basis_vec[7,:] = np.array([1, -1])
for basis_vec_element in basis_vec:
    wall_neig = tuple( x + bv for x, bv in zip(basis_vec_element, solid_wall_ind))
    geo_padded[wall_neig] = (geo_padded[wall_neig]> 0)*(2 - geo_padded[wall_neig]) + geo_padded[wall_neig]

pressure_bnd_ind = np.where(pressure_boundary>0)
pressure_bnd_ind = tuple(x+1 for x in pressure_bnd_ind)

geo_padded[pressure_bnd_ind] = 2

# ------------------------------------------------------------------- Find normal vectors
boundary_x_pos = np.zeros(system_size, dtype=int)
boundary_y_pos = np.zeros(system_size, dtype=int)
fluid_wall_ind = np.where(geo_padded == 2)
dist_wall = 10*np.ones(len(fluid_wall_ind[0]), dtype=int)
for count, point in enumerate(np.array(fluid_wall_ind).transpose()):
    for bv in basis_vec:
        if (geo_padded[tuple(point + bv)] == 0) | (geo_padded[tuple(point + bv)] == -1):
            dw = np.sum(np.abs(bv))
            if dw < dist_wall[count]:
                dist_wall[count] = dw
                boundary_x_pos[tuple(point-1)] = 0
                boundary_y_pos[tuple(point-1)] = 0
            if dw == dist_wall[count]:
                boundary_x_pos[tuple(point-1)] -= bv[0]
                boundary_y_pos[tuple(point-1)] -= bv[1]

# -------------------------------------------------------------------
#             W R I T E   V T K L B   A T T R I B U T E S
# -------------------------------------------------------------------
# ------------------------------------------------------------------- Set pressure boundary
#pressure_boundary = np.zeros(system_size, dtype=int)
#pressure_boundary[0, 1:-1] = 1
#pressure_boundary[-1, 1:-1] = 2
#pressure_boundary[-1, 1:15] = 3

# pressure_boundary[0, :15] = 1
# pressure_boundary[0, 15:] = 2
vtk.append_data_set("pressure_boundary", pressure_boundary)
# ------------------------------------------------------------------- Set boundary intersection
normal = np.sqrt(boundary_x_pos**2 + boundary_y_pos**2) 
normal[normal == 0] = 1
normal_vec_x = boundary_x_pos/normal
normal_vec_y = boundary_y_pos/normal
# normal_vec_x[0, 20] = 0
# normal_vec_y[0, 20] = 1
# normal_vec_x[0, -21] = 0
# normal_vec_y[0, -21] = -1
# normal_vec_x[-1, 1] = 0
# normal_vec_y[-1, 1] = 1
# normal_vec_x[-1, -2] = 0
# normal_vec_y[-1, -2] = -1
normal_vec_x[0, 1] = 0
normal_vec_y[0, 1] = 1

normal_vec_x[0, 30] = 0
normal_vec_y[0, 30] = -1

normal_vec_x[0, 41] = 0
normal_vec_y[0, 41] = 1

normal_vec_x[0, 70] = 0
normal_vec_y[0, 70] = -1


vtk.append_data_set("normal_x", normal_vec_x)
vtk.append_data_set("normal_y", normal_vec_y)
# ------------------------------------------------------------------- Set relative fraction
boundary_q = np.zeros(system_size)
boundary_q[tuple(i - 1 for i in fluid_wall_ind)] = 0.5
vtk.append_data_set("q", boundary_q)
# ------------------------------------------------------------------- Set position of neighbor
vtk.append_data_set("neighbor_x", boundary_x_pos)
vtk.append_data_set("neighbor_y", boundary_y_pos)

# -------------------------------------------------------------------
#                   P L O T   V A R I A B L E S
# -------------------------------------------------------------------
plt.figure()
plt.title("geometry")
ax = plt.pcolormesh(geo.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("rho")
ax = plt.pcolormesh(rho.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("boundary number")
ax = plt.pcolormesh(pressure_boundary.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("normal_x pos")
ax = plt.pcolormesh((boundary_x_pos).transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("normal_y pos")
ax = plt.pcolormesh((boundary_y_pos).transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("q")
ax = plt.pcolormesh(boundary_q.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()


plt.show()

