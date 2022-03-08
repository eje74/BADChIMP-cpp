#!/usr/bin/env python3
from operator import ge
import sys
#from this import d
from vtklb import vtklb
import numpy as np
import matplotlib.pyplot as plt

def expand_in_z(vals, system_size_xyz):
    tmp = np.zeros(system_size_xyz, dtype=vals.dtype)
    for nz in np.arange(system_size_xyz[-1]):
        tmp[:,:,nz] = vals
    return tmp

# -------------------------------------------------------------------
#                S Y S T E M   G E O M E T R Y
# -------------------------------------------------------------------
# ------------------------------------------------------------------- System size xy: 
N = 150
system_size = (N, N) # nx, ny
system_size_xyz = (N, N, 240) # nx, ny, nz
geo = np.ones(system_size, dtype=int)
# ------------------------------------------------------------------- Add solid nodes xy
geo[75:,:75] = 2
geo[:75,75:] = 3
geo[75:,75:] = 4
# ------------------------------------------------------------------- Add solid nodes xy
X = np.linspace(-1.1, 1.1, N)
print([X[0], X[-1]])
xx, yy = np.meshgrid(X, X, indexing='ij')

geo[xx**2 + yy**2 >= 1] = 0
# ------------------------------------------------------------------- Set pressure boundary
pressure_boundary = np.zeros(system_size, dtype=int)
# ------------------------------------------------------------------- Badchimp folder path
path_badchimp = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/"
# ------------------------------------------------------------------- Generate geometry input file(s)
tmp = expand_in_z(geo, system_size_xyz)
geo[xx**2 + yy**2 >= 0.8] = 0
for dz in np.arange(20):
    tmp[:,:, 110+dz] = geo 

for n in np.arange(1,3):
    tmp[:,:,60*n:60*(n+1)] += (tmp[:,:,60*n:60*(n+1)] > 0)*(4*n) 
tmp[:,:,60*3:] += (tmp[:,:,60*3:] > 0)*4*3

vtk = vtklb(tmp, "D3Q19", "z", "tmp", path_badchimp + "input/mpi/") 

# ------------------------------------------------------------------- Plot geometry xz
plt.figure()
plt.title("geometry x z")
plt.pcolormesh(tmp[:,75,:])


# ------------------------------------------------------------------- Set initial rho
rho = np.ones(system_size)
tmp = expand_in_z(rho, system_size_xyz)
vtk.append_data_set("init_rho", tmp)

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
# ------------------------------------------------------------------- Relative fraction
boundary_q = np.zeros(system_size)
# ------------------------------------------------------------------- Normal vectors
normal_vec_x = np.zeros(system_size)
normal_vec_y = np.zeros(system_size)
# ------------------------------------------------------------------- Find bulk fluid neighbors
boundary_x_pos = np.zeros(system_size, dtype=int)
boundary_y_pos = np.zeros(system_size, dtype=int)
fluid_wall_ind = np.where(geo_padded == 2)
dist_wall = 10*np.ones(len(fluid_wall_ind[0]), dtype=int)
for count, point in enumerate(np.array(fluid_wall_ind).transpose()):
    for bv in basis_vec:
        if (geo_padded[tuple(point + bv)] == 0):
            dw = np.sum(np.abs(bv))
            if dw < dist_wall[count]:
                dist_wall[count] = dw
                boundary_x_pos[tuple(point-1)] = 0
                boundary_y_pos[tuple(point-1)] = 0
            if dw == dist_wall[count]:
                boundary_x_pos[tuple(point-1)] -= bv[0]
                boundary_y_pos[tuple(point-1)] -= bv[1]
                # # find q
                # x = xx[tuple(point-1)]
                # y = yy[tuple(point-1)]
                # vx = xx[tuple(point+bv-1)] - x
                # vy = yy[tuple(point+bv-1)] - y
                # c = x**2 + y**2 -1
                # b = 2*(x*vx + y*vy)
                # a = vx**2 + vy**2
                # q = (-b + np.sqrt(b**2 - 4*a*c))/(2*a) 
                # boundary_q[tuple(point-1)] = q
                # n_norm = (x + q*vx)**2 + (y + q*vy)**2
                # normal_vec_x[tuple(point-1)] = -(x + q*vx)/n_norm 
                # normal_vec_y[tuple(point-1)] = -(y + q*vy)/n_norm 
    # find q
    bv = np.array([-boundary_x_pos[tuple(point-1)], -boundary_y_pos[tuple(point-1)]], dtype=int)
    x = xx[tuple(point-1)]
    y = yy[tuple(point-1)]
    vx = xx[tuple(point+bv-1)] - x
    vy = yy[tuple(point+bv-1)] - y
    c = x**2 + y**2 -1
    b = 2*(x*vx + y*vy)
    a = vx**2 + vy**2
    q = (-b + np.sqrt(b**2 - 4*a*c))/(2*a) 
    boundary_q[tuple(point-1)] = q
    n_norm = (x + q*vx)**2 + (y + q*vy)**2
    normal_vec_x[tuple(point-1)] = -(x + q*vx)/n_norm 
    normal_vec_y[tuple(point-1)] = -(y + q*vy)/n_norm 

# -------------------------------------------------------------------
#             W R I T E   V T K L B   A T T R I B U T E S
# -------------------------------------------------------------------
# ------------------------------------------------------------------- Set pressure boundary
tmp = expand_in_z(pressure_boundary, system_size_xyz)
vtk.append_data_set("pressure_boundary", tmp)
# -------------------------------------------------------------------S et boundary intersection

tmp = expand_in_z(normal_vec_x, system_size_xyz)
vtk.append_data_set("normal_x", tmp)
tmp = expand_in_z(normal_vec_y, system_size_xyz)
vtk.append_data_set("normal_y", tmp)
tmp[:] = 0 
vtk.append_data_set("normal_z", tmp)
# ------------------------------------------------------------------- Set relative fraction
tmp = expand_in_z(boundary_q, system_size_xyz)
vtk.append_data_set("q", tmp)
# ------------------------------------------------------------------- Set position of neighbor
tmp = expand_in_z(boundary_x_pos, system_size_xyz)
vtk.append_data_set("neighbor_x", tmp)
tmp = expand_in_z(boundary_y_pos, system_size_xyz)
vtk.append_data_set("neighbor_y", tmp)
tmp[:] = 0
vtk.append_data_set("neighbor_z", tmp)

# -------------------------------------------------------------------
#                   P L O T   V A R I A B L E S
# -------------------------------------------------------------------
plt.figure()
plt.title("geometry xy")
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
plt.title("boundary_x pos")
ax = plt.pcolormesh((boundary_x_pos).transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("boundary_y pos")
ax = plt.pcolormesh((boundary_y_pos).transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("normal_vec_x")
ax = plt.pcolormesh((normal_vec_x).transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("normal_vec_y")
ax = plt.pcolormesh((normal_vec_y).transpose())
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

print(f"N = {2/(X[1]-X[0])}")