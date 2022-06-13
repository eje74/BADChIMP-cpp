#!/usr/bin/env python3

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
#geo[0,:]  = 0
geo[:50,31:41] = 0

# ------------------------------------------------------------------- Badchimp folder path
path_badchimp = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/"
# ------------------------------------------------------------------- Generate geometry input file(s)
vtk = vtklb(geo, "D2Q9", "", "tmp", path_badchimp + "input/mpi/") 

# -------------------------------------------------------------------
#             W R I T E   V T K L B   A T T R I B U T E S
# -------------------------------------------------------------------
# ------------------------------------------------------------------- Set pressure boundary
pressure_boundary = np.zeros(system_size, dtype=int)
pressure_boundary[0, 0:32] = 1
pressure_boundary[0, 40:] = 2

vtk.append_data_set("pressure_boundary", pressure_boundary)

# ------------------------------------------------------------------- Signed distance function
xx, yy = np.meshgrid(np.array(np.arange(81)), np.array(np.arange(72)), indexing='ij')

# signed_dist = np.abs(xx-0.5)
# m_add = (-xx + np.max(xx)) - 0.5
# m_add = np.abs( m_add)
# signed_dist = np.minimum(m_add, signed_dist)

signed_dist =  np.abs((-xx + np.max(xx)) - 0.5)


m_add = np.abs(yy-0.5)
signed_dist = np.minimum(m_add, signed_dist)

m_add = (-yy + np.max(yy)) - 0.5
m_add = np.abs( m_add)
signed_dist = np.minimum(m_add, signed_dist)

m_add = np.abs(yy-30.5)
signed_dist[1:50,:] = np.minimum(m_add[1:50,:], signed_dist[1:50,:])
m_add = np.abs(yy-40.5)
signed_dist[1:50,:] = np.minimum(m_add[1:50,:], signed_dist[1:50,:])

m_add = np.abs(xx-49.5)
signed_dist[:, 31:41] = np.minimum(m_add[:, 31:41], signed_dist[:, 31:41])

# Set courners
x0, y0 = 0.5, 0.5
m_add = np.abs(np.sqrt( (xx-x0)**2 + (yy-y0)**2 ))
signed_dist = np.minimum(m_add, signed_dist)

x0, y0 = 79.5, 0.5
m_add = np.abs(np.sqrt( (xx-x0)**2 + (yy-y0)**2 ))
signed_dist = np.minimum(m_add, signed_dist)

x0, y0 = 79.5, 70.5
m_add = np.abs(np.sqrt( (xx-x0)**2 + (yy-y0)**2 ))
signed_dist = np.minimum(m_add, signed_dist)

x0, y0 = 0.5, 70.5
m_add = np.abs(np.sqrt( (xx-x0)**2 + (yy-y0)**2 ))
signed_dist = np.minimum(m_add, signed_dist)

x0, y0 = 0.5, 30.5
m_add = np.abs(np.sqrt( (xx-x0)**2 + (yy-y0)**2 ))
signed_dist = np.minimum(m_add, signed_dist)

x0, y0 = 0.5, 40.5
m_add = np.abs(np.sqrt( (xx-x0)**2 + (yy-y0)**2 ))
signed_dist = np.minimum(m_add, signed_dist)

x0, y0 = 49.5, 30.5
m_add = np.abs(np.sqrt( (xx-x0)**2 + (yy-y0)**2 ))
signed_dist = np.minimum(m_add, signed_dist)

x0, y0 = 49.5, 40.5
m_add = np.abs(np.sqrt( (xx-x0)**2 + (yy-y0)**2 ))
signed_dist = np.minimum(m_add, signed_dist)

signed_dist[geo==0] = -signed_dist[geo==0]


vtk.append_data_set("signed_distance", signed_dist)

# -------------------------------------------------------------------
#      S E T   T H E   I N L E T / O U T L E T   G E O M E T R Y
# -------------------------------------------------------------------
vtk.begin_point_data_subset_section()
# ------------------------------------------------------------------- Set mask
mask = pressure_boundary > 0
# ------------------------------------------------------------------- Set normal
normal = np.zeros((2,) + mask.shape)
normal[0, pressure_boundary==1] = 1
normal[0, pressure_boundary==2] = 1

char_index = ["_x", "_y", "_z"]
for i in np.arange(2):
    vtk.append_data_subset("boundary_normal" + char_index[i], normal[i, ...], mask)
# ------------------------------------------------------------------- Set distance from boundary
s_dist = np.zeros(mask.shape)
s_dist[pressure_boundary==1] = 0.0
s_dist[pressure_boundary==2] = 0
vtk.append_data_subset("boundary_distance", s_dist, mask)

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
plt.title("boundary number")
ax = plt.pcolormesh(pressure_boundary.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()

plt.figure()
plt.title("signed dist")
ax = plt.pcolormesh(signed_dist.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()





plt.show()

