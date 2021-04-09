#!/usr/bin/env python3

from read_unstructured_grid import get_geo
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import vtklb
#import vtk_tools

directory="/cluster/home/janlv/BADChIMP-cpp/"
# Get the wing profile
#geo = get_geo(directory+"PythonScripts/liege1_2.vtk")
#data = vtk_tools.read_vtk(directory+"PythonScripts/liege1.vtk")
#file = open(directory+"PythonScripts/liege1_2.vtk", 'rb')
#geo1 = np.fromfile(file, dtype=np.float32)
#file.read(8)
file2 = "Porer-70kv2.modif.raw"  # void = 1, solid = 0
geo1 = np.fromfile(file2, dtype=np.int8)
geo2 = np.transpose(geo1.reshape((1422, 459, 455)))
geo3 = geo2[:100,:100,:100]
#geo3 = geo2[:100,:100,:1]
#geo3 = geo2[:50,:50,:1]

#geo3 = np.ones(geo3.shape, dtype=float)
#geo3[:,-1,:] = 0
#geo3[:,0,:] = 0
#geo3[:,:,-1] = 0
#geo3[:,:,0] = 0
#geo3[-1,:,:] = 0
#geo3[0,:,:] = 0


# Mirror geometry
geo4=np.zeros(geo3.shape, dtype=int)
geo4=geo3[::-1,:,:]
geo5=np.concatenate((geo3,geo4),axis=0)

geo3=geo5

# Distribute water evenly over solid surface
dist_transf = ndimage.distance_transform_edt(geo3)
dist_transf = np.array(dist_transf)
dist_transf[dist_transf==0]=1000
dist_ind = np.sort(dist_transf.flatten())
num_el = sum(dist_ind<1000)

# Water saturation
Sw=0.9
num_water = int(Sw*num_el)

#number of nodes less than cut
num_less_cut = sum(dist_ind<dist_ind[num_water])
S = np.zeros(geo3.shape, dtype=float)
S[dist_transf < dist_ind[num_water]] = 1.0

num_on_cut = sum(dist_ind == dist_ind[num_water])
# Local water concentration field 
S[dist_transf == dist_ind[num_water]] = (1.0*(num_water - num_less_cut))/num_on_cut

# BADChIMP require solid = 0, void >= 1
geo3=1-geo3  # solid = 1, void = 0
#geo3[:,:,:]=0

#geo3[:,-1,:] = 1
#geo3[:,0,:] = 1

#geo3[-1,:,:] = 1
#geo3[0,:,:] = 1






#s = (140, 92, 100)
#geo = data['geo']
#geo2 = np.ndarray(shape=(140, 92, 100), geo1)
#geo=np.array([v.replace(',', '') for v in geo2], dtype=np.float32)
geo=geo3

s = geo.shape
#geo_test = np.zeros( (s[0], s[1]*3, s[2]), dtype=int)
#for y in range(s[1]):
#    geo_test[:,3*y:(3*y+3),:] = geo[:,y,:].reshape(s[0], 1, s[2])
#geo = geo_test

#geo_large = np.zeros((geo.shape[0]+40, geo.shape[1]+6, geo.shape[2]+20), dtype=int)
#geo_large[20:(20+geo.shape[0]), :geo.shape[1], 10:(10 + geo.shape[2])] = geo
#geo_large[:,:5,:] = 0
#geo_large[: ,:2, :] = 0
#val = np.zeros((geo_large.shape[0]+100, geo_large.shape[1], geo_large.shape[2]), dtype=int)
#val = np.zeros((70, 200, 70), dtype=int)
#val[:,:96 ,:] = 1
#val[:,96:192,:] = 2
#val[:, 192:, :] = 3
#val[:geo_large.shape[0], :, :] = val[:geo_large.shape[0], :, :]*(geo_large == 0) 


#geo[:,-1,:] = 1
#geo[:,0,:] = 1
#geo[:,:,-1] = 1
#geo[:,:,0] = 1
#geo[-1,:,:] = 1
#geo[0,:,:] = 1

#val = np.zeros((50, 50, 50), dtype=int)
# MPI load partitioning
val = np.zeros((s[0], s[1], s[2]), dtype=int)
val[:int(s[0]/3.),: ,:] = 1
val[int(s[0]/3.):int(s[0]*2./3.),:,:] = 2
val[int(s[0]*2/3.):, :, :] = 3
val=val*(geo == 0)  # val = 0 for solid nodes

#Solid nodes
#val[:,0,:] = 0 #y = 0
#val[:,-1,:] = 0 #y = Y_max


#val[np.int(val.shape[0]*0.55):, np.int(val.shape[1]*0.75), :] = 0
#val[:np.int(val.shape[0]*0.45), np.int(val.shape[1]*0.75), :] = 0
#val[:, np.int(val.shape[1]*0.75), np.int(val.shape[2]*0.55):] = 0
#val[:, np.int(val.shape[1]*0.75), :np.int(val.shape[2]*0.45)] = 0

# Periodic in x,y,z
vtk = vtklb.vtklb(val, "D3Q19", "xyz", "tmp", directory+"input/mpi/") 

# Setup boundary marker
# Do we need this?
bnd = np.zeros(val.shape, dtype=int)
#x boundaries
#bnd[0, 1:-1, 1:-1] = 1 # Inlet x = 0
#bnd[-1, 1:-1, 1:-1] = 2 # Outlet x = -1
#y boundaries
bnd[:, 0, :] = 3 # Right hand boundary y = 0
bnd[:, -1, :] = 4 # Left hand boundary y = -1
#z boundaries 
bnd[:, 1:-1, 0] = 5 # Bottom boundary z = 0
bnd[:, 1:-1, -1] = 6 # Top boundary z = 0
vtk.append_data_set("boundary", bnd)

qSrc = np.zeros(val.shape, dtype=int)
##y boundaries
#qSrc[:, 1:3, :] = 1 # sink
#qSrc[:, -3:-2, :] = 2 # constant density
#qSrc=qSrc*(geo == 0)
#
vtk.append_data_set("source", qSrc)

#rho0 = np.ones(val.shape, dtype=float)
#phiInd = np.zeros(val.shape, dtype=float)
#phiDiff0 = np.zeros(val.shape, dtype=float)
#phiDiff1 = np.zeros(val.shape, dtype=float)
#phiDiff2 = np.zeros(val.shape, dtype=float)

#nx=rho0.shape[0]
#ny=rho0.shape[1]

#r0=0.25*ny
#x0=0.25*nx
#y0=0.5*ny

#r01=0.17*ny
#x01=0.75*nx
#y01=0.5*ny
#H=0.08
#sigma=1.e-1

# #for x in xrange(nx):
# #    for y in xrange(ny):
# for x in range(nx):
#     for y in range(ny):
#         rSq= (x-x0)*(x-x0) + (y-y0)*(y-y0)
#         r1Sq = (x-x01)*(x-x01) + (y-y01)*(y-y01)
#         if rSq<=r0*r0:
#             rho0[x,y,0] = 1.0+ sigma*3/r0
#             phiDiff0[x,y,0] = sigma*3/r0/rho0[x,y,0]#H/rho0[x,y,0]
#             phiInd[x,y,0] = 1-phiDiff0[x,y,0]
#         elif x<2*x0:     
#             rho0[x,y,0] = 1.0
#             phiDiff1[x,y,0] = H/rho0[x,y,0] #0
#             phiDiff2[x,y,0] = 1-phiDiff1[x,y,0]
#         elif r1Sq<=r01*r01:
#             rho0[x,y,0] = 1.0 + sigma*3/r01
#             phiDiff1[x,y,0] = H/rho0[x,y,0]
#             phiDiff2[x,y,0] = 1-phiDiff1[x,y,0]
#         else:
#             rho0[x,y,0] = 1.0
#             phiDiff0[x,y,0] = 0
#             phiInd[x,y,0] = 1-phiDiff0[x,y,0]


#rho0 = S
#rho0[10:40,10:40,10:40] = 0.5

#rho0[:, :50, :] = 1
#rho0[:, 1:-10, :] = 1
#rho0[:np.int(val.shape[0]*0.35), np.int(val.shape[1]*0.7):np.int(val.shape[1]*0.75), :np.int(val.shape[0]*0.35)] = 0
#vtk.append_data_set("rho0", rho0)

#rho1 = np.zeros(val.shape, dtype=float)
#rho1 = 1 - rho0
#rho1[:np.int(val.shape[0]*0.35), np.int(val.shape[1]*0.7):np.int(val.shape[1]*0.75), :np.int(val.shape[0]*0.35)] = 1

#vtk.append_data_set("phiInd", phiInd)
#vtk.append_data_set("phiDiff0", phiDiff0)
#vtk.append_data_set("phiDiff1", phiDiff1)
#vtk.append_data_set("phiDiff2", phiDiff2)

#phase 0 wetting
wet = np.ones(val.shape, dtype=float)
wet[:,:,:]=0.5*wet[:,:,:]
wet = wet*(geo == 1)
#wet[np.int(val.shape[0]*0.5):,:,:] = 0.;

vtk.append_data_set("wettability", wet)

#smooth_surf = 5*np.ones(val.shape, dtype=int)
#val[val > 0] = 1
#smooth_surf[:, 2:-2, :] = (val[:, 4:, :] + val[:, 3:-1, :] + val[:, 2:-2, :]  + val[:, 1:-3, :] + val[:, :-4, :] )
#vtk.append_data_set("smoothGeo", smooth_surf)

#plt.ion()

plt.figure(1)
plt.pcolormesh(wet[:, :, 10].T)
plt.gca().set_aspect('equal')
plt.title('wettability')
plt.colorbar()


#plt.figure(2)
#plt.pcolormesh(rho0[:, :, 0].T)
#plt.gca().set_aspect('equal')
#plt.title('rho0')
#plt.colorbar()


# plt.figure(3)
# plt.pcolormesh(phiDiff0[:, :, 0].T)
# plt.gca().set_aspect('equal')
# plt.title('phiDiff0')
# plt.colorbar()


# plt.figure(4)
# plt.pcolormesh(phiDiff1[:, :, 0].T)
# plt.gca().set_aspect('equal')
# plt.title('phiDiff1')
# plt.colorbar()


# plt.figure(5)
# plt.pcolormesh(phiDiff2[:, :, 0].T)
# plt.gca().set_aspect('equal')
# plt.title('phiDiff2')
# plt.colorbar()


# plt.figure(6)
# plt.pcolormesh(phiInd[:, :, 0].T)
# plt.gca().set_aspect('equal')
# plt.title('phiInd2')
# plt.colorbar()
plt.show()