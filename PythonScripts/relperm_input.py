#!/usr/bin/env python3

from read_unstructured_grid import get_geo
from numpy import zeros, load, concatenate, sort, array, ones, prod
#import matplotlib.pyplot as plt
from scipy import ndimage
import vtklb
#import vtk_tools
import sys

size = [int(a) for a in sys.argv[1:4]]
nproc = [int(a) for a in sys.argv[4:7]]
Sw = float(sys.argv[7])
print('Size:',size,', Nproc:',nproc,', Sw:',Sw) 

def mpi_process_dist(size=None, nproc=None):
    pdist = zeros(size, dtype=int)
    for rank in range(prod(nproc)):
        ind = zeros(3,dtype=int)
        ind[0] = rank%nproc[0]
        ind[1] = ((rank - ind[0]) / nproc[0]) % nproc[1]
        ind[2] = int(rank / (nproc[1] * nproc[0]))
        lb = zeros(3,dtype=int)
        ub = zeros(3,dtype=int)
        for i in range(3):
            n_tmp = size[i] / nproc[i]
            n_rest = size[i] % nproc[i]
            lb[i] = ind[i] * n_tmp
            lb[i] += ind[i] if (ind[i] <= n_rest) else n_rest
            if ind[i] + 1 <= n_rest:
                n_tmp += 1
            ub[i] = lb[i] + n_tmp
        pdist[ lb[0]:ub[0], lb[1]:ub[1], lb[2]:ub[2] ] = rank+1
        #print(lb, ub)
    return pdist


directory="/cluster/home/janlv/BADChIMP-cpp/"
#file2 = "Porer-70kv2.modif.raw"  # void = 1, solid = 0
file2 = "Porer-70kv2.npy"  # void = 1, solid = 0
#geo1 = fromfile(file2, dtype=int8)
#geo2 = transpose(geo1.reshape((1422, 459, 455)))
geo2 = load(file2)
geo3 = geo2[:size[0],:size[1],:size[2]]

# Mirror geometry
geo4=zeros(geo3.shape, dtype=int)
geo4=geo3[::-1,:,:]
geo5=concatenate((geo3,geo4),axis=0)

geo3=geo5

#print('Distribute water evenly over solid surface')
# Distribute water evenly over solid surface
dist_transf = ndimage.distance_transform_edt(geo3)
dist_transf = array(dist_transf)
dist_transf[dist_transf==0]=1000
dist_ind = sort(dist_transf.flatten())
num_el = sum(dist_ind<1000)

# Water saturation
#Sw=0.9
print('Water saturation: '+str(Sw))
num_water = int(Sw*num_el)

#number of nodes less than cut
num_less_cut = sum(dist_ind<dist_ind[num_water])
S = zeros(geo3.shape, dtype=float)
S[dist_transf < dist_ind[num_water]] = 1.0

num_on_cut = sum(dist_ind == dist_ind[num_water])
# Local water concentration field 
S[dist_transf == dist_ind[num_water]] = (1.0*(num_water - num_less_cut))/num_on_cut

# BADChIMP require solid = 0, void >= 1
geo3=1-geo3  # solid = 1, void = 0
geo=geo3

# MPI load partitioning
rank = mpi_process_dist(size=geo.shape, nproc=nproc) 
val=rank*(geo == 0)  # val = 0 for solid nodes
#for n in range(1,prod(nproc)+1):
#    print(str(n),': ',sum(val==n))

# Periodic in x,y,z
vtk = vtklb.vtklb(val, "D3Q19", "xyz", "tmp", directory+"input/mpi/") 

# Setup boundary marker
# Do we need this?
bnd = zeros(val.shape, dtype=int)
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

qSrc = zeros(val.shape, dtype=int)
vtk.append_data_set("source", qSrc)

rho0 = S
vtk.append_data_set("rho0", rho0)

rho1 = zeros(val.shape, dtype=float)
rho1 = 1 - rho0
vtk.append_data_set("rho1", rho1)

# phase 0 wetting
wet = ones(val.shape, dtype=float)
wet[:,:,:]=0.5*wet[:,:,:]
wet = wet*(geo == 1)
#wet[int(val.shape[0]*0.5):,:,:] = 0.;
vtk.append_data_set("wettability", wet)

print()

#plt.ion()

#plt.figure(1)
#plt.pcolormesh(wet[:, :, 10].T)
#plt.gca().set_aspect('equal')
#plt.title('wettability')
#plt.colorbar()


#plt.figure(2)
#plt.pcolormesh(rho0[:, :, 0].T)
#plt.gca().set_aspect('equal')
#plt.title('rho0')
#plt.colorbar()

#plt.show()
