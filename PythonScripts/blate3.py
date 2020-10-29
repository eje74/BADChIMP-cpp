from read_unstructured_grid import get_geo
import numpy as np
import matplotlib.pyplot
import vtklb


# Get the wing profile
geo = get_geo("/home/ejette/Programs/Python/blateTest3.vtu")


s = geo.shape
geo_test = np.zeros( (s[0], s[1]*3, s[2]), dtype=int)
for y in range(s[1]):
    geo_test[:,3*y:(3*y+3),:] = geo[:,y,:].reshape(s[0], 1, s[2])
geo = geo_test

geo_large = np.zeros((geo.shape[0]+40, geo.shape[1]+6, geo.shape[2]+20), dtype=int)
geo_large[20:(20+geo.shape[0]), :geo.shape[1], 10:(10 + geo.shape[2])] = geo
#geo_large[:,:5,:] = 0
geo_large[: ,:2, :] = 0
val = np.zeros((geo_large.shape[0]+100, geo_large.shape[1], geo_large.shape[2]), dtype=int)
val[:,:96 ,:] = 1
val[:,96:192,:] = 2
val[:, 192:, :] = 3
val[:geo_large.shape[0], :, :] = val[:geo_large.shape[0], :, :]*(geo_large == 0) 


vtk = vtklb.vtklb(val, "D3Q19", "x", "tmp", "/home/ejette/Programs/GitHub/BADChIMP-cpp/input/mpi/") 

# Setup boundary marker
bnd = np.zeros(val.shape, dtype=int)
#bnd[0, 1:-1, 1:-1] = 1 # Inlet x = 0
#bnd[-1, 1:-1, 1:-1] = 2 # Outlet x = -1
bnd[:, 0, :] = 3 # Right hand boundary y = 0
bnd[:, -1, :] = 4 # Left hand boundary y = -1
bnd[:, 1:-1, 0] = 5 # Bottom boundary z = 0
bnd[:, 1:-1, -1] = 6 # Top boundary z = 0
vtk.append_data_set("boundary", bnd)

smooth_surf = 5*np.ones(val.shape, dtype=int)
val[val > 0] = 1
smooth_surf[:, 2:-2, :] = (val[:, 4:, :] + val[:, 3:-1, :] + val[:, 2:-2, :]  + val[:, 1:-3, :] + val[:, :-4, :] )
vtk.append_data_set("smoothGeo", smooth_surf)

