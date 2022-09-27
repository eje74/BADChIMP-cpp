#!/usr/bin/env python3
from vtklb import vtklb
import sys
from pathlib import Path
from Blood import Geo

import matplotlib.pyplot as plt
import matplotlib

tube_list = [2, 3, 5, 1] # Geometry specific

narg = len(sys.argv) - 1
if narg < 2:
    raise SystemExit('\n  Usage:\n    blood_geo.py <grid resolution> <nproc> [connect=0] [workers=2]\n')
dx       = float(sys.argv[1])
nproc    =   int(sys.argv[2])
connect  = narg>2 and int(sys.argv[3]) or 0
workers  = narg>3 and int(sys.argv[4]) or 2

scriptdir = Path(sys.argv[0]).resolve().parent # Abolute path of script-folder
chimpdir = scriptdir.parent                    # Assume script-folder in BADChIMP-folder

bfiles = scriptdir/'bloodfiles'
mesh = bfiles/'fontan_stor_v5_L1_bl050_sl3.vtu'
centerline = bfiles/'centerlines.vtp'
geo = Geo(dx=dx, surface=mesh, centerline=centerline, tube_list=tube_list, workers=workers, nproc=nproc, connected=connect, echo=True)
#geo.create(use_file=True)
geo.create(use_file=False)
geo.save('blood_geo_test.vtk')
fluid_nodes = geo.copy().threshold(value=0.5, scalars='wall', invert=True)
fluid_nodes.save('fluid_nodes.vtk')
print()

### Geometry:
###     process number at fluid nodes, 0 at solid nodes
vtk = vtklb(geo.array('proc'), 'D3Q19', 'xyz', 'tmp', f"{chimpdir/'input'/'mpi'}/") 
#vtk.append_data_set('boundary', geo.array('boundary'))
###  Pressure boundary:
###     0 everywhere execpt at the pressure boundaries and behind them 
###     The boundaries are numbered from 1 and up. The non-fluid side
###     of the boundary is numbered -1. The boundary nodes must give the
###     distance to the boundary plane, and the normal vector pointing inwards

p_boundaries = (fluid_nodes.array('boundary')).astype(int)

vtk.append_data_set('pressure_boundary', p_boundaries)

vtk.append_data_set("signed_distance", fluid_nodes.array('distance'))

# -------------------------------------------------------------------
#      S E T   T H E   I N L E T / O U T L E T   G E O M E T R Y
# -------------------------------------------------------------------
vtk.begin_point_data_subset_section()
# ------------------------------------------------------------------- Set mask
mask = fluid_nodes.array('boundary') > 0

# ------------------------------------------------------------------- Set normal
normal = -fluid_nodes.array('normal')
vtk.append_data_subset('boundary_normal_x', normal[..., 0], mask)
vtk.append_data_subset('boundary_normal_y', normal[..., 1], mask)
vtk.append_data_subset('boundary_normal_z', normal[..., 2], mask)


# ------------------------------------------------------------------- Set distance from pressure boundary
vtk.append_data_subset("boundary_distance", fluid_nodes.array('distance'), mask)


"""
vtk.append_data_set('pressure_boundary', geo.array('boundary'))
### Add the node to surface distance map 
vtk.append_data_set('distance', geo.array('distance'))
normal = geo.array('normal')
vtk.append_data_set('boundary_normal_x', normal[...,0])
vtk.append_data_set('boundary_normal_y', normal[...,1])
vtk.append_data_set('boundary_normal_z', normal[...,2])
vtk.append_data_set('init_rho', geo.array(1.0))
"""

cm= matplotlib.cm.viridis #'viridis'#matplotlib.cm.YlOrBr_r

dist = fluid_nodes.array('distance')
boundary = fluid_nodes.array('boundary')

plt.figure(1, figsize=(8, 8))
#ax = plt.pcolormesh(phi[:,:,0].transpose(), cmap=grayify_cmap(cm))
ax = plt.pcolormesh(boundary[:,:,5], cmap=cm, rasterized=True)
cbar=plt.colorbar(shrink=0.8, pad= 0.06)


plt.show()