#!/usr/bin/env python3
from vtklb import vtklb
import sys
from pathlib import Path
from Blood import Geo

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
geo.create(use_file=True)
geo.save('blood_geo.vtk')
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
vtk.append_data_set('pressure_boundary', geo.array('boundary'))
### Add the node to surface distance map 
vtk.append_data_set('distance', geo.array('distance'))
normal = geo.array('normal')
vtk.append_data_set('n_x', normal[...,0])
vtk.append_data_set('n_y', normal[...,1])
vtk.append_data_set('n_z', normal[...,2])
vtk.append_data_set('init_rho', geo.array(1.0))


