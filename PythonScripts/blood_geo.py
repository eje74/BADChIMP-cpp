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
geo.create()
print()

vtk = vtklb(geo.array('proc'), 'D3Q19', 'xyz', 'tmp', f"{chimpdir/'input'/'mpi'}/") 
vtk.append_data_set('boundary', geo.array('boundary'))
vtk.append_data_set('init_rho', geo.array(1.0))


