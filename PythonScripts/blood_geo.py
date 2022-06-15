#!/usr/bin/env python3
from vtklb import vtklb
import sys
from pathlib import Path
from Blood import Geo

tube_list = [2, 3, 5, 1] # Geometry specific


narg = len(sys.argv) - 1
if narg < 4:
    raise SystemExit('\n  Usage:\n    blood_geo.py <path to surface mesh file> <path to centerline file> <grid resolution> <nproc> [workers=2]\n')
meshpath =  Path(sys.argv[1])
clpath   =  Path(sys.argv[2])
dx       = float(sys.argv[3])
nproc    =   int(sys.argv[4])
workers  = narg>4 and int(sys.argv[5]) or 2

geo = Geo(dx=dx, surface=meshpath, centerline=clpath, tube_list=tube_list, workers=workers, nproc=nproc, connected=False, echo=True)
geo.create()

# Assume script in a sub-folder of BADChIMP-folder
chimpdir = list(Path(sys.argv[0]).resolve().parents)[1]

vtk = vtklb(geo.array('proc'), 'D3Q19', 'xyz', 'tmp', f"{chimpdir/'input'/'mpi'}/") 
vtk.append_data_set('boundary', geo.array('boundary'))
vtk.append_data_set('init_rho', geo.array(1.0))


