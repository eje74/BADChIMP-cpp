#!/usr/bin/env python3
from vtklb import vtklb
import numpy as np
#import matplotlib.pyplot as plt
import sys
from pathlib import Path
from blood_grid import Grid

dim = 3
connected = True
echo = True
tube_list = [2, 3, 5, 1]

chimpdir = Path.home()/"github/BADChIMP-cpp"

narg = len(sys.argv) - 1
if narg < 4:
    raise SystemExit('\n  Usage:\n    blood_geo.py <path to surface mesh file> <path to centerline file> <grid resolution> <nproc> [workers=1] [filename=blood_geo.vtk]\n')
meshpath =  Path(sys.argv[1])
clpath   =  Path(sys.argv[2])
dx       = float(sys.argv[3])
nproc    =   int(sys.argv[4])
workers  = narg>4 and int(sys.argv[5]) or 1
savename = meshpath.parent/(narg>5 and sys.argv[6] or 'blood_geo.vtk')

#print(f'surface:{meshpath}, centerline:{clpath}, tube_list:{tube_list}, nproc:{nproc}, dx:{dx}, connected:{connected}, echo:{echo}')

grid = Grid(dx=dx, surface=meshpath, centerline=clpath, tube_list=tube_list, nproc=nproc, workers=workers, echo=echo, connected=connected)
grid.create()
grid.save(f'blood_geo__np-{nproc}.vtk')
vtk = vtklb(grid.array('proc'), 'D3Q19', 'xyz', 'tmp', f"{chimpdir/'input'/'mpi'}/") 
vtk.append_data_set('boundary', grid.array('boundary'))
vtk.append_data_set('init_rho', grid.array(1.0))


