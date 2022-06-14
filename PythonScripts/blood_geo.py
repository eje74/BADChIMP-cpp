#!/usr/bin/env python3
from vtklb import vtklb
import numpy as np
#import matplotlib.pyplot as plt
import sys
from pathlib import Path
import Blood
#from blood_grid import Grid, Process_map, ijk_array
#import pickle

loadfile = True
dim = 3
connected = False #True
echo = True
tube_list = [2, 3, 5, 1]

chimpdir = Path.home()/"github/BADChIMP-cpp"

narg = len(sys.argv) - 1
if narg < 4:
    raise SystemExit('\n  Usage:\n    blood_geo.py <path to surface mesh file> <path to centerline file> <grid resolution> <nproc> [workers=2] [filename=blood_geo.vtk]\n')
meshpath =  Path(sys.argv[1])
clpath   =  Path(sys.argv[2])
dx       = float(sys.argv[3])
nproc    =   int(sys.argv[4])
workers  = narg>4 and int(sys.argv[5]) or 2

grid = Blood.Grid(dx=dx, surface=meshpath, centerline=clpath, tube_list=tube_list, workers=workers, echo=echo)
pmap = Blood.Process_map(nproc=nproc, connected=connected, echo=echo)

dx_str = f'{int(dx)}-{round((dx%1)*100):2d}'
dim_str = f'{"x".join([str(d) for d in grid.dim])}'
gridfile = Path(f'{meshpath.parent/("LB__"+meshpath.stem)}__dx_{dx_str}__dim_{dim_str}.vtk')

if loadfile and gridfile.is_file():
    pmap.read_grid(gridfile)
else:
    grid.create()
    grid.save(gridfile)
    pmap.grid = grid.grid

grid.grid = pmap.distribute_load()

vtk = vtklb(grid.array('proc'), 'D3Q19', 'xyz', 'tmp', f"{chimpdir/'input'/'mpi'}/") 
vtk.append_data_set('boundary', grid.array('boundary'))
vtk.append_data_set('init_rho', grid.array(1.0))


