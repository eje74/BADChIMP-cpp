#!/usr/bin/env python3

from pyvista import read as pvread, UniformGrid, PolyData
from numpy import array, ceil, sum as npsum, ones, double as npdouble, mean
from scipy.spatial import KDTree
from pathlib import Path
import matplotlib.pyplot as pl
import sys
#import centerline

narg = len(sys.argv) - 1
if narg < 3:
    raise SystemExit('\n  Usage:\n    distance_map <path to surface mesh file> <path to centerline file> <grid resolution> [number of processors]\n')
meshpath = Path(sys.argv[1])
clpath = Path(sys.argv[2])
dx = float(sys.argv[3])
workers = 1
if narg > 3:
    workers = int(sys.argv[4])

wall = 2  # Number of wall nodes

# Read surface mesh
mesh = pvread(meshpath)
mesh = mesh.extract_surface()
# Only use boundary cells, i.e. cells with 'CellEntityIds' > 0 (created by VMTK)
cell_id = 'CellEntityIds' 
id_scalar = mesh.cell_data.get(cell_id)
if id_scalar is None:
    raise SystemExit(f'The mesh-file is missing the {cell_id} scalar that identify boundary cells, aborting...\n')
bndry_cells = id_scalar > 0
mesh_cells   = mesh.cell_centers().points.astype(npdouble)[bndry_cells]
mesh_normals = mesh.cell_normals[bndry_cells]
mesh_cell_id = mesh[cell_id][bndry_cells]

# Create uniform grid
xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
size = array((xmax-xmin, ymax-ymin, zmax-zmin))
size += 2*wall*dx
dim = ceil(size/dx).astype(int)
grid = UniformGrid()
grid.dimensions = dim + 1 # Because we want to add cell data (not point data)
grid.origin = array((xmin, ymin, zmin))-wall*dx
grid.spacing = (dx, dx, dx) # Evenly spaced grids
grid_cells = grid.cell_centers().points.astype(npdouble)
print(dim)

# Return shortest distance from grid-cell to mesh surface cell together with mesh index
dist_idx = array(KDTree(mesh_cells).query(grid_cells, workers=workers))

# Inside or outside surface? 
idx = dist_idx[1].astype(int) # Mesh-surface index
norm_dot_vec = mesh_normals[idx] * (mesh_cells[idx]-grid_cells)
inside = npsum(norm_dot_vec, axis=1) > 0
sign = -ones(inside.shape)
sign[inside] = 1  # Positive values inside surface 

# Add distance from grid node to surface as cell_data
distance = 'distance'
grid[distance] = sign * dist_idx[0]

# Add boundary markers
marker = 'boundary'
grid[marker] = mesh_cell_id[idx]
grid[marker][grid[distance]>0] = 0  # Mark interior (fluid) nodes as boundary 0

# The fluid nodes is confined by a wall of boundary nodes. 
# Threshold is used to remove obsolete nodes
# Threshold returns an unstructured grid
ugrid = grid.threshold(value=-wall*dx, scalars=distance)

# Return shortest distance from grid-cells to centerline together with centerline index
#cl = centerline.read(clpath)
cl = pvread(clpath)
#lines = [PolyData(cl.cell_points(i)) for i in range(cl.n_cells)]
#npoints = [cl.cell_n_points(i) for i in range(cl.n_cells)]
#start = cumsum([0]+npoints)
#lines = [arange(start[i],start[i+1]) for i in range(len(start)-1)]
#cl.compute_arc_length()
#length = []
#cl_points = cl.points.astype(npdouble)
ugrid_cells = ugrid.cell_centers().points.astype(npdouble)
for n in range(cl.n_cells):
    cl_points = cl.cell_points(n).astype(npdouble)
    dist_idx = array(KDTree(cl_points).query(ugrid_cells, workers=workers))
    idx = dist_idx[1].astype(int) # centerline point index
    ugrid[f'line_{n}_distance'] = dist_idx[0]
    pl.figure(n)
    m = mean(dist_idx[0])
    pl.hist(dist_idx[0], 100)
    print(f'line {n}: mean = {m}')
    pl.draw()
pl.show()

#ugrid['slice'] = idx

# Save final grid
map = meshpath.parent/'distance_map.vtk'
ugrid.save(map)
print(f'  Saved {map}')


