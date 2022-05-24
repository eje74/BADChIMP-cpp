#!/usr/bin/env python3

from re import U
from pyvista import read as pvread, UniformGrid, UnstructuredGrid
from numpy import argsort, array, ceil, min as npmin, max as npmax, nonzero, sum as npsum, ones, double as npdouble, mean, zeros, sum as npsum
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

savename = meshpath.parent/'distance_map.vtk'
wall = 2  # Number of wall nodes

# Read surface mesh
mesh = pvread(meshpath)
mesh = mesh.extract_surface()
# Only use boundary cells, i.e. cells with 'CellEntityIds' > 0 (created by VMTK)
cell_id = 'CellEntityIds' 
id_scalar = mesh.cell_data.get(cell_id)
if id_scalar is None:
    raise SystemExit(f'  ERROR! The mesh-file is missing the {cell_id} scalar that identify boundary cells, aborting...\n')
bndry_cells = id_scalar > 0
mesh_cells   = mesh.cell_centers().points.astype(npdouble)[bndry_cells]
mesh_normals = mesh.cell_normals[bndry_cells]
mesh_cell_id = mesh[cell_id][bndry_cells]

# Create uniform grid
print('  Creating uniform grid of dimension ', end='')
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
print('  Calculating distance from grid nodes to the surface ...')
dist_idx = array(KDTree(mesh_cells).query(grid_cells, workers=workers))

# Inside or outside surface? 
print(f'  Creating unstructured grid with outer wall of {wall} voxels ...')
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

# # Join centerline points
# cl = pvread(clpath)
# pts = []
# [pts.extend(cl.cell_points(i)) for i in range(cl.n_cells)]
# pts = array(pts)

# # Remove duplicate points from the centerlines
# pairs = KDTree(pts).query_pairs(r=0.98)
# remove = list(set([b for a,b in pairs]))
# keep = ones(len(pts))
# keep[remove] = 0
# cl_points = pts[keep.astype(bool)]

# # Return shortest distance from grid-cells to centerline together with centerline index
# ugrid_cells = ugrid.cell_centers().points.astype(npdouble)
# dist_idx = array(KDTree(cl_points).query(ugrid_cells, workers=workers))
# ugrid['slice'] = dist_idx[1]

# ugrid.save(savename)
# print(f'  Saved {savename}')

print('  Creating tubes ...')
cl = pvread(clpath)
ugrid_cells = ugrid.cell_centers().points.astype(npdouble)
tubes = ones((cl.n_cells, ugrid.n_cells), dtype=int)
# slices = ones((cl.n_cells, ugrid.n_cells), dtype=int)
#length2 = zeros(cl.n_cells)
for n in range(cl.n_cells):
    cl_points = cl.cell_points(n).astype(npdouble)
    # length2[n] = npsum( npsum((cl_points[1:,:]-cl_points[:-1,:])**2, axis=1) )
    dist_idx = array(KDTree(cl_points).query(ugrid_cells, workers=workers))
    # idx = dist_idx[1].astype(int) # centerline point index
    dist = dist_idx[0]
    # ugrid[f'line_{n}_distance'] = dist
    m = mean(dist)
    remove = dist>m
    tubes[n,:] = n
    tubes[n,remove] = -1
    # slices[n,:] = idx 
    # slices[n,:][remove] = -1

# # Join tubes and select slices 
# is_tube = tubes > 0
# sizes = npsum(is_tube, axis=1)
tube_list = [2, 3, 5, 1]
tube = -1*ones(ugrid.n_cells, dtype=int)
# slice = -1*ones(ugrid.n_cells, dtype=int)
c = 0
for n in tube_list:
    mask = (tubes[n,:]>=0)*(tube<0)
    tube[mask] = tubes[n,mask]
    # new_slice[mask] -= npmin(new_slice[mask])
    # slice[mask] = new_slice[mask] + c
    # c = npmax(slice)+1
ugrid['tube'] = tube
# ugrid['slice'] = slice

# Check if tubes are unconnected
print('  Check for unconnected tubes ...')
split_tubes = []
solid_tubes = []
GRID, BODY, N = 0, 1, 2
for i,n in enumerate(tube_list):
    th_grid = ugrid.threshold(value=(n,n), scalars='tube')
    bodies = th_grid.connectivity().split_bodies()
    gbn_list = [th_grid, bodies, n]
    if len(bodies) > 1:
        split_tubes.append(gbn_list)
        #print(f'Tube {n}(i:{i}) is split in {len(bodies)} bodies')
    else:
        solid_tubes.append(gbn_list)

# Join unconnected tubes with solid tubes
new_solid = []
for split in split_tubes:
    for solid in solid_tubes:
        joined = split[GRID] + solid[GRID] 
        bodies = joined.connectivity().split_bodies()
        if len(bodies) == 2:
            print(f'  Unconnected tube {split[N]} is merged into solid tube {solid[N]} ...')
            # List of split body volumes. Merge all but the largest volume  
            vol = [split[BODY][k].volume for k in split[BODY].keys()]
            # List of split bodies sorted by increasing volume
            body_list = [split[BODY][split[BODY].keys()[i]] for i in argsort(vol)]
            # Merge the split bodies with the solid grid (except for the largest one)
            to_merge = UnstructuredGrid()
            solid[GRID].merge(body_list[:-1], merge_points=True, inplace=True)
            solid[GRID]['tube'] = solid[N]*ones(solid[GRID].n_cells, dtype=int)
            cl_points = cl.cell_points(solid[N]).astype(npdouble)
            dist_idx = array(KDTree(cl_points).query(solid[GRID].cell_centers().points, workers=workers))
            solid[GRID]['slice'] = dist_idx[1]
            # Create new solid grid from the largest split grid
            new_grid = body_list[-1]
            new_bodies = new_grid.connectivity().split_bodies()
            new_solid = [new_grid, new_bodies, split[N]]
            if len(new_bodies) > 1:
                msg =   '  ERROR!\n'
                msg += f'  Joining {len(body_list)-1} pieces of split tube {split[N]} into solid tube {solid[N]} did not succeed!\n'
                msg += f'  Tube {split[N]} has {len(new_bodies)} bodies, and {solid[N]} has {len(solid[GRID].connectivity().split_bodies())} bodies\n'
                raise SystemError(msg)
solid_tubes.append(new_solid)

# Make slices
print('  Calculate slice indices ...')
c = 0
for tube in solid_tubes:
    cl_points = cl.cell_points(tube[N]).astype(npdouble)
    tube_cells = tube[GRID].cell_centers().points.astype(npdouble)
    dist_idx = array(KDTree(cl_points).query(tube_cells, workers=workers))
    idx = dist_idx[1].astype(int) # centerline point index
    tube[GRID]['slice'] = 1 + c + idx - npmin(idx)
    c = npmax(tube[GRID]['slice'])
    print(f'min, max, c : {npmin(idx)}, {npmax(idx)}, {c}')

ugrid2 = UnstructuredGrid().merge([solid[GRID] for solid in solid_tubes])
#ugrid2.save(savename.with_name('ugrid2.vtk'))
print(f'ugrid.n_cells: {ugrid.n_cells}, ugrid2.n_cells: {ugrid2.n_cells}')
# Join slices into process-blocks
#nproc = 10
#proc_size = ceil(tube.size/nproc)


# Save final grid
ugrid2.save(savename)
print(f'  Saved {savename}')


