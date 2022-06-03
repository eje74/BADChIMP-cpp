#!/usr/bin/env python3

from pyvista import read as pvread, UniformGrid, UnstructuredGrid
from numpy import argsort, array, ceil, min as npmin, max as npmax, nonzero, sort, sum as npsum, ones, double as npdouble, mean, zeros, sum as npsum
from scipy.spatial import KDTree
from pathlib import Path
import matplotlib.pyplot as pl
import sys

group_processes = True
nproc = 10
tube_list = [2, 3, 5, 1]
echo = True

class Cluster:
    def __init__(self, grid, n=-1) -> None:
        self.grid = grid
        self.n = n
        self.body = self.grid.connectivity().split_bodies()

    def __str__(self) -> str:
        return f'<grid:{self.grid.n_cells}, body:{self.num_bodies()}, n:{self.n}, is_solid:{self.is_solid()}>'

    def __repr__(self) -> str:
        return self.__str__()

    def is_solid(self):
        if len(self.body) == 1:
            return True
        return False
            
    def num_bodies(self):
        return len(self.body)

    def bodies(self):
        return [self.body[k] for k in self.body.keys()]

    def volume(self):
        return [body.volume for body in self.bodies()]
        #return [self.body[k].volume for k in self.body.keys()]

    def bodies_by_volume(self):
        vol = self.volume()
        keys = self.body.keys()
        return [self.body[keys[i]] for i in argsort(vol)]

    def merge(self, bodies):
        self.grid = self.grid.merge(bodies, merge_points=True)
        self.body = self.grid.connectivity().split_bodies()

    def distance_to_centerline(self, centerline): 
        cl_points = centerline.cell_points(self.n).astype(npdouble)
        dist_idx = array(KDTree(cl_points).query(self.grid.cell_centers().points, workers=workers))
        return dist_idx[0], dist_idx[1].astype(int)

#    def cell_centers(self):
#        return self.grid.cell_centers().points.astype(npdouble)

def reconnect_clusters(grid, values, scalar='tube', echo=True):
    # Check if tubes are unconnected
    #echo and print(f'  Check for unconnected {scalar} clusters ...')
    split_clusters = []
    solid_clusters = []
    for n in values:
        cluster = Cluster(grid.threshold(value=(n,n), scalars=scalar), n=n)
        if cluster.is_solid():
            solid_clusters.append(cluster)
        else:
            split_clusters.append(cluster)

    # Join unconnected clusters with solid clusters
    for split in split_clusters:
        joined = [(Cluster(split.grid + solid.grid), solid.n) for solid in solid_clusters]
        joined = [j for j in joined if j[0].num_bodies()==2]
        for j,n in joined:
            echo and print(f'  Unconnected {scalar} cluster {split.n} is merged into solid {scalar} cluster {n} ...')
            bodies = j.bodies()
            mrg_ind = [i for i,b in enumerate(bodies) if b[scalar].max()-b[scalar].min()!=0][0]
            merged = bodies[mrg_ind]
            solid_rest = bodies[int(not mrg_ind)]
            merged[scalar] = n * ones(merged.n_cells, dtype=int)            
            ind = [i for i,s in enumerate(solid_clusters) if s.n==n][0]
            solid_clusters[ind] = Cluster(merged, n=n)
            solid_clusters.append(Cluster(solid_rest, n=split.n))
            break
    return solid_clusters


# Merge slices into processes
def group_slices(nproc, slices):
    size = len(slices)
    proc_size = ceil(size/nproc)
    proc = -1*ones(size, dtype=int)
    stop = npmax(slices)
    a = b = 0
    p = 0
    while p < nproc:
        mask = (slices >= a) * (slices <= b)
        maskb = (slices >= a) * (slices <= b+1)
        if abs(proc_size-npsum(mask)) < abs(proc_size-npsum(maskb)) or b == stop:
            #print(f'p: {p}, a,b: {a},{b}  {npsum(mask)}, {proc_size}, {npsum(mask)-proc_size}')
            a = b + 1 
            proc[mask] = p
            p += 1
        b += 1
    if npsum(proc==-1) > 0:
        print(f'  WARNING! {npsum(proc==-1)} nodes not assigned to a process!')
    return proc



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
#wall = 2  # Number of wall nodes

def grid_from_surface(surface, wall=2, distance='distance', marker='boundary', echo=True):
    # Read surface mesh
    mesh = pvread(surface)
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
    echo and print('  Creating uniform grid of dimension ', end='')
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
    size = array((xmax-xmin, ymax-ymin, zmax-zmin))
    size += 2*wall*dx
    dim = ceil(size/dx).astype(int)
    grid = UniformGrid()
    grid.dimensions = dim + 1 # Because we want to add cell data (not point data)
    grid.origin = array((xmin, ymin, zmin))-wall*dx
    grid.spacing = (dx, dx, dx) # Evenly spaced grids
    grid_cells = grid.cell_centers().points.astype(npdouble)
    echo and print(dim)

    # Return shortest distance from grid-cell to mesh surface cell together with mesh index
    echo and print('  Calculating distance from grid nodes to the surface ...')
    dist_idx = array(KDTree(mesh_cells).query(grid_cells, workers=workers))

    # Inside or outside surface? 
    echo and print(f'  Creating unstructured grid with outer wall of {wall} voxels ...')
    idx = dist_idx[1].astype(int) # Mesh-surface index
    norm_dot_vec = mesh_normals[idx] * (mesh_cells[idx]-grid_cells)
    inside = npsum(norm_dot_vec, axis=1) > 0
    sign = -ones(inside.shape)
    sign[inside] = 1  # Positive values inside surface 

    # Add distance from grid node to surface as cell_data
    grid[distance] = sign * dist_idx[0]

    # Add boundary markers
    grid[marker] = mesh_cell_id[idx]
    grid[marker][grid[distance]>0] = 0  # Mark interior (fluid) nodes as boundary 0

    # The fluid nodes is confined by a wall of boundary nodes. 
    # Threshold is used to remove obsolete nodes
    # Threshold returns an unstructured grid
    ugrid = grid.threshold(value=-wall*dx, scalars=distance)

    return ugrid

ugrid = grid_from_surface(meshpath, echo=echo)

def create_tube(ugrid, centerline, selection=[], echo=True):
    echo and print('  Creating tubes ...')
    cl = pvread(centerline)
    ugrid_cells = ugrid.cell_centers().points.astype(npdouble)
    tubes = ones((cl.n_cells, ugrid.n_cells), dtype=int)
    for n in range(cl.n_cells):
        cl_points = cl.cell_points(n).astype(npdouble)
        dist = array(KDTree(cl_points).query(ugrid_cells, workers=workers))[0]
        m = mean(dist)
        remove = dist>m
        tubes[n,:] = n
        tubes[n,remove] = -1

    # Join tubes and select slices 
    tube = -1*ones(ugrid.n_cells, dtype=int)
    c = 0
    for tb in tubes[selection]:
        mask = (tb>=0)*(tube<0)
        tube[mask] = tb[mask]
    #ugrid['tube'] = tube

    return tube

ugrid['tube'] = create_tube(ugrid, clpath, selection=tube_list)

solid_tubes = reconnect_clusters(ugrid, tube_list, scalar='tube')

def create_slices(centerline, tubes, echo=True):
    # Make slices
    echo and print('  Calculating slice indices ...')
    cl = pvread(centerline)
    c = 0
    for tube in tubes:
        #tube.grid['tube'] = tube.n * ones(tube.grid.n_cells, dtype=int)
        dist, idx = tube.distance_to_centerline(cl)
        tube.grid['slice'] = 1 + c + idx - npmin(idx)
        c = npmax(tube.grid['slice'])
        #print(f'min, max, c : {npmin(idx)}, {npmax(idx)}, {c}')
    ugrid = UnstructuredGrid().merge([tb.grid for tb in tubes])
    return ugrid

ugrid = create_slices(clpath, solid_tubes)

echo and print('  Make process groups ...')
ugrid['proc'] = group_slices(nproc, ugrid['slice'])
ugrid['proc2'] = ugrid['proc']

if group_processes: 
    echo and print('  Group split process groups ...')
    proc_solid = reconnect_clusters(ugrid, list(range(nproc)), scalar='proc')
    ugrid = UnstructuredGrid().merge([solid.grid for solid in proc_solid], merge_points=True)

# Save final grid
ugrid.save(savename)
print(f'  Saved {savename}')

#print(f'ugrid.n_cells: {ugrid.n_cells}, ugrid2.n_cells: {ugrid2.n_cells}, ugrid3.n_cells: {ugrid3.n_cells}')



