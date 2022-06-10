#!/usr/bin/env python3

from pyvista import read as pvread, UniformGrid, UnstructuredGrid
from numpy import argsort, array, ceil, min as npmin, max as npmax, nonzero, sort, sum as npsum, ones, double as npdouble, mean, zeros, sum as npsum
from scipy.spatial import KDTree
from pathlib import Path
import matplotlib.pyplot as pl

#====================================================================================
class Cluster:
#====================================================================================
    #--------------------------------------------------------------------------------
    def __init__(self, grid, n=-1) -> None:
    #--------------------------------------------------------------------------------
        self.grid = grid
        self.n = n
        self.body = self.grid.connectivity().split_bodies()

    #--------------------------------------------------------------------------------
    def __str__(self) -> str:
    #--------------------------------------------------------------------------------
        return f'<grid:{self.grid.n_cells}, body:{self.num_bodies()}, n:{self.n}, is_solid:{self.is_solid()}>'

    #--------------------------------------------------------------------------------
    def __repr__(self) -> str:
    #--------------------------------------------------------------------------------
        return self.__str__()

    #--------------------------------------------------------------------------------
    def is_solid(self):
    #--------------------------------------------------------------------------------
        if len(self.body) == 1:
            return True
        return False
            
    #--------------------------------------------------------------------------------
    def num_bodies(self):
    #--------------------------------------------------------------------------------
        return len(self.body)

    #--------------------------------------------------------------------------------
    def bodies(self):
    #--------------------------------------------------------------------------------
        return [self.body[k] for k in self.body.keys()]

    #--------------------------------------------------------------------------------
    def volume(self):
    #--------------------------------------------------------------------------------
        return [body.volume for body in self.bodies()]

    #--------------------------------------------------------------------------------
    def bodies_by_volume(self):
    #--------------------------------------------------------------------------------
        vol = self.volume()
        keys = self.body.keys()
        return [self.body[keys[i]] for i in argsort(vol)]

    #--------------------------------------------------------------------------------
    def merge(self, bodies):
    #--------------------------------------------------------------------------------
        self.grid = self.grid.merge(bodies, merge_points=True)
        self.body = self.grid.connectivity().split_bodies()

    #--------------------------------------------------------------------------------
    def distance_to_centerline(self, centerline, workers=1): 
    #--------------------------------------------------------------------------------
        cl_points = centerline.cell_points(self.n).astype(npdouble)
        dist_idx = array(KDTree(cl_points).query(self.grid.cell_centers().points, workers=workers))
        return dist_idx[0], dist_idx[1].astype(int)

#    def cell_centers(self):
#        return self.grid.cell_centers().points.astype(npdouble)


#====================================================================================
class Grid:
#====================================================================================
    #--------------------------------------------------------------------------------
    def __init__(self, dx=None, surface=None, centerline=None, wall=2, echo=True, 
                 tube_list=None, connected=True, nproc=1, workers=1) -> None:
    #--------------------------------------------------------------------------------
        self.echo = echo
        self.surface = surface
        self.dx = dx
        self.wall = wall
        self.tube_list = tube_list
        self.connected = connected
        self.nproc = nproc
        self.centerline = centerline
        self.workers = workers
        self.grid = None
        self.ijk = None
        self.solid = {}

    #--------------------------------------------------------------------------------
    def create(self): 
    #--------------------------------------------------------------------------------
        tube = 'tube'
        slice = 'slice'
        proc = 'proc'
        self.centerline = pvread(self.centerline)
        self.add_distance_map('distance', marker='boundary', cell_id='CellEntityIds')
        self.add_tubes(tube)
        self.reconnect_clusters(values=self.tube_list, scalar=tube)
        self.add_slices(slice, tube=tube)
        self.add_processes(proc, start=1, slice=slice)
        if self.connected: 
            self.reconnect_clusters(values=list(range(1,self.nproc+1)), scalar=proc)
            self.grid = UnstructuredGrid().merge([solid.grid for solid in self.solid[proc]], merge_points=True)
        return self

    #--------------------------------------------------------------------------------
    def save(self, *args, **kwargs):
    #--------------------------------------------------------------------------------
        self.grid.save(*args, **kwargs)

    #--------------------------------------------------------------------------------
    def array(self, value):
    #--------------------------------------------------------------------------------
        arr = zeros(self.dim)
        if isinstance(value, str):
            data = self.grid[value]
        else:
            # Assume number...
            data = value*ones(self.grid.n_cells)
        ijk = self.grid['ijk']
        arr[ijk[:,0], ijk[:,1], ijk[:,2]] = data
        return arr

    # #--------------------------------------------------------------------------------
    # def ones(self):
    # #--------------------------------------------------------------------------------
    #     arr = zeros(self.dim)
    #     arr[self.ijk[:,0], self.ijk[:,1], self.ijk[:,2]] = ones(self.grid.n_cells)
    #     return arr

    #--------------------------------------------------------------------------------
    def add_distance_map(self, name, marker=None, cell_id=None):
    #--------------------------------------------------------------------------------
        # Read surface mesh
        mesh = pvread(self.surface)
        mesh = mesh.extract_surface()
        # Only use boundary cells, i.e. cells with 'CellEntityIds' > 0 (created by VMTK)
        # cell_id = 'CellEntityIds' 
        id_scalar = mesh.cell_data.get(cell_id)
        if id_scalar is None:
            raise SystemExit(f'  ERROR! The mesh-file is missing the {cell_id} scalar that identify boundary cells, aborting...\n')
        bndry_cells = id_scalar > 0
        mesh_cells   = mesh.cell_centers().points.astype(npdouble)[bndry_cells]
        mesh_normals = mesh.cell_normals[bndry_cells]
        mesh_cell_id = mesh[cell_id][bndry_cells]

        # Create uniform grid
        self.echo and print('  Creating uniform grid of dimension ', end='')
        xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
        size = array((xmax-xmin, ymax-ymin, zmax-zmin))
        dx = self.dx
        size += 2*self.wall*dx
        self.dim = ceil(size/dx).astype(int)
        grid = UniformGrid()
        grid.dimensions = self.dim + 1 # Because we want to add cell data (not point data)
        grid.origin = array((xmin, ymin, zmin))-self.wall*dx
        grid.spacing = (dx, dx, dx) # Evenly spaced grids
        grid_cells = grid.cell_centers().points.astype(npdouble)
        ijk = ( (grid_cells-grid.origin)/dx ).astype(int)
        grid['ijk'] = ijk
        #grid['BLOCK_I'] = ijk[:,0]
        #grid['BLOCK_J'] = ijk[:,1]
        #grid['BLOCK_K'] = ijk[:,2]
        self.echo and print(self.dim)

        # Return shortest distance from grid-cell to mesh surface cell together with mesh index
        self.echo and print(f'  Calculating distance from grid nodes to the surface using {self.workers} workers ...')
        dist_idx = array(KDTree(mesh_cells).query(grid_cells, workers=self.workers))

        # Inside or outside surface? 
        self.echo and print(f'  Creating unstructured grid with outer wall of {self.wall} voxels ...')
        idx = dist_idx[1].astype(int) # Mesh-surface index
        norm_dot_vec = mesh_normals[idx] * (mesh_cells[idx]-grid_cells)
        inside = npsum(norm_dot_vec, axis=1) > 0
        sign = -ones(inside.shape)
        sign[inside] = 1  # Positive values inside surface 

        # Add distance from grid node to surface as cell_data
        grid[name] = sign * dist_idx[0]

        # Add boundary markers
        grid[marker] = mesh_cell_id[idx]
        grid[marker][grid[name]>0] = 0  # Mark interior (fluid) nodes as boundary 0
        self.uniform = grid

        # The fluid nodes is confined by a wall of boundary nodes. 
        # Threshold is used to remove obsolete nodes
        # Threshold returns an unstructured grid
        self.grid = grid.threshold(value=-self.wall*dx, scalars=name)
        return self.grid


    #--------------------------------------------------------------------------------
    def add_tubes(self, name):
    #--------------------------------------------------------------------------------
        self.echo and print('  Creating tubes ...')
        cells = self.grid.cell_centers().points.astype(npdouble)
        tubes = ones((self.centerline.n_cells, self.grid.n_cells), dtype=int)
        for n in range(self.centerline.n_cells):
            cl_points = self.centerline.cell_points(n).astype(npdouble)
            dist = array(KDTree(cl_points).query(cells, workers=1))[0]
            m = mean(dist)
            remove = dist>m
            tubes[n,:] = n
            tubes[n,remove] = -1

        # Join tubes and select slices 
        tube = -1*ones(self.grid.n_cells, dtype=int)
        c = 0
        for tb in tubes[self.tube_list]:
            mask = (tb>=0)*(tube<0)
            tube[mask] = tb[mask]
        self.grid[name] = tube
        return tube


    #--------------------------------------------------------------------------------
    def reconnect_clusters(self, values=None, scalar=None):
    #--------------------------------------------------------------------------------
        # Check if tubes are unconnected
        #echo and print(f'  Check for unconnected {scalar} clusters ...')
        split_clusters = []
        solid_clusters = []
        for n in values:
            cluster = Cluster(self.grid.threshold(value=(n,n), scalars=scalar), n=n)
            if cluster.is_solid():
                solid_clusters.append(cluster)
            else:
                split_clusters.append(cluster)

        # Join unconnected clusters with solid clusters
        for split in split_clusters:
            joined = [(Cluster(split.grid + solid.grid), solid.n) for solid in solid_clusters]
            joined = [j for j in joined if j[0].num_bodies()==2]
            for j,n in joined:
                self.echo and print(f'  Unconnected {scalar} cluster {split.n} is merged into solid {scalar} cluster {n} ...')
                bodies = j.bodies()
                mrg_ind = [i for i,b in enumerate(bodies) if b[scalar].max()-b[scalar].min()!=0][0]
                merged = bodies[mrg_ind]
                solid_rest = bodies[int(not mrg_ind)]
                merged[scalar] = n * ones(merged.n_cells, dtype=int)            
                ind = [i for i,s in enumerate(solid_clusters) if s.n==n][0]
                solid_clusters[ind] = Cluster(merged, n=n)
                solid_clusters.append(Cluster(solid_rest, n=split.n))
                break
        self.solid[scalar] = solid_clusters
        return solid_clusters


    #--------------------------------------------------------------------------------
    def add_slices(self, name, tube=None):
    #--------------------------------------------------------------------------------
        # Make slices
        self.echo and print('  Calculating slice indices ...')
        c = 0
        tubes = self.solid[tube]
        for tb in tubes:
            dist, idx = tb.distance_to_centerline(self.centerline, workers=1)
            tb.grid[name] = 1 + c + idx - npmin(idx)
            c = npmax(tb.grid[name])
            #print(f'min, max, c : {npmin(idx)}, {npmax(idx)}, {c}')
        self.grid = UnstructuredGrid().merge([tb.grid for tb in tubes])


    #--------------------------------------------------------------------------------
    def add_processes(self, name, start=0, slice=None):
    #--------------------------------------------------------------------------------
        # Merge slices into processes
        self.echo and print(f'  Split grid in {self.nproc} process groups ...')
        slices = self.grid[slice]
        size = len(slices)
        proc_size = ceil(size/self.nproc)
        proc = -1*ones(size, dtype=int)
        stop = npmax(slices)
        a = b = 0
        p = start
        while p < self.nproc + start:
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
        self.grid[name] = proc
        return proc



#####################################  M A I N  ########################################

if __name__ == '__main__':
    import sys
    
    connected = True
    nproc = 10
    tube_list = [2, 3, 5, 1]
    echo = True

    narg = len(sys.argv) - 1
    if narg < 3:
        raise SystemExit('\n  Usage:\n    blood_grid <path to surface mesh file> <path to centerline file> <grid resolution> [number of processors]\n')
    meshpath = Path(sys.argv[1])
    clpath = Path(sys.argv[2])
    dx = float(sys.argv[3])
    workers = narg>3 and int(sys.argv[4]) or 1

    savename = meshpath.parent/'blood_grid.vtk'

    grid = Grid(dx=dx, surface=meshpath, centerline=clpath, tube_list=tube_list, nproc=nproc, echo=echo, connected=connected)
    grid.create()
    grid.save(savename)
    #ugrid = Grid(surface=meshpath, centerline=clpath)
    # ugrid.add_distance_map()

    # ugrid.add_tubes(select=tube_list)
    # ugrid.reconnect_clusters(values=tube_list, scalar='tube')
    # ugrid.add_slices()

    # ugrid.add_process_groups(nproc)

    # if group_processes: 
    #     echo and print('  Group split process groups ...')
    #     ugrid.reconnect_clusters(values=list(range(nproc)), scalar='proc')
    #     ugrid.grid = UnstructuredGrid().merge([solid.grid for solid in ugrid.solid['proc']], merge_points=True)

    # # Save final grid
    # ugrid.save(savename)
    # print(f'  Saved {savename}')

    #print(f'ugrid.n_cells: {ugrid.n_cells}, ugrid2.n_cells: {ugrid2.n_cells}, ugrid3.n_cells: {ugrid3.n_cells}')



