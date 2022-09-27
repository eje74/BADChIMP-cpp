#!/usr/bin/env python3

from pyvista import read as pvread, UniformGrid, UnstructuredGrid
#from numpy import newaxis, sqrt, argsort, array as nparray, ceil, min as npmin, max as npmax, sum as npsum, ones, double as npdouble, mean, zeros, sum as npsum
import numpy as np
from scipy.spatial import KDTree
from pathlib import Path
from copy import deepcopy
#import matplotlib.pyplot as pl

PROC_NAME = 'proc'
SLICE_NAME = 'slice'
TUBE_NAME = 'tube'


#--------------------------------------------------------------------------------
def distance_map(surface_points, grid_points, normals=None, workers=2, **kwargs):
#--------------------------------------------------------------------------------
    ### dist and idx are the same length as grid_points
    ### dist[0] : Distance from grid_points[0] to mesh_points[idx[0]]
    dist, idx = KDTree(surface_points).query(grid_points, workers=workers, **kwargs)
    if normals is None:
        return dist, idx
    norm_dot_vec = normals[idx] * (surface_points[idx]-grid_points)
    inside = np.sum(norm_dot_vec, axis=1) > 0
    sign = -np.ones(inside.shape)
    sign[inside] = 1  # Positive values inside surface 
    return sign * dist, idx


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
    def n_cells(self):
    #--------------------------------------------------------------------------------
        return [body.n_cells for body in self.bodies()]

    #--------------------------------------------------------------------------------
    def bodies_by_volume(self):
    #--------------------------------------------------------------------------------
        vol = self.volume()
        keys = self.body.keys()
        return [self.body[keys[i]] for i in np.argsort(vol)]

    #--------------------------------------------------------------------------------
    def merge(self, bodies):
    #--------------------------------------------------------------------------------
        self.grid = self.grid.merge(bodies, merge_points=True)
        self.body = self.grid.connectivity().split_bodies()

    #--------------------------------------------------------------------------------
    def distance_to_centerline(self, centerline, workers=1): 
    #--------------------------------------------------------------------------------
        cl_points = centerline.cell_points(self.n).astype(np.double)
        dist_idx = np.array(KDTree(cl_points).query(self.grid.cell_centers().points, workers=workers))
        return dist_idx[0], dist_idx[1].astype(int)

#    def cell_centers(self):
#        return self.grid.cell_centers().points.astype(npdouble)

#====================================================================================
class Clusters:
#====================================================================================
    #--------------------------------------------------------------------------------
    def __init__(self, grid, values=None, scalar=None) -> None:
    #--------------------------------------------------------------------------------
        self.scalar = scalar
        self.split = []
        self.solid = []
        for n in values:
            cluster = Cluster(grid.threshold(value=(n,n), scalars=scalar), n=n)
            if cluster.is_solid():
                self.solid.append(cluster)
            else:
                self.split.append(cluster)


    #--------------------------------------------------------------------------------
    def get(self):
    #--------------------------------------------------------------------------------
        return self.solid + self.split


    #--------------------------------------------------------------------------------
    def reconnect(self, echo=True):
    #--------------------------------------------------------------------------------
        # Join unconnected clusters with solid clusters
        for split in self.split:
            cluster_list = [(Cluster(split.grid + solid.grid), solid.n) for solid in self.solid]
            joined = [c for c in cluster_list if c[0].num_bodies()==2]
            if not joined:
                print(f'  WARNING! Unable to merge split {self.scalar} cluster {split.n}')
                return self.solid + self.split
            for j,n in joined:
                echo and print(f'  Unconnected {self.scalar} cluster {split.n} is merged into solid {self.scalar} cluster {n} ...')
                bodies = j.bodies()
                mrg_ind = [i for i,b in enumerate(bodies) if b[self.scalar].max()-b[self.scalar].min()!=0][0]
                merged = bodies[mrg_ind]
                merged[self.scalar] = n * np.ones(merged.n_cells, dtype=int)            
                # Merge connected part of split cluster to solid n 
                ind = [i for i,s in enumerate(self.solid) if s.n==n][0]
                self.solid[ind] = Cluster(merged, n=n)
                # Rest of split cluster
                bodies.pop(mrg_ind) 
                if len(bodies) == 1:
                    # Create new solid cluster and stop 
                    self.solid.append(Cluster(bodies[0], n=split.n))
                    break
        return self.solid


#====================================================================================
class Geo:
#====================================================================================
    #--------------------------------------------------------------------------------
    def __init__(self, dx=None, surface=None, centerline=None, wall=3, echo=True, 
                 tube_list=None, workers=1, nproc=1, connected=False) -> None:
    #--------------------------------------------------------------------------------
        self.proc_map = Process_map(nproc=nproc, connected=connected, echo=echo)
        self.echo = echo
        self.surface = None
        self.dx = dx 
        self.wall = wall
        self.tube_list = tube_list
        self.centerline = centerline
        self.workers = workers
        self.grid = None
        self.solid = {}
        self.mesh = None
        self.min, self.max = 0, 0
        self.size = 0
        self.dim = 0
        self.file = ''
        if surface:
            self.surface = Path(surface).resolve()
            self.echo and print(f'  Reading surface mesh file {self.surface.relative_to(Path.cwd())} ... ', end='', flush=True)
            self.mesh = pvread(surface).extract_surface()
            self.echo and print('done')        
            self.min, self.max = np.hsplit(np.array(self.mesh.bounds).reshape((3,2)), 2)
            self.size = (self.max - self.min).squeeze() 
            self.size += 2*self.wall*self.dx
            self.dim = np.ceil(self.size/self.dx).astype(int)
            # Filename
            dx_str = f'{int(dx)}-{round((dx%1)*100):2d}'
            dim_str = f'{"x".join([str(d) for d in self.dim])}'
            self.file = Path(f'{self.surface.parent/self.surface.stem}__LB__dx_{dx_str}__dim_{dim_str}.vtk')

    #--------------------------------------------------------------------------------
    def copy(self):
    #--------------------------------------------------------------------------------
        geo = Geo()
        geo.dim = deepcopy(self.dim)
        geo.grid = deepcopy(self.grid)
        return geo

    #--------------------------------------------------------------------------------
    def threshold(self, *args, **kwargs):
    #--------------------------------------------------------------------------------
        self.grid = self.grid.threshold(*args, **kwargs)
        return self

    #--------------------------------------------------------------------------------
    def voxelize(self): 
    #--------------------------------------------------------------------------------
        self.centerline = pvread(self.centerline)
        self.add_distance_map('distance', boundary='boundary', cell_id='CellEntityIds')
        self.add_tubes(TUBE_NAME)
        tubes = Clusters(self.grid, values=self.tube_list, scalar=TUBE_NAME).reconnect()
        self.add_slices(SLICE_NAME, tubes=tubes)
        return self


    #--------------------------------------------------------------------------------
    def create(self, use_file=True): #, name, save=True): 
    #--------------------------------------------------------------------------------
        if use_file and self.file.is_file():
            self.read_grid()
        else:
            self.voxelize()
            self.save_grid()
        self.grid = self.proc_map.distribute_load(self.grid)
        ### Only fluid and boundary nodes have a process number.
        ### We need to set the wall nodes to 0
        self.grid[PROC_NAME][self.grid['wall']] = 0
        # if save:
        #     self.save_grid(name=name)

    #--------------------------------------------------------------------------------
    def save(self, name, **kwargs):
    #--------------------------------------------------------------------------------
        self.save_grid(name=name, **kwargs)

    #--------------------------------------------------------------------------------
    def save_grid(self, name=None, echo=True, **kwargs):
    #--------------------------------------------------------------------------------
        if not name:
            name = self.file
        name = Path(name).with_suffix('.vtk').resolve()
        echo and print(f'  Saving grid as {name.relative_to(Path.cwd())} ... ', end='', flush=True)
        self.grid.save(name, **kwargs)
        echo and print('done')


    #--------------------------------------------------------------------------------
    def read_grid(self, name=None, echo=True):
    #--------------------------------------------------------------------------------
        if not name:
            name = self.file
        path = Path(name).resolve()
        if path.is_file():
            echo and print(f'  Reading grid file {path.relative_to(Path.cwd())} ... ', end='', flush=True)
            self.grid = pvread(path)
            echo and print('done')
        else:
            raise FileNotFoundError(f'{path}')

    #--------------------------------------------------------------------------------
    def array(self, value):
    #--------------------------------------------------------------------------------
        # arr = zeros(self.dim)
        if isinstance(value, str):
            data = self.grid[value]
        else:
            # Assume number...
            data = value*np.ones(self.grid.n_cells)
        ijk = self.grid['IJK']
        data_dim = int(data.size/data.shape[0])
        arr = np.zeros(list(self.dim) + [data_dim]).squeeze()
        #print(f'value: {value}, shape: {arr.shape}')
        arr[ijk[:,0], ijk[:,1], ijk[:,2]] = data
        return arr

    #--------------------------------------------------------------------------------
    def add_distance_map(self, distance, boundary=None, cell_id=None):
    #--------------------------------------------------------------------------------
        # Only use boundary cells, i.e. cells with 'CellEntityIds' > 0 (created by VMTK)
        # cell_id = 'CellEntityIds' 
        id_scalar = self.mesh.cell_data.get(cell_id)
        if id_scalar is None:
            raise SystemExit(f'  ERROR! The mesh-file is missing the {cell_id} scalar that identify boundary cells, aborting...\n')

        # Create uniform grid
        self.echo and print(f'  Creating uniform grid of dimension {self.dim} ...')
        grid = UniformGrid()
        grid.dimensions = self.dim + 1 # Because we want to add cell data (not point data)
        dx = self.dx
        #grid.origin = nparray((xmin, ymin, zmin))-self.wall*dx
        grid.origin = self.min - self.wall*dx
        grid.spacing = (dx, dx, dx) # Evenly spaced grids
        grid_cells = grid.cell_centers().points.astype(np.double)
        ijk = ( (grid_cells-grid.origin)/dx ).astype(int)
        grid['IJK'] = ijk
        #grid['BLOCK_I'] = ijk[:,0]
        #grid['BLOCK_J'] = ijk[:,1]
        #grid['BLOCK_K'] = ijk[:,2]

        # Return shortest distance from grid-cell to mesh surface cell together with mesh index
        mesh_cells = self.mesh.cell_centers().points.astype(np.double)
        mesh_normals = self.mesh.cell_normals
        mesh_points = np.array([self.mesh.cell_points(i) for i in range(self.mesh.n_cells)])
        mesh_id = self.mesh[cell_id]
        bndry = id_scalar > 0
        # bndry_mesh_cells   = mesh_cells[bndry]
        # bndry_mesh_normals = mesh_normals[bndry]
        # bndry_mesh_cell_id = mesh_id[bndry]
        
        # self.echo and print(f'  Calculating distance from grid nodes to the surface using {self.workers} workers ...')
        # k = 1 # k=1, return only nearest neighbour
        # dist, idx = KDTree(bndry_mesh_cells).query(grid_cells, k=k, workers=self.workers)

        # ### Inside or outside surface? 
        # self.echo and print(f'  Creating unstructured grid with outer wall of {self.wall} voxels ...')
        # norm_dot_vec = bndry_mesh_normals[idx] * (bndry_mesh_cells[idx]-grid_cells)
        # inside = npsum(norm_dot_vec, axis=1) > 0
        # sign = -ones(inside.shape)
        # sign[inside] = 1  # Positive values inside surface 

        self.echo and print(f'  Calculating distance from grid nodes to the surface using {self.workers} workers ...')
        dist, idx = distance_map(mesh_cells[bndry], grid_cells, normals=mesh_normals[bndry], workers=self.workers)
        mesh_id_bndry = mesh_id[bndry][idx]

        ### Add distance from grid node to surface as cell_data
        grid[distance] = dist/dx
        outside = grid[distance] < 0
        # grid['wall'] = (mesh_id_bndry==1) * outside 
        grid['wall'] = np.array(outside, dtype=int) 

        ### Add boundary markers
        P_bndry = mesh_id_bndry > 1
        fluid_at_surface = (grid[distance]>=0) * (grid[distance]<=np.sqrt(2))
        ### cell_id for in/out boundaries starts at 2, so we subtract 1 to make it start at 1. 
        grid[boundary] = (mesh_id_bndry-1) * P_bndry * fluid_at_surface
        ### The boundary caps are identified by -1
        grid[boundary] -= P_bndry * outside * (grid[distance]>=-np.sqrt(2))
        grid[boundary] = grid[boundary].astype(int)
        ###
        #grid['wall'][grid[boundary]<0] = 1

        ### Add normal vectors for all nodes touching the surface/wall
        grid['normal'] = mesh_normals[bndry][idx] * (P_bndry*fluid_at_surface)[:, None]
        self.uniform = grid

        ### The fluid nodes is confined by a wall of boundary nodes. 
        ### Threshold is used to remove obsolete nodes
        ### Threshold returns an unstructured grid
        self.grid = grid.threshold(value=-self.wall, scalars=distance)

        ### Use vertex-points (3) instead of cell-centers to 
        ### improve the location of the pressure boundary
        geo_grid_cells = self.grid.cell_centers().points.astype(np.double)
        ### Array to map from vertex-point to cell number
        point_to_cell = np.repeat(np.arange(self.mesh.n_cells), 3)
        ### Loop over all pressure boundaries
        for i in range(2, np.max(mesh_id)+1):
            id_mask = mesh_id == i
            points = mesh_points[id_mask].reshape((-1,3)) # points must be a 2D-matrix
            dist, ind = distance_map(points, geo_grid_cells)
            grid_mask = (dist<=np.sqrt(3)*dx) * (self.grid['wall']!=1) * (self.grid[boundary]>=0)
            mesh_ind = point_to_cell[ind[grid_mask]]
            self.grid[boundary][grid_mask] = mesh_id[id_mask][mesh_ind]-1
            self.grid['normal'][grid_mask] = mesh_normals[id_mask][mesh_ind]

        return self.grid

    #--------------------------------------------------------------------------------
    def add_tubes(self, name):
    #--------------------------------------------------------------------------------
        self.echo and print('  Creating tubes ...')
        cells = self.grid.cell_centers().points.astype(np.double)
        tubes = np.ones((self.centerline.n_cells, self.grid.n_cells), dtype=int)
        for n in range(self.centerline.n_cells):
            cl_points = self.centerline.cell_points(n).astype(np.double)
            dist = np.array(KDTree(cl_points).query(cells, workers=1))[0]
            m = np.mean(dist)
            remove = dist>m
            tubes[n,:] = n
            tubes[n,remove] = -1

        # Join tubes 
        tube = -1*np.ones(self.grid.n_cells, dtype=int)
        c = 0
        for tb in tubes[self.tube_list]:
            mask = (tb>=0)*(tube<0)
            tube[mask] = tb[mask]
        self.grid[name] = tube
        return tube


    #--------------------------------------------------------------------------------
    def add_slices(self, name, tubes=None):
    #--------------------------------------------------------------------------------
        # Make slices
        self.echo and print('  Calculating slice indices ...')
        c = 0
        #tubes = self.solid[tube]
        for tb in tubes:
            dist, idx = tb.distance_to_centerline(self.centerline, workers=1)
            tb.grid[name] = 1 + c + idx - np.min(idx)
            c = np.max(tb.grid[name])
            #print(f'min, max, c : {npmin(idx)}, {npmax(idx)}, {c}')
        self.grid = UnstructuredGrid().merge([tb.grid for tb in tubes])

    # #--------------------------------------------------------------------------------
    # def threshold(self, *args, **kwargs):
    # #--------------------------------------------------------------------------------
    #     new_geo = deepcopy(self)
    #     new_geo.grid = new_geo.grid.threshold(*args, **kwargs)
    #     return new_geo


#====================================================================================
class Process_map:
#====================================================================================
    #--------------------------------------------------------------------------------
    def __init__(self, nproc=1, start=1, connected=True, grid=None, echo=True) -> None:
    #--------------------------------------------------------------------------------
        self.nproc = nproc
        self.connected = connected
        self.grid = grid
        self.start = start
        self.echo = echo
        self.clusters = None

    #--------------------------------------------------------------------------------
    def distribute_load(self, grid, stat=True): 
    #--------------------------------------------------------------------------------
        self.echo and print(f'  Creating process-map with {self.nproc} processes ...')
        self.grid = grid
        values = list(range(1,self.nproc+1))
        self.add_process_array()
        # Split grid in clusters based on process number
        self.clusters = Clusters(self.grid, values=values, scalar=PROC_NAME)
        if self.connected: 
            solid = self.clusters.reconnect()
            self.grid = UnstructuredGrid().merge([s.grid for s in solid], merge_points=True)
            self.clusters = Clusters(self.grid, values=values, scalar=PROC_NAME)
        stat and self.statistics()
        return self.grid

    #--------------------------------------------------------------------------------
    def statistics(self):
    #--------------------------------------------------------------------------------
        print()
        print('    Process |   Size   |  Deviation  |  Parts')
        print('    ------------------------------------------')
        size = [c.n_cells() for c in self.clusters.get()]
        parts = [len(s) for s in size]
        totsize = np.array([sum(s) for s in size])
        share = 1/self.nproc
        dev = 100*(totsize/np.sum(totsize) - share)/share
        for i, v in enumerate(totsize): #range(self.nproc): 
            txt = f'    {i+1: 7d} |  {v:.1e} |    {dev[i]:5.1f} %  |  {parts[i]: 2d}'
            if parts[i] > 1:
                txt += ' (' + ', '.join([f'{100*s/v:.0f}%' for s in size[i]]) + ')'
            print(txt)
        print('    ------------------------------------------')
        print()


    #--------------------------------------------------------------------------------
    def add_process_array(self):
    #--------------------------------------------------------------------------------
        # Merge slices into processes
        slices = self.grid[SLICE_NAME]
        size = len(slices)
        proc_size = np.ceil(size/self.nproc)
        proc = -1*np.ones(size, dtype=int)
        stop = np.max(slices)
        a = b = 0
        p = self.start
        while p < self.nproc + self.start:
            mask = (slices >= a) * (slices <= b)
            maskb = (slices >= a) * (slices <= b+1)
            if abs(proc_size-np.sum(mask)) < abs(proc_size-np.sum(maskb)) or b == stop:
                #print(f'p: {p}, a,b: {a},{b}  {npsum(mask)}, {proc_size}, {npsum(mask)-proc_size}')
                a = b + 1 
                proc[mask] = p
                p += 1
            b += 1
        if np.sum(proc==-1) > 0:
            print(f'  WARNING! {np.sum(proc==-1)} nodes not assigned to a process!')
        self.grid[PROC_NAME] = proc
        


