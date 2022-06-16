#!/usr/bin/env python3

from pyvista import read as pvread, UniformGrid, UnstructuredGrid
from numpy import argsort, array as nparray, ceil, min as npmin, max as npmax, sum as npsum, ones, double as npdouble, mean, zeros, sum as npsum
from scipy.spatial import KDTree
from pathlib import Path
import matplotlib.pyplot as pl

PROC_NAME = 'proc'
SLICE_NAME = 'slice'
TUBE_NAME = 'tube'

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
        dist_idx = nparray(KDTree(cl_points).query(self.grid.cell_centers().points, workers=workers))
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
                merged[self.scalar] = n * ones(merged.n_cells, dtype=int)            
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
    def __init__(self, dx=None, surface=None, centerline=None, wall=2, echo=True, 
                 tube_list=None, workers=1, nproc=1, connected=False) -> None:
    #--------------------------------------------------------------------------------
        self.proc_map = Process_map(nproc=nproc, connected=connected, echo=echo)
        self.echo = echo
        self.surface = Path(surface)
        self.dx = dx
        self.wall = wall
        self.tube_list = tube_list
        self.centerline = centerline
        self.workers = workers
        self.grid = None
        self.solid = {}
        self.echo and print(f'  Reading surface mesh file {surface} ... ', end='', flush=True)
        self.mesh = pvread(surface).extract_surface()
        self.echo and print('done')        
        bnds = self.mesh.bounds
        self.min, self.max = nparray(bnds[::2]), nparray(bnds[1::2])
        #xmin, xmax, ymin, ymax, zmin, zmax = self.mesh.bounds
        self.size = self.max - self.min # nparray((xmax-xmin, ymax-ymin, zmax-zmin))
        self.size += 2*self.wall*self.dx
        self.dim = ceil(self.size/self.dx).astype(int)
        # Filename
        dx_str = f'{int(dx)}-{round((dx%1)*100):2d}'
        dim_str = f'{"x".join([str(d) for d in self.dim])}'
        self.file = Path(f'{self.surface.parent/self.surface.stem}__LB__dx_{dx_str}__dim_{dim_str}.vtk')


    #--------------------------------------------------------------------------------
    def voxelize(self): 
    #--------------------------------------------------------------------------------
        self.centerline = pvread(self.centerline)
        self.add_distance_map('distance', marker='boundary', cell_id='CellEntityIds')
        self.add_tubes(TUBE_NAME)
        tubes = Clusters(self.grid, values=self.tube_list, scalar=TUBE_NAME).reconnect()
        self.add_slices(SLICE_NAME, tubes=tubes)
        return self


    #--------------------------------------------------------------------------------
    def create(self, save=True): 
    #--------------------------------------------------------------------------------
        if self.file.is_file():
            self.read_grid()
        else:
            self.voxelize()
            self.save_grid()
        self.grid = self.proc_map.distribute_load(self.grid)
        if save:
            self.save_grid(name='blood_geo.vtk')


    #--------------------------------------------------------------------------------
    def save_grid(self, name=None, echo=True, **kwargs):
    #--------------------------------------------------------------------------------
        if not name:
            name = self.file
        echo and print(f'  Saving grid as {name} ... ', end='', flush=True)
        self.grid.save(Path(name).with_suffix('.vtk'), **kwargs)
        echo and print('done')


    #--------------------------------------------------------------------------------
    def read_grid(self, name=None, echo=True):
    #--------------------------------------------------------------------------------
        if not name:
            name = self.file
        path = Path(name)
        if path.is_file():
            echo and print(f'  Reading grid file {path} ... ', end='', flush=True)
            self.grid = pvread(path)
            echo and print('done')
        else:
            raise FileNotFoundError(f'{path}')

    #--------------------------------------------------------------------------------
    def array(self, value):
    #--------------------------------------------------------------------------------
        arr = zeros(self.dim)
        if isinstance(value, str):
            data = self.grid[value]
        else:
            # Assume number...
            data = value*ones(self.grid.n_cells)
        ijk = self.grid['IJK']
        arr[ijk[:,0], ijk[:,1], ijk[:,2]] = data
        return arr

    #--------------------------------------------------------------------------------
    def add_distance_map(self, name, marker=None, cell_id=None):
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
        grid_cells = grid.cell_centers().points.astype(npdouble)
        ijk = ( (grid_cells-grid.origin)/dx ).astype(int)
        grid['IJK'] = ijk
        #grid['BLOCK_I'] = ijk[:,0]
        #grid['BLOCK_J'] = ijk[:,1]
        #grid['BLOCK_K'] = ijk[:,2]

        # Return shortest distance from grid-cell to mesh surface cell together with mesh index
        bndry_cells = id_scalar > 0
        mesh_cells   = self.mesh.cell_centers().points.astype(npdouble)[bndry_cells]
        mesh_normals = self.mesh.cell_normals[bndry_cells]
        mesh_cell_id = self.mesh[cell_id][bndry_cells]
        self.echo and print(f'  Calculating distance from grid nodes to the surface using {self.workers} workers ...')
        dist_idx = nparray(KDTree(mesh_cells).query(grid_cells, workers=self.workers))

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
            dist = nparray(KDTree(cl_points).query(cells, workers=1))[0]
            m = mean(dist)
            remove = dist>m
            tubes[n,:] = n
            tubes[n,remove] = -1

        # Join tubes 
        tube = -1*ones(self.grid.n_cells, dtype=int)
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
            tb.grid[name] = 1 + c + idx - npmin(idx)
            c = npmax(tb.grid[name])
            #print(f'min, max, c : {npmin(idx)}, {npmax(idx)}, {c}')
        self.grid = UnstructuredGrid().merge([tb.grid for tb in tubes])


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
        self.volume = zeros(nproc)

    #--------------------------------------------------------------------------------
    def distribute_load(self, grid, stat=True): 
    #--------------------------------------------------------------------------------
        self.echo and print(f'  Creating process-map with {self.nproc} processes ...')
        self.grid = grid
        self.add_process_array()
        stat and self.statistics()
        if self.connected: 
            conn = Clusters(self.grid, values=list(range(1,self.nproc+1)), scalar=PROC_NAME).reconnect()
            self.grid = UnstructuredGrid().merge([c.grid for c in conn], merge_points=True)
        return self.grid

    #--------------------------------------------------------------------------------
    def statistics(self):
    #--------------------------------------------------------------------------------
        print()
        print('    Process |  Volume  |  Deviation ')
        print('    --------------------------------')
        tot = npsum(self.volume)
        share = 1/self.nproc
        for i in range(self.nproc):
            print(f'    {i+1: 7d} |  {self.volume[i]:.1e} |  {100*(self.volume[i]/tot - share)/share: .1f} %' )
        print('    --------------------------------')
        print()


    #--------------------------------------------------------------------------------
    def add_process_array(self):
    #--------------------------------------------------------------------------------
        # Merge slices into processes
        slices = self.grid[SLICE_NAME]
        size = len(slices)
        proc_size = ceil(size/self.nproc)
        proc = -1*ones(size, dtype=int)
        stop = npmax(slices)
        a = b = 0
        p = self.start
        while p < self.nproc + self.start:
            mask = (slices >= a) * (slices <= b)
            maskb = (slices >= a) * (slices <= b+1)
            if abs(proc_size-npsum(mask)) < abs(proc_size-npsum(maskb)) or b == stop:
                #print(f'p: {p}, a,b: {a},{b}  {npsum(mask)}, {proc_size}, {npsum(mask)-proc_size}')
                a = b + 1 
                proc[mask] = p
                self.volume[p-self.start] = npsum(mask)
                p += 1
            b += 1
        if npsum(proc==-1) > 0:
            print(f'  WARNING! {npsum(proc==-1)} nodes not assigned to a process!')
        self.grid[PROC_NAME] = proc
        


