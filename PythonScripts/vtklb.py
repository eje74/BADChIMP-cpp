#from numba import jit
import numpy as np


#@jit # Set "nopython" mode for best performance, equivalent to @njit    
def write_neighbors_fast(self_ind_local, self_basis, self_local_label, self_geo, self_label,self_file, rank):
    # Neighbor header
    self_file.write( "NEIGHBORS int\n" )
    # list of neighbors for each node
    num_points_global = np.shape(self_ind_local)[1]
    num_dim = np.shape(self_ind_local)[0]
    num_basis = np.shape(self_basis)[0]
    sys_size = self_local_label.shape
    
    # Strategy to handle big systems. Devide the sytem into smaller boxes
    # -- box size
    box_size = int(np.ceil(num_points_global/(num_basis)))
    beginBox = 0
    endBox = min(beginBox + box_size, num_points_global)
    num_points = endBox - beginBox
    while beginBox < num_points_global: 
        neig_list = np.zeros((num_basis, num_points), dtype=int)
        for cnt,  c in enumerate(self_basis):
            pos_tmp = np.array(np.arange(num_points), dtype=int)
            ind_neig = np.array(self_ind_local)[:, beginBox:endBox] + c.reshape((num_dim, 1))
            # Remove neighbor points outside the domain
            for dim in np.arange(num_dim):
                pos_tmp = pos_tmp[ind_neig[dim, :]>= 0]
                ind_neig = ind_neig[:, ind_neig[dim, :]>= 0]
                pos_tmp = pos_tmp[ sys_size[dim] > ind_neig[dim, :]]
                ind_neig = ind_neig[:, sys_size[dim] > ind_neig[dim, :]]
            # Update neiglist
            neig_list[cnt, pos_tmp] = self_local_label[tuple(ind_neig)]
            # Add periodic boundaries on same process
            pos_tmp = pos_tmp[self_local_label[tuple(ind_neig)] == 0]
            ind_neig = ind_neig[:, self_local_label[tuple(ind_neig)] == 0]
            pos_tmp = pos_tmp[self_geo[tuple(ind_neig)]==rank]
            ind_neig = ind_neig[:, self_geo[tuple(ind_neig)]==rank]
            neig_list[cnt, pos_tmp] = self_label[tuple(ind_neig)]
        # Write to file
        np.savetxt(self_file, neig_list.transpose(), fmt="%d", delimiter=' ', newline='\n')    
        beginBox = endBox  
        endBox = min(beginBox + box_size, num_points_global)
        num_points = endBox - beginBox



class vtklb:
    # The list of indecies follows the labeling so that
    # the node with label 1 is the first in the point list.
    #  The fluid nodes on a given rank are the first one labled
    # then the solid wall and finally the nodes on different 
    # processors
    #
    # Row major Column major?
    #  Here we use the convention that the first index gives
    #  the x position, the second gives the y-direction and  
    #  the third gives the z-direction.  
    #  Hence, A[x,y,z]
    #  For Vector and tensor we have that
    #  * v[0] = v_x, v[1] = v_y, ...
    #  * S[0,0] = S_xx, S[1,0] = S_yx, ...
    #  
    #
    # TODO:
    # - Make a sparse version of this class
    # - Check how to use numpy write to file functions
    # - Add periodic boundaries
    #
    def __init__(self, geo, basis, periodic="", name='tmp', path='~/', version='na'):
        # geo :  ndarray with int values from 0 and up
        # basis : Either ndarray or string
        # periodic: string that indicate periodicity by x, y and z.
        #           If the system is periodic in x and z write periodic = "xz"
        # Number of spatial dimensions
        self.nd = geo.ndim  
        # Sytem size
        self.n = geo.shape
        # number of processes
        self.np = np.amax(geo)
        # rim size
        self.rim_width = 1
        # Add a ghost node rim to the global geometry
        self.geo = np.zeros(tuple(n+2*self.rim_width for n in self.n), dtype=int)
        self.geo[:] = -1  # Using -1 as an marker for unset values
        self.geo[tuple([slice(self.rim_width,-self.rim_width)]*self.nd)] = geo
        # Boolean array that define the bulk of the system
        #  helps to added information at the ghost nodes
        self.bulk = np.zeros(self.geo.shape, dtype=bool)
        self.bulk[tuple([slice(self.rim_width,-self.rim_width)]*self.nd)] = True
        # Setup the global node numbers
        self.label = np.zeros(self.geo.shape, dtype=int)
        for rank in np.arange(1,np.amax(geo)+1):
            self.label[(self.geo == rank) & self.bulk] = np.arange(1, 1 + np.count_nonzero((self.geo == rank) & self.bulk))   
        # Setup periodic boundary conditioins
        self.periodic = periodic.lower()
        if 'x' in periodic.lower():
            print("PERIODIC IN X\n");
            self.geo[0, ...] = self.geo[-2, ...]
            self.geo[-1, ...] = self.geo[1, ...]
            self.label[0, ...] = self.label[-2, ...]
            self.label[-1, ...] = self.label[1, ...]
        if 'y' in periodic.lower():
            print("PERIODIC IN Y\n");
            self.geo[:,0, ...] = self.geo[:,-2, ...]
            self.geo[:,-1, ...] = self.geo[:,1, ...]
            self.label[:, 0, ...] = self.label[:,-2, ...]
            self.label[:, -1, ...] = self.label[:,1, ...]
        if 'z' in periodic.lower():
            print("PERIODIC IN Z\n");
            self.geo[:,:,0, ...] = self.geo[:,:,-2, ...]
            self.geo[:,:,-1, ...] = self.geo[:,:,1, ...]
            self.label[:,:,0, ...] = self.label[:,:,-2, ...]
            self.label[:,:,-1, ...] = self.label[:,:,1, ...]                
        # Boolean array used to mark nodes as added to a list
        self.added = np.zeros(self.geo.shape, dtype=bool)    
        # Preamable information_     
        self.version = version
        # File information
        self.filename = name
        self.path = path
        self.basis = self.get_basis(basis)
        self.write()
        

    def set_path(self, path):
        self.path = path

                
    def set_filename(self, name):
        self.filename = name


    def open(self, rank):
        self.local_filename = self.filename + str(rank-1) + ".vtklb"
        try:
            self.file = open(self.path+self.local_filename, "w")
        except IOError:
            #print("Could not open file: {}".format(self.path+self.filename)) 
            raise SystemExit("ERROR! Could not open file: {}\n".format(self.path+self.filename))


    def append(self, rank):        
        self.local_filename = self.filename + str(rank) + ".vtklb"
        try:
            self.file = open(self.path+self.local_filename, "a")
        except IOError:
            print("Could not open file: {}".format(self.path+self.filename)) 

        
    def read(self, rank):
        # Opens file for reading
        self.local_filename = self.filename + str(rank) + ".vtklb"
        try:
            self.file = open(self.path+self.local_filename, "r")
        except IOError:
            print("Could not open file: {}".format(self.path+self.filename)) 

                                                
    def close(self):
        try:
            self.file.close()
        except IOError:
            print("Error when closing file {}".format(self.path+self.filename))

                        
    def get_basis(self, basis):  
        if type(basis) == str:
            if basis == "D2Q9":
                print("SYSTEM DEFINED BASIS D2Q9\n")                
                return np.array([[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1], [0, 0]], dtype=int)
            if basis == "D3Q19":
                print("SYSTEM DEFINED BASIS D3Q19\n")
                return np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], 
                                 [0, 1, 1], [0, 1, -1], [-1, 0, 0], [0, -1, 0], [0, 0, -1], [-1, -1, 0], [-1, 1, 0], 
                                 [-1, 0, -1], [-1, 0, 1], [0, -1, -1], [0, -1, 1], [0, 0, 0]], dtype=int)   
        if type(basis) == np.ndarray:
            print("USER DEFINED BASIS OF TYPE D{}Q{}\n".format(basis.shape[1], basis.shape[0]))            
            return np.copy(basis)

                
    def setup_processor_labels(self, rank):
        # Find the indecies and use them as position
        self.ind = np.where( (self.geo == rank) & self.bulk)        
        # Ensure that element positions folows the labeling
        ind_sort = np.argsort(self.label[self.ind])
        for axis_ind in self.ind:
            axis_ind[:] = axis_ind[ind_sort]
        # List of local indecies
        self.ind_local = tuple(np.copy(self.ind))        
        # Find solids        
        self.added[:] = False
        for v in self.basis:
            ind_tmp = tuple( ind + dind for dind, ind in zip(v, self.ind) )
            self.added[ind_tmp] = (self.geo[ind_tmp] == 0) | self.added[ind_tmp] 
        # Insert solid wall nodes into self.ind
        self.ind_local = tuple(np.concatenate((ind1, ind2)) for ind1, ind2 in zip(self.ind_local, np.where(self.added)))
        # Generate a list of neighboring processors
        self.rank_set = set()
        for v in self.basis:
            ind_tmp = tuple(ind + dind for dind, ind in zip(v, self.ind) )
            self.rank_set = self.rank_set.union(set(self.geo[ind_tmp]))
        #   remove ghost, solid and current rank from the set
        self.rank_set.discard(-1)
        self.rank_set.discard(0)
        self.rank_set.discard(rank)
        # Find nodes on neighboring processors
        for r in self.rank_set:
            # reset added matrix
            self.added[:] = False
            for v in self.basis:
                ind_tmp = tuple( ind + dind for dind, ind in zip(v, self.ind) )
                self.added[ind_tmp] = (self.geo[ind_tmp] == r) | self.added[ind_tmp]
            # Insert neighboring nodes to self.ind            
            self.ind_local = tuple(np.concatenate((ind1, ind2)) for ind1, ind2 in zip(self.ind_local, np.where(self.added)))
        # Setup a local label matrix    
        self.local_label = np.zeros(self.label.shape, dtype=int) 
        self.local_label[self.ind_local] = np.arange(1, 1 + len(self.ind_local[0]))   


                
    def write_preamble(self, rank):
        self.file.write( "# BADChIMP vtklb Version {}\n".format(self.version) )
        self.file.write( "Geometry file for process {}\n".format(rank-1) )      
        self.file.write( "ASCII\n")

                
    def write_dataset(self):
        self.file.write( "DATASET UNSTRUCTURED_LB_GRID\n" )
        self.file.write( "NUM_DIMENSIONS {}\n".format(self.geo.ndim) )
        line = "GLOBAL_DIMENSIONS";
        for nd in range(self.geo.ndim):
            line += " {}".format(self.geo.shape[nd])
        self.file.write(line + "\n");
        self.file.write( "USE_ZERO_GHOST_NODE\n" )

                
    def write_points(self):
        self.file.write( "POINTS {} int\n".format(len(self.ind_local[0])) )
        np.savetxt(self.file, (np.array(self.ind_local)-self.rim_width).transpose(), fmt="%d", delimiter=' ', newline='\n')    
        
                        
    def read_points(self, rank):
        self.read(rank)
        line = self.file.readline()
        words = line.split()
        # Read file until point pos block
        while words[0] != "POINTS":
            line = self.file.readline()
            words = line.split()
        num_points = int(words[1])        
        self.ind_local = np.zeros((self.nd, num_points),dtype=int)
        for n in range(num_points):
            line = self.file.readline()
            words = line.split()
            self.ind_local[:, n] = [int(w) for w in words]      
        self.ind_local = tuple(self.ind_local + self.rim_width)          
        self.close()
        
        
    def write_lattice(self):
        self.file.write( "LATTICE {} int\n".format(len(self.basis)) )
        np.savetxt(self.file, self.basis, fmt="%d", delimiter=' ', newline='\n')    
      

    def write_neighbors(self, rank):
        write_neighbors_fast(self.ind_local, self.basis, self.local_label, self.geo, self.label,self.file, rank)
                                                                                                                                                                    
                        
    def write_mpi(self, rank):
        # MPI header
        self.file.write( "PARALLEL_COMPUTING {}\n".format(rank-1) )
        # Go through the neighboring ranks
        for neig_rank in self.rank_set:
            # Find local labels overlapping with neighboring processes
            neig_ind = tuple( np.array(self.ind_local)[:, self.geo[self.ind_local] == neig_rank] )
            # Rank header
            self.file.write( "PROCESSOR {} {}\n".format(len(neig_ind[0]), neig_rank-1) )
            # Write list of overlapping nodes [local label, neigh process label]
            np.savetxt(self.file, np.array([self.local_label[neig_ind], self.label[neig_ind]]).transpose(), fmt="%d", delimiter=' ', newline='\n')            

        
    def write_point_data(self):
        # Header for point data
        self.file.write( "POINT_DATA {}\n".format(len(self.ind_local[0])) )

                
    def write_node_type(self):
        # Write the node type:
        # 0 : solid, that is geo = 0
        # 1 : fluid, that is geo > 0
        # dataset atribute header
        self.file.write( "SCALARS nodetype int\n" )
        np.savetxt(self.file, (self.geo[self.ind_local]>0).transpose().astype(int), fmt="%d", delimiter=' ', newline='\n')

                        
    def write_data_set_attribute(self, name, val):
#        np.issubdtype(some_dtype, np.integer)
        if np.issubdtype(val.dtype, np.integer):
            self.file.write( "SCALARS {} int\n".format(name) )
            np.savetxt(self.file, val[self.ind_local].transpose(), fmt="%d", delimiter=' ', newline='\n')
        else:
            self.file.write( "SCALARS {} float\n".format(name) )
            np.savetxt(self.file, val[self.ind_local].transpose(), delimiter=' ', newline='\n')
                 
             
    def write_proc(self, rank):
        self.open(rank)
        self.setup_processor_labels(rank)
        self.write_preamble(rank)
        self.write_dataset() 
        self.write_points() 
        self.write_lattice()
        self.write_neighbors(rank)
        self.write_mpi(rank)
        self.write_point_data()
        self.write_node_type()
        self.close()
        print( "\rWrote file: {}".format(self.local_filename), end='' )

                
    def write(self):
        print( "Writing files to folder {}".format(self.path))
        for rank in np.arange(self.np):
            self.write_proc(rank + 1)

                
    def append_data_set(self, name, val_input):
        val = np.zeros(tuple(n+2 for n in val_input.shape), dtype=val_input.dtype)
        val[tuple([slice(1,-1)]*self.nd)] = val_input
        if 'x' in self.periodic:
            val[0, ...] = val[-2, ...]
            val[-1, ...] = val[1, ...]
        if 'y' in self.periodic:
            val[:,0, ...] = val[:,-2, ...]
            val[:,-1, ...] = val[:,1, ...]
        if 'z' in self.periodic:
            val[:,:,0, ...] = val[:,:,-2, ...]
            val[:,:,-1, ...] = val[:,:,1, ...]

        print()
        for rank in np.arange(self.np):
            self.read_points(rank)
            self.append(rank)
            self.write_data_set_attribute(name, val)
            self.close()
            print( "\rAppended data: {} to file {}".format(name, self.local_filename), end='' )
        
