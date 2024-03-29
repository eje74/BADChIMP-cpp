import numpy as np


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
        self.label = -np.ones(self.geo.shape, dtype=int)
        for rank in np.arange(1,np.amax(geo)+1):
            # Find the bounding box for a given rank
            ind_tmp_tuble = np.where(self.geo == rank)
            ind_tmp = np.array(ind_tmp_tuble, dtype=int) 
            inds_min = np.amin(ind_tmp, axis=1)
            inds_max = np.amax(ind_tmp, axis=1)
            sys_size = (inds_max - inds_min + 1) + 2*self.rim_width
            # The + (1,) is to uphold the broadcasting rules
            ind_tmp = ind_tmp - inds_min.reshape( inds_min.shape + (1,) ) + self.rim_width
            # Generate and linear indexing ind = nx + ny*NX + nz*NX*NY
            prod_lin_index = np.roll(sys_size, 1)
            prod_lin_index[0] = 1
            prod_lin_index = np.cumprod(prod_lin_index).reshape(prod_lin_index.shape + (1,))
            self.label[ind_tmp_tuble] = np.sum(ind_tmp*prod_lin_index, axis=0)
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
            print("Could not open file: {}".format(self.path+self.filename)) 


    def append(self, rank):        
        self.local_filename = self.filename + str(rank) + ".vtklb"
        try:
            self.file = open(self.path+self.local_filename, "a")
        except IOError:
            print("Could not open file: {}".format(self.path+self.filename)) 

        
    def read(self, rank):
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

                
    def setup_processor_labels_reg(self, rank):
        # input: 'rank' as written in the geo file
        # output: 'local_label' labeling of nodes on the processor
        #           Needs to match the global labeling. 
        #         'rank_set' list of neighboring processors
        #           Need to check if they can be reached from the bulk nodes 
        #           , that is, nodes that are not part of the rim.
        # Find the indecies and use them as position
        self.ind = np.where( (self.geo == rank) & self.bulk)        
        # List of local indecies
        self.ind_local = tuple(np.copy(self.ind))        
        # Find solids        
        for v in self.basis:
            ind_tmp = tuple( ind + dind for dind, ind in zip(v, self.ind) )
        # Insert solid wall nodes into self.ind
        self.ind_local = tuple(np.concatenate((ind1, ind2)) for ind1, ind2 in zip(self.ind_local, np.where(self.added)))
        # Generate a list of neighboring processors
        self.rank_set = set()
        for v in self.basis:
            ind_tmp = tuple(ind + dind for dind, ind in zip(v, self.ind))
            self.rank_set = self.rank_set.union(set(self.geo[ind_tmp]))
        #   remove ghost, solid and current rank from the set     
        self.rank_set.discard(-1)
        self.rank_set.discard(0)
        self.rank_set.discard(rank)
        # set the local labels (use -1 as marking regions outside of the current processor)
        self.local_label = -np.ones(self.label.shape, dtype=int) 
        ind_tmp_tuble = np.where((self.geo == rank) & self.bulk)
        ind_tmp = np.array(ind_tmp_tuble, dtype=int) 
        inds_min = np.amin(ind_tmp, axis=1) - self.rim_width
        inds_max = np.amax(ind_tmp, axis=1) + 1 + self.rim_width
        self.sys_size = (inds_max - inds_min) 
        slice_list =  tuple(slice(sb, se) for sb, se in zip(inds_min, inds_max))   
        label_shape = tuple(se-sb for sb, se in zip(inds_min, inds_max))
        self.local_label[slice_list] = np.arange(np.prod(self.sys_size)).reshape(label_shape, order='F')
        self.global_ind = np.copy(inds_min)        

                
    def write_preamble(self, rank):
        self.file.write( "# BADChIMP vtklb Version {}\n".format(self.version) )
        self.file.write( "Geometry file for process {}\n".format(rank-1) )      
        self.file.write( "ASCII\n")


    def read_local_data(self, rank):
        """
        Reads the prosessor spesific data. 
        Assumptions: 
            - self.nd is set to the correct number of dimensions
        Actions: 
            - reads 'self.sys_size', 'self.global_ind', 'self.rim_width' from
              file and set corresponding variables to these values            
        """
        self.read(rank)
        line = self.file.readline()
        words = line.split()
        while words[0] != "LOCAL_DIMENSIONS":
            line = self.file.readline()
            words = line.split()
        # Set local dimensions
        for dim, word in enumerate(words[1:]):
            self.sys_size[dim] = int(word)
        while words[0] != "LOCAL_RIM_WIDTH":
            line = self.file.readline()
            words = line.split()
        self.rim_width = int(words[1])
        while words[0] != "LOCAL_TO_GLOBAL":
            line = self.file.readline()
            words = line.split()
        for dim, word in enumerate(words[1:]):
            self.global_ind[dim] = int(word)
        # Add the rim width to the size
        self.sys_size += 2*self.rim_width
        self.close()

                
    def write_dataset_reg(self):
        self.file.write( "DATASET STRUCTURED_LB_GRID\n" )
        self.file.write( "NUM_DIMENSIONS {}\n".format(self.geo.ndim) )
        line = "GLOBAL_DIMENSIONS"
        for nd in range(self.geo.ndim):
            line += " {}".format(self.geo.shape[nd])
        self.file.write(line + "\n");
        line = "LOCAL_DIMENSIONS"
        for nd in range(self.geo.ndim):
            line += " {}".format(self.sys_size[nd] - 2*self.rim_width)
        self.file.write(line + "\n");
        line = "LOCAL_RIM_WIDTH"
        line += " {}".format(self.rim_width)
        self.file.write(line + "\n");
        line = "LOCAL_TO_GLOBAL"
        for nd in range(self.geo.ndim):
            line += " {}".format(self.global_ind[nd])
        self.file.write(line + "\n");
        
        
    def write_lattice(self):
        self.file.write( "LATTICE {} int\n".format(len(self.basis)) )
        np.savetxt(self.file, self.basis, fmt="%d", delimiter=' ', newline='\n')    
                                                                                                                                                                     
                        
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


    def write_mpi_reg(self, rank):
        # Periodic nodes
        neig_ind = np.where( (self.local_label > -1) & (~self.bulk) & (self.geo == rank) )
        self.file.write( "PERIODIC_NODES {} {}\n".format(len(neig_ind[0]), rank-1 ))
        np.savetxt(self.file, np.array([self.local_label[neig_ind], self.label[neig_ind]]).transpose(), fmt="%d", delimiter=' ', newline='\n')            
        # MPI header
        self.file.write( "PARALLEL_COMPUTING {}\n".format(rank-1) )
        # Go through the neighboring ranks
        for neig_rank in self.rank_set:
            neig_ind = np.where( (self.local_label > -1) & (self.geo == neig_rank) )
            # Rank header
            self.file.write( "PROCESSOR {} {}\n".format(len(neig_ind[0]), neig_rank-1) )
            # Write list of overlapping nodes [local label, neigh process label]
            np.savetxt(self.file, np.array([self.local_label[neig_ind], self.label[neig_ind]]).transpose(), fmt="%d", delimiter=' ', newline='\n')            

        
    def write_point_data(self):
        # Header for point data
        self.file.write( "POINT_DATA {}\n".format( len(self.ind_local[0])) )

                                        
    def write_data_set_attribute_reg(self, name, val):
        slice_list =  tuple(slice(sb, sb + ds) for sb, ds in zip(self.global_ind, self.sys_size))
        if np.issubdtype(val.dtype, np.integer):
            self.file.write( "SCALARS {} int\n".format(name) )
            np.savetxt(self.file, val[slice_list].flatten(order='F').transpose(), fmt="%d", delimiter=' ', newline='\n')
        else:
            self.file.write( "SCALARS {} float\n".format(name) )
            np.savetxt(self.file, val[slice_list].flatten(order='F').transpose(), delimiter=' ', newline='\n')
                 
             
    def write_proc(self, rank):
        self.open(rank)
        self.setup_processor_labels_reg(rank)
        self.write_preamble(rank)
        self.write_dataset_reg() 
        self.write_lattice()
        self.write_mpi_reg(rank)
        self.write_point_data()
        # Write the node types to file
        self.write_data_set_attribute_reg('nodetype', self.geo)
        self.close()
        print( "Wrote file: {}".format(self.local_filename) )

                
    def write(self):
        print( "Writing files to folder {}".format(self.path) )
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

        for rank in np.arange(self.np):
            self.read_local_data(rank)
            self.append(rank)
            self.write_data_set_attribute_reg(name, val)
            self.close()
            print( "Appended data: {} to file {}".format(name, self.local_filename) )
        
