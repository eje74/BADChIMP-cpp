#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from sys import argv

#c = np.array([[1, 0],   [1, 1],  [0, 1],  [-1, 1],
#              [-1, 0], [-1, -1], [0, -1], [1, -1], [0, 0]])

c = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0],
              [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1],
	      [0, 1, -1], [-1, 0, 0], [0, -1, 0], [0, 0, -1],
              [-1, -1, 0], [-1, 1, 0], [-1, 0, -1], [-1, 0, 1],
              [0, -1, -1], [0, -1, 1], [0, 0, 0]])


# Simple utility functions
def getNumProc(geo):
    """Finds the number of processes used in geo
    Assumes that 0 is used to give a non standard fluid type
    and the ranks are number form 1.
    """
    # max = number of ranks, since rank identifiers are index from 1 and not 0.
    return np.amax(geo)


def getRimWidth(c):
    """Sets the width of the system based on size of the LB lattice
    c : are the velocity basis vectors
    """
    return np.amax(np.abs(c))


def setNodeLabels(geo, num_proc):
    # Sets the node labels for the fluid nodes
    node_labels = np.zeros(geo.shape, dtype=np.int)
    for rank in np.arange(num_proc):
        ind = np.where(geo == rank + 1)
        # Remember that the 0 is the default node number
        node_labels[ind] = 1 + np.array(np.arange(ind[0].size))
    return node_labels


def addRim(geo, rim_width):
    """Add a rim to the geometry.
    -1: markes the rim
    """
    geo_rim = -np.ones(np.array(geo.shape) + 2*rim_width, dtype=np.int)
    geo_rim[(slice(rim_width,-rim_width),)*geo.ndim] = geo
    return geo_rim


def addPeriodicBoundary(ind_of_periodic_nodes, geo_rim, rim_width):
    """Add periodic boundary for PM values in geo
    We assume that geo has been extended with rim.

    ind_of_periodic_nodes : give by np.array(np.where(geo_rim == PM))
        where 'PM' is the periodic node marker.
    geo_rim : map included size for the rim.
    """
    # shape without rim (NB important the shape is (ndim, 1) and not (ndim,))
    # when using np.mod
    geo_shape = np.array(geo_rim.shape).reshape((geo_rim.ndim,1)) - 2*rim_width
    # from tuple of array to array of array
    ind = np.array(ind_of_periodic_nodes)
    ind_mod =  np.mod(ind-rim_width, geo_shape) + rim_width
    geo_rim[tuple(ind)] = geo_rim[tuple(ind_mod)]
    return geo_rim


def setNodeLabelsLocal(geo, node_labels, ind_periodic, myRank, c, rim_width):
    """Sets the local node labels. Must assign local
    labels to the solid boundary, and to mpi boundary
    """
    fluid_markerA = geo == myRank
    fluid_markerB = np.copy(fluid_markerA)

    for q in np.arange(c.shape[0]):
        slicerA = ()
        slicerB = ()
        for d in np.arange(geo.ndim):
            cval = c[q, d]
            if cval > 0:
                sliceA = slice(cval, None, None)
                sliceB = slice(0, -cval, None)
            elif cval < 0:
                sliceA = slice(0, cval, None)
                sliceB = slice(-cval, None, None)
            else:
                sliceA = slice(0, None, None)
                sliceB = slice(0, None, None)
            slicerA = (sliceA,) + slicerA
            slicerB = (sliceB,) + slicerB
        fluid_markerA[slicerA] = np.logical_or(fluid_markerA[slicerA], fluid_markerB[slicerB])
    # Nodes that need new labels
    fluid_markerA = np.logical_and(fluid_markerA, np.logical_not(fluid_markerB))
    # Periodic nodes are not given a new label
    fluid_markerA[ind_periodic] = False
    # Set values bulk labels
    node_labels_local = np.zeros(fluid_markerA.shape, dtype=node_labels.dtype)
    node_labels_local[fluid_markerB] = node_labels[fluid_markerB]
    # Find largest label
    max_label = np.amax(node_labels_local)
    # Find the boundary nodes that will be given a local label
    ind_loc_lab = np.where(fluid_markerA)
    # Set labels local to the given process
    node_labels_local[ind_loc_lab] = max_label + 1 + np.arange(ind_loc_lab[0].size)
    # Add the periodic boundaries
    addPeriodicBoundary(ind_periodic, node_labels_local, rim_width)

    return node_labels_local

def writeFile(filename, geo, fieldtype, global_dim, origo_index, rim_width):
    f = open(filename, "w")
    ## Write dimensions x y
    f.write("dim")
    for dim in np.arange(geo.ndim, 0, -1):
        f.write(" " + str(geo.shape[dim-1]))
    f.write("\n")
    ## Write the dimension for the whole sysmte
    f.write("dim_global")
    for dim in np.arange(geo.ndim, 0, -1):
        f.write(" " + str(global_dim[dim-1]))
    f.write("\n")
    f.write("origo")
    for dim in np.arange(geo.ndim, 0, -1):
        f.write(" " + str(origo_index[dim-1]))
    f.write("\n")
    f.write("rim " + str(rim_width))
    f.write("\n")
    f.write("<" + fieldtype + ">\n")

    def writeLoop(A):
        if A.ndim > 1:
            for x in np.arange(A.shape[0]):
                writeLoop(A[x])
        else :
            for x in np.arange(A.shape[0]):
                f.write(str(A[x]) + " ")
            f.write("\n")

    writeLoop(geo)

    f.write("<end>\n")
    f.close()


def readGeoFile(file_name):
    f = open(file_name)
    line = f.readline().rstrip('\n')
    dim = []
    for i in line.split(" "):
        if i.isdigit():
            dim.append(int(i))
    line = f.readline().rstrip('\n')
    line = f.readline().rstrip('\n')
    ret = np.zeros((np.prod(dim), ), dtype=np.int)
    cnt = 0
    while cnt < ret.size:
    #for line in f:
        line = f.readline().rstrip('\n')
        #line = line.rstrip('\n')
        for c in line:
            ret[cnt] = int(c)
            cnt += 1
    f.close()
    dim = dim[-1::-1]
    ret = ret.reshape(dim)

    ret[ret == 1] = 2
    ret[ret == 0] = 1
    ret[ret == 2] = 0

    return ret

write_dir = "/home/olau/Programs/Git/BADChIMP-cpp/PythonScripts/"
#write_dir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/PythonScripts/"  # Home
#write_dir = "/home/ejette/Programs/GITHUB/badchimpp/PythonScripts/"  # Work
#write_dir = "PythonScripts/"  # Work

# geo.shape = (NZ, NY, NX)
# SETUP GEOMETRY with rank (0: SOLID, 1:RANK0, 2:RANK1, ...)
#geo_input = readGeoFile(write_dir + "test.dat") # assumes this shape of geo_input [(nZ, )nY, nX]
file_name = "test.dat"
#file_name = "frac.dat"

procs = np.array((1,3,1))

if len(argv)>1:
    if '/' in argv[1]:
        file_name = argv[1].split('/')[-1]
        write_dir = '/'.join(argv[1].split('/')[:-1]) + '/'
    else:
        file_name = argv[1]
print('Reading ' + write_dir + file_name + ' ...')
geo_input = readGeoFile(write_dir + file_name) # assumes this shape of geo_input [(nZ, )nY, nX]

# -- option to provide domain decomposition via commandline
if len(argv)>2:
    if len(argv)<4:
        print(' Missing arguments!')
        print(' Usage: ')
        print('    mpiGrid.py  filename ')
        print('    mpiGrid.py  filename proc_x proc_y proc_z')
    procs = np.asarray([int(i) for i in argv[2:5]])[::-1] # -- reverse input so that first index is x

# -- setup domain-decomposition (this could be written in the geo-file)
dim = np.asarray(geo_input.shape)
print(dim)
step = (dim-1)//procs + 1 # does ceil division using floor division (//)
nproc = 2
for d in range(len(dim)):
    ax = [np.arange(0,dim[i]) for i in range(len(dim))]
    for n in range(1, procs[d]):
        if n<procs[d]-1:
            ax[d] = np.arange(step[d])+n*step[d]
            print(str(nproc)+': ax['+str(d)+'] = np.arange('+str(step[d])+')+'+str(n*step[d]))
        else:
            # -- open end-range for last process
            ax[d] = np.arange(n*step[d],dim[d])
            print(str(nproc)+': ax['+str(d)+'] = np.arange('+str(n*step[d])+','+str(dim[d])+')')
        geo_input[ax[0][:,np.newaxis,np.newaxis], ax[1][np.newaxis,:,np.newaxis], ax[2][np.newaxis,np.newaxis,:]] *= nproc
        nproc += 1
       
plt.figure(10)
pltmat = geo_input
print 'whoop'
print geo_input.shape
plt.pcolormesh(pltmat[0, :,:])#plt.pcolormesh(pltmat[4, :,:])
plt.colorbar()       
                  
#geo_input[:, 67:134, :] = 2*geo_input[:, 67:134, :]
#geo_input[:, 134:, :] = 3*geo_input[:, 134:, :]

#ind_zero = np.where(geo_input == 0)
#geo_input = geo_input - 1
#geo_input[:,:,0:50] = geo_input[:,:,0:50] + 3
#geo_input = geo_input + 1
#geo_input[ind_zero] = 0

# -- derive the number of processors and the rim width
rim_width = getRimWidth(c)
# -- add rim
geo = addRim(geo_input, rim_width)

print geo.shape

# SETUP RIM VALUES
## Here we need to set aditional values for the rim
# set as periodic -1 (default value)
# set as solid = -2
# set as fluid = -<#rank> - 3

# -- BEGIN user input

# -- close system 
#geo[:, 0, :] = -2
#geo[:, -1, :] = -2
#geo[:, :, 0] = -2
#geo[:, :, -1] = -2


# -- END user input


ind_periodic = np.where(geo == -1)
ind_solid = np.where(geo == -2)
ind_fluid = np.where(geo < - 2)

# Set the rim solid nodes
geo[ind_solid] = 0
# Set the rim fluid nodes
geo[ind_fluid] = -geo[ind_fluid] - 2


# Set number of processes
num_proc = getNumProc(geo_input)


# Make node labels
node_labels = setNodeLabels(geo, num_proc)

geo = addPeriodicBoundary(ind_periodic, geo, rim_width)
node_labels[ind_solid] = 0
node_labels = addPeriodicBoundary(ind_periodic, node_labels, rim_width)




# Write files GEO and NODE LABELS

plt.figure(1)
pltmat = geo
plt.pcolormesh(pltmat[1, :,:])#plt.pcolormesh(pltmat[4, :,:])
plt.colorbar()
plt.figure(2)
pltmat = node_labels
plt.pcolormesh(pltmat[1, :, :])#plt.pcolormesh(pltmat[4, :,:])
plt.colorbar()




for my_rank in np.arange(1, num_proc + 1):

    print("AT: " + str(my_rank-1))

    node_labels_local = setNodeLabelsLocal(geo, node_labels, ind_periodic, my_rank, c, rim_width)

    # The length of this tuple gives the number of dimensions
    #  remember [0]:list of z, [1]: list of y and [2] list of x
    ind_labels_rank = np.where(node_labels_local > 0)
    slicer = tuple(slice(imin, imax+1) for imin, imax in zip(np.min(ind_labels_rank, axis=1), np.max(ind_labels_rank, axis=1)))

    origo_index = np.min(ind_labels_rank, axis=1)

    file_name = "rank_" + str(my_rank-1) + "_rank.mpi";
    writeFile(write_dir + file_name, geo[slicer], "rank int", geo.shape, origo_index, rim_width)
    file_name = "rank_" + str(my_rank-1) + "_global_labels.mpi";
    writeFile(write_dir + file_name, node_labels[slicer], "label int", geo.shape, origo_index, rim_width)
    file_name = "rank_" + str(my_rank-1) + "_local_labels.mpi";
    writeFile(write_dir + file_name, node_labels_local[slicer], "local label int", geo.shape, origo_index, rim_width)

    plt.figure(3)
    pltmat = node_labels_local[slicer]
    plt.pcolormesh(pltmat[0, :, :])#plt.pcolormesh(pltmat[4, :,:])


plt.show()
