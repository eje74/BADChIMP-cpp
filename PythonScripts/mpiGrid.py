#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

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
    #geo_shape = geo.shape # System size
    # Flatten the array
    #geo = geo.reshape([geo.size, ])
    node_labels = np.zeros(geo.shape, dtype=np.int)
    label_counter = np.zeros(num_proc, dtype=np.int)
    for rank in np.arange(num_proc):
        ind = np.where(geo == rank + 1)
        label_counter[rank] = ind[0].size
        node_labels[ind] = np.array(np.arange(label_counter[rank]))
    return (node_labels, label_counter)


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
    """
    # shape without rim (NB important the shape is (ndim, 1) and not (ndim,))
    # when using np.mod
    geo_shape = np.array(geo_rim.shape).reshape((geo_rim.ndim,1)) - 2*rim_width
    # from tuple of array to array of array
    ind = np.array(ind_of_periodic_nodes)
    ind_mod =  np.mod(ind-rim_width, geo_shape) + rim_width
    geo_rim[tuple(ind)] = geo_rim[tuple(ind_mod)]
    return geo_rim


def setNodeLabelsLocal(geo, node_labels, myRank, c):
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
    return fluid_markerA

def writeFile(filename, geo, fieldtype, origo_index, rim_width):
    f = open(filename, "w")
    ## Write dimensions x y
    f.write("dim")
    for dim in np.arange(geo.ndim, 0, -1):
        f.write(" " + str(geo.shape[dim-1]))
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



def setNodeType(rank, c, my_rank):
    # Note: Only 2D
    # NAN (SOLID MPI)      = 0
    # SOLID BOUNDARY       = 1
    # FLUID                = 2
    # MPI BOUNDARY         = 3

    node_type = np.zeros(np.array(rank.shape) - 2, dtype=np.int)
    for y in np.arange(rank.shape[0]-2):
        for x in np.arange(rank.shape[1]-2):
            x_node = x + 1
            y_node = y + 1
            rank_node = rank[y_node, x_node]
            neig_node = np.array([rank[y_node + cq[1], x_node + cq[0]] for cq in c[:-1]])
            if rank_node == 0: # either NAN (solid) or solid boundary
                if np.any(neig_node == my_rank):
                    node_type[y, x] = 1
            elif rank_node == my_rank: # fluid
                node_type[y, x] = 2
            else: # mpi_rank != my_rank, either NAN (mpi) or MPI BOUNDARY
                if np.any(neig_node == my_rank):
                    node_type[y, x] = 3
    return node_type


def setNodeType3D(rank, c, my_rank):
    # Note: Only 3D
    # NAN (SOLID MPI)      = 0
    # SOLID BOUNDARY       = 1
    # FLUID                = 2
    # MPI BOUNDARY         = 3

    node_type = np.zeros(np.array(rank.shape) - 2, dtype=np.int)
    for z in  np.arange(rank.shape[0]-2):
        for y in np.arange(rank.shape[1]-2):
            for x in np.arange(rank.shape[2]-2):
                x_node = x + 1
                y_node = y + 1
                z_node = z + 1
                rank_node = rank[z_node, y_node, x_node]
                neig_node = np.array([rank[z_node +cq[2], y_node + cq[1], x_node + cq[0]] for cq in c[:-1]])
                if rank_node == 0: # either NAN (solid) or solid boundary
                    if np.any(neig_node == my_rank):
                        node_type[z, y, x] = 1
                elif rank_node == my_rank: # fluid
                    node_type[z, y, x] = 2
                else: # mpi_rank != my_rank, either NAN (mpi) or MPI BOUNDARY
                    if np.any(neig_node == my_rank):
                        node_type[z, y, x] = 3
    return node_type




def setBoundaryLabels(node_types, labels, my_rank):
    # Note: Only 2D
    # Need to set SOLID BOUNDARY LABELS
    # Need to set MPI BOUNDARY LABELS (for each rank != my_rank)
    # Possible challanges: ?
    label_counter = np.max(labels[node_types == 2])
    labels_including_bnd = np.zeros(np.array(labels.shape) - 2, dtype=np.int)
    for y in np.arange(labels.shape[0]-2):
        for x in np.arange(labels.shape[1]-2):
            x_node = x + 1
            y_node = y + 1
            node_type = node_types[y_node, x_node]
            if node_type == 2:  # FLUID
                labels_including_bnd[y, x] = labels[y_node, x_node]
            elif node_type == 1: #SOLID BND
                label_counter += 1
                labels_including_bnd[y, x] = label_counter
            elif node_type == 3: # MPI BND
                label_counter += 1
                labels_including_bnd[y, x] = label_counter
    return labels_including_bnd


def setBoundaryLabels3D(node_types, labels, my_rank):
    # Note: Only 3D
    # Need to set SOLID BOUNDARY LABELS
    # Need to set MPI BOUNDARY LABELS (for each rank != my_rank)
    # Possible challanges: ?
    label_counter = np.max(labels[node_types == 2])
    labels_including_bnd = np.zeros(np.array(labels.shape) - 2, dtype=np.int)
    for z in np.arange(labels.shape[0]-2):
        for y in np.arange(labels.shape[1]-2):
            for x in np.arange(labels.shape[2]-2):
                x_node = x + 1
                y_node = y + 1
                z_node = z + 1
                node_type = node_types[z_node, y_node, x_node]
                if node_type == 2:  # FLUID
                    labels_including_bnd[z, y, x] = labels[z_node, y_node, x_node]
                elif node_type == 1: #SOLID BND
                    label_counter += 1
                    labels_including_bnd[z, y, x] = label_counter
                elif node_type == 3: # MPI BND
                    label_counter += 1
                    labels_including_bnd[z, y, x] = label_counter
    return labels_including_bnd


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

#write_dir = "/home/olau/Programs/Git/BADChIMP-cpp/PythonScripts/"
#write_dir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/PythonScripts/"
write_dir = "/home/ejette/Programs/GITHUB/badchimpp/PythonScripts/"


# SETUP GEOMETRY with rank (0: SOLID, 1:RANK0, 2:RANK1, ...)
geo_input = readGeoFile(write_dir + "test.dat") # assumes this shape of geo_input [(nZ, )nY, nX]
# -- setup domain-decomposition (this could be written in the geo-file)
geo_input[:, 67:134, :] = 2*geo_input[:, 67:134, :]
geo_input[:, 134:, :] = 3*geo_input[:, 134:, :]
# -- derive the number of processors and the rim width
num_proc = getNumProc(geo_input)
rim_width = getRimWidth(c)

# Add rim
geo = addRim(geo_input, rim_width)


# Create node labels
node_labels, num_labels = setNodeLabels(geo_input, num_proc)
node_labels = addRim(node_labels, rim_width)

## Here we need to set aditional values for the rim
# set as periodic -1 (default value)
# set as solid = -2
# set as fluid = assign -#rank - 3


geo[:, 0, :] = -2
geo[:, -1, :] = -2


ind_periodic = np.where(geo == -1)
ind_solid = np.where(geo == -2)
ind_fluid = np.where(geo < - 2)

geo[ind_solid] = 0
geo = addPeriodicBoundary(ind_periodic, geo, rim_width)
node_labels[ind_solid] = 0
node_labels = addPeriodicBoundary(ind_periodic, node_labels, rim_width)

#geo = addPeriodicBoundary(geo, rim_width)
#node_labels = addRim(node_labels, rim_width)
#node_labels = addPeriodicBoundary(node_labels, rim_width)
#node_labels = addPeriodicRim(node_labels)



node_labels_local = setNodeLabelsLocal(geo, node_labels, 1, c)

plt.figure(1)
plt.pcolormesh(geo[4, :,:])
plt.figure(2)
plt.pcolormesh(node_labels[4, :, :])
plt.figure(3)
plt.pcolormesh(node_labels_local[4, :, :])

plt.show()


# Write files GEO and NODE LABELS
# origo_index = np.array([0, 0])
origo_index = np.array([0, 0, 0])
#rim_width = 1

writeFile(write_dir + "rank.mpi", geo, "rank int", origo_index, rim_width)
writeFile(write_dir + "node_labels.mpi", node_labels, "label int", origo_index, rim_width)



for my_rank in np.arange(1, num_proc + 1):
    node_types = setNodeType3D(geo, c, my_rank)
    node_types = addPeriodicRim(node_types)

    node_labels_extended = setBoundaryLabels3D(node_types, node_labels, geo)
    node_labels_extended = addPeriodicRim(node_labels_extended)


    rank_label_file_name = "rank_" + str(my_rank-1) + "_labels.mpi";
    writeFile(write_dir + rank_label_file_name, node_labels_extended, "label int", origo_index, rim_width)
