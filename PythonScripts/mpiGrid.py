import numpy as np
import matplotlib.pyplot as plt

c = np.array([[1, 0],   [1, 1],  [0, 1],  [-1, 1],
              [-1, 0], [-1, -1], [0, -1], [1, -1], [0, 0]])


def setNodeLabels(geo, num_proc):
    geo_shape = geo.shape # System size
    # Flatten the array
    geo = geo.reshape([geo.size, ])
    node_labels = np.zeros(geo.shape, dtype=np.int)
    label_counter = np.zeros(num_proc, dtype=np.int)
    for current_node in np.arange(geo.size):
        rank = geo[current_node]
        if rank > 0:
            rank -= 1
            label_counter[rank] += 1
            node_labels[current_node] = label_counter[rank]
    return (node_labels.reshape(geo_shape), label_counter)


def addPeriodicRim(A):
    A_shape = np.array(A.shape)
    A_periodic = np.zeros(A_shape + 2, dtype=np.int)
    A_periodic[(slice(1,-1),)*A.ndim] = 1
    ind = np.array(np.where(A_periodic == 0))
    A_periodic[(slice(1,-1),)*A.ndim] = A
    for bnd_ind in np.arange(ind.shape[1]):
        pnt = tuple(ind[:, bnd_ind])
        pnt_mod = [np.mod(pnt[d]-1, A_shape[d]) for d in np.arange(A.ndim)]
        pnt_mod = tuple(pnt_mod)
        A_periodic[pnt] = A[pnt_mod]
    return A_periodic



def writeFile(filename, geo, fieldtype):

    f = open(filename, "w")
    ## Write dimensions x y
    f.write("dim")
    for dim in np.arange(geo.ndim, 0, -1):
        f.write(" " + str(geo.shape[dim-1]))
    f.write("\n")
    f.write("<" + fieldtype + ">\n")

    def writeLoop(A):
        if A.ndim > 1:
            for x in np.arange(geo.shape[0]):
                writeLoop(A[x])
        else :
            for x in np.arange(A.shape[0]):
                f.write(str(A[x]) + " ")
            f.write("\n")

    writeLoop(geo)

    f.write("<end>\n")
    f.close()


def setNodeType(rank, my_rank):
    # NAN (SOLID MPI)      = 0
    # SOLID BOUNDARY       = 1
    # FLUID                = 2
    # MPI BOUNDARY         = 3

    node_type = np.zeros(rank.shape, dtype=np.int)
    for y in np.arange(rank.shape[0]-2):
        for x in np.arange(rank.shape[1]-2):
            x_node = x + 1
            y_node = y + 1
            rank_node = rank[y_node, x_node]
            neig_node = np.array([rank[y_node + cq[1], x_node + cq[0]] for cq in c[:-1]])
            if rank_node == 0: # either NAN (solid) or solid boundary
                if np.any(neig_node == my_rank):
                    node_type[y_node, x_node] = 1
            elif rank_node == my_rank: # fluid
                node_type[y_node, x_node] = 2
            else: # mpi_rank != my_rank, either NAN (mpi) or MPI BOUNDARY
                if np.any(neig_node == my_rank):
                    node_type[y_node, x_node] = 3

    return node_type

def setBoundaryLabels(types, labels):
    pass

write_dir = "/home/ejette/Programs/GITHUB/badchimpp/PythonScripts/"

NN = [7, 12]
geo_input = np.zeros(NN, dtype=np.int)

# setup basic geomtry with rank (0: SOLID, 1:RANK0, 2:RANK1, ...)
# Row 0, 1
geo_input[:2, :2]  = 2
geo_input[:2, 2:6] = 1
geo_input[:2, 6:]  = 2
# Row 2, 3
geo_input[2:4, :2]  = 1
geo_input[2:4, 4:8] = 1
geo_input[2:4, 8:]  = 2
# Row 4
geo_input[4,:9] = 1
geo_input[4, 9] = 2
# Row 5
geo_input[5, :7]  = 1
geo_input[5,7:10] = 2
# Row 6
geo_input[6,:5] = 1
geo_input[6,6:] = 2


# ANALYSE GEOMETRY
# Create node labels
node_labels, num_labels = setNodeLabels(geo_input, 2)
# Add periodic rim
geo = addPeriodicRim(geo_input)
node_labels = addPeriodicRim(node_labels)
# Write files
writeFile(write_dir + "rank.mpi", geo, "rank int")
writeFile(write_dir + "node_labels.mpi", node_labels, "label int")

node_types = setNodeType(geo, 2)
print(node_labels)
print(geo)
print(node_types)
