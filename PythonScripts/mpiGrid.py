import numpy as np
import matplotlib.pyplot as plt

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


def writeNDarray(A, f):
    if A.ndim > 1:
        for x in np.arange(geo.shape[0]):
            writeNDarray(A[x], f)
    else :
        for x in np.arange(A.shape[0]):
            f.write(str(A[x]) + " ")
        f.write("\n")

def writeFile(filename, geo, fieldtype):
    f = open(filename, "w")
    ## Write dimensions x y
    f.write("dim")
    for dim in np.arange(geo.ndim, 0, -1):
        f.write(" " + str(geo.shape[dim-1]))
    f.write("\n")
    f.write("<" + fieldtype + ">\n")

    writeNDarray(geo, f)
    f.write("<end>\n")
    f.close()


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
