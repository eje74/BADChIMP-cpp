import numpy as np
import matplotlib.pyplot as plt

#c = np.array([[1, 0],   [1, 1],  [0, 1],  [-1, 1],
#              [-1, 0], [-1, -1], [0, -1], [1, -1], [0, 0]])

c = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0],
              [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1],
	      [0, 1, -1], [-1, 0, 0], [0, -1, 0], [0, 0, -1],
              [-1, -1, 0], [-1, 1, 0], [-1, 0, -1], [-1, 0, 1],
              [0, -1, -1], [0, -1, 1], [0, 0, 0]])

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
    for bnd_ind in np.arange(ind.shape[1]): # ind.shape[1] : number of boundary points
        pnt = tuple(ind[:, bnd_ind])
        pnt_mod = [np.mod(pnt[d]-1, A_shape[d]) for d in np.arange(A.ndim)]
        pnt_mod = tuple(pnt_mod)
        A_periodic[pnt] = A[pnt_mod]
    return A_periodic



def writeFileOld(filename, geo, fieldtype, origo_index, rim_width):
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
            for x in np.arange(geo.shape[0]):
                writeLoop(A[x])
        else :
            for x in np.arange(A.shape[0]):
                f.write(str(A[x]) + " ")
            f.write("\n")

    writeLoop(geo)

    f.write("<end>\n")
    f.close()


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


def readGeoFile3D(file_name):
    f = open(file_name)
    line = f.readline().rstrip('\n')
    dim = []
    for i in line.split(" "):
        if i.isdigit():
            dim.append(int(i))
    line = f.readline().rstrip('\n')
    line = f.readline().rstrip('\n')
    ret = np.zeros((dim[0]*dim[1]*dim[2], ), dtype=np.int)
    cnt = 0
    for i in np.arange(dim[1]*dim[2]):
        line = f.readline().rstrip('\n')
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

write_dir = "/home/ejette/Programs/GITHUB/badchimpp/PythonScripts/"
#write_dir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/PythonScripts/"

geo_input = readGeoFile3D(write_dir + "test.dat")


# NN = [nY, nX]
# NN = [7, 12]
NN = list(geo_input.shape)

# Setup geomtry with rank (0: SOLID, 1:RANK0, 2:RANK1, ...)

geo_input[:, 266:534, :] = 2*geo_input[:, 266:534, :]
geo_input[:, 534:, :] = 3*geo_input[:, 534:, :]


# ANALYSE GEOMETRY
num_proc = 3



# Create node labels
node_labels, num_labels = setNodeLabels(geo_input, num_proc)
# Add periodic rim
geo = addPeriodicRim(geo_input)
node_labels = addPeriodicRim(node_labels)

# Write files GEO and NODE LABELS
# origo_index = np.array([0, 0])
origo_index = np.array([0, 0, 0])
rim_width = 1

writeFile(write_dir + "rank.mpi", geo, "rank int", origo_index, rim_width)
writeFile(write_dir + "node_labels.mpi", node_labels, "label int", origo_index, rim_width)



for my_rank in np.arange(1, num_proc + 1):
    node_types = setNodeType3D(geo, c, my_rank)
    node_types = addPeriodicRim(node_types)

    node_labels_extended = setBoundaryLabels3D(node_types, node_labels, geo)
    node_labels_extended = addPeriodicRim(node_labels_extended)


    rank_label_file_name = "rank_" + str(my_rank-1) + "_labels.mpi";
    writeFile(write_dir + rank_label_file_name, node_labels_extended, "label int", origo_index, rim_width)
