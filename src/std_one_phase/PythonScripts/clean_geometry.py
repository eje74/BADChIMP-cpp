import numpy as np
from numba import njit
from numba import prange
from numba.typed import List
from itertools import product

@njit
def bwlabel3D_set(cc, geo, marker):
    """
    """
    I, J, K = geo.shape # System size
    # ---------------------------------------------------------------------- Find set of known neighborhood
    tmp = []
    for c in cc:
        if c[0] < 0:
            tmp.append(c)
        elif (c[0] == 0) & (c[1] < 0):
            tmp.append(c)
        elif (c[0] == 0) & (c[1] == 0) & (c[2] < 0):
            tmp.append(c)
    known_neighbors = np.zeros((len(tmp), 3), dtype=np.int32)
    for n, c in enumerate(tmp):
        known_neighbors[n, :] = c
    # ---------------------------------------------------------------------- Find connected clusters
    label = 1
    domains = np.zeros(geo.shape, dtype=np.int32)
    label_pairs = List()
    for i in prange(I):
        for j in prange(J):
            for k in prange(K):
                if geo[i, j, k] in marker:    
                    neig_labels = List() 
                    for c in known_neighbors:
                        ni = i + c[0]
                        nj = j + c[1]
                        nk = k + c[2]
                        if (ni >= 0) and (ni < I) and (nj >= 0) and (nj < J) and (nk >= 0) and (nk < K) and (geo[ni, nj, nk] in marker):
                            neig_labels.append(domains[ni, nj, nk])
                    if len(neig_labels) == 0:
                        domains[i, j, k] = label
                        label += 1
                    else:
                        min_label = min(neig_labels)
                        domains[i, j, k] = min_label
                        for a in neig_labels:
                            if a != min_label:
                                label_pairs.append(List([a, min_label]))
    # ---------------------------------------------------------------------- Relabel clusters
    labels = np.zeros(label, dtype=np.int64)
    for n in prange(label):
        labels[n] = n

    for a, b in label_pairs:
        amin = a
        while labels[amin] != amin:
            amin = labels[amin]
        bmin = b
        while labels[bmin] != bmin:
            bmin = labels[bmin]
        if bmin < amin:
            labels[amin] = bmin
        else:
            labels[bmin] = amin

    highest_value_set = -1
    new_label = -1
    new_labels = np.zeros(label, dtype=np.int64)
    for n, val in enumerate(labels):
        valmin = val
        while labels[valmin] != valmin:
            valmin = labels[valmin]

        if valmin > highest_value_set:
            highest_value_set = valmin
            new_label += 1
            new_labels[n] = new_label
        else:
            new_labels[n] = new_labels[valmin]
    # ---------------------------------------------------------------------- Relabel clusters
    for i in prange(I):
        for j in prange(J):
            for k in prange(K):
                domains[i, j, k] = new_labels[domains[i, j, k]]   
                         
    return domains, new_label

@njit
def find_boundary_nodes_set_3D_njit(cc, geo, fluid_marker_set, solid_marker_set):
    I, J, K = geo.shape # System size
    boundary_nodes = np.zeros(geo.shape, dtype=np.int32)
    # ---------------------------------------------------------------------- Find boundary nodes
    for i in prange(I):
        for j in prange(J):
            for k in prange(K):
                if geo[i, j, k] in fluid_marker_set:    
                    for c in cc:
                        ni = i + c[0]
                        nj = j + c[1]
                        nk = k + c[2]
                        # if (ni >= 0) and (ni < I) and (nj >= 0) and (nj < J) and (nk >= 0) and (nk < K) and (geo[ni, nj, nk] in solid_marker_set):
                        #     boundary_nodes[i, j, k] = 1
                        if (ni < 0) or (ni >= I) or (nj < 0) or (nj >= J) or (nk < 0) or (nk >= K) or (geo[ni, nj, nk] in solid_marker_set):
                            boundary_nodes[i, j, k] = 1
    return boundary_nodes

def find_boundary_nodes_set_3d(cc, geo, fluid_marker_set, solid_marker_set):
    nodes = np.transpose(np.nonzero(find_boundary_nodes_set_3D_njit(cc, geo, fluid_marker_set, solid_marker_set)))
    return set([tuple(x) for x in nodes])

def clean_geometry(cc, cc_boundary, geo_org, cluster_cutoff_size=10):
    """
    Cleans the input geometry:
    - Removes the non-percolating clusters of pore space
    - Removes isolated fluid phases that contains less than 
        'cluster_cutoff_size' number of bulk nodes, where
        bulk nodes means nodes that are near a boundary. 
    - Fluid that is removed is set to solid, marked with a '0'.

    Parameters
    ----------
    cc : numpy.ndarray[int]
        List of lattice directions with shape (Q, D). Defines 
        the connectivity of the clusters

    cc_boundary : numpy.ndarray[int]
        List of lattice directions with shape (Q, D), that identifies 
        the boundary nodes.
                

    geo_org : numpy.ndarray[int]
        matrix that that contains the markers:
        SOLID : 0
        FLUID PHASES : 1 or 2

    cluster_cutoff_size : int, default = 10
        minimum size of cluster to keep

    Return
    ------
    list : numpy.ndarray[int]
        Cleaned geometry
    """       
    geo = np.copy(geo_org)
    domains, max_label = bwlabel3D_set(cc, geo, (1, 2))
    # ========================================================================== Remove non-percolating pore space
    percolating_domains = set(domains[:, :, 0].flatten()) & set(domains[:, :, -1].flatten())
    non_percolating_domains = set(range(1,max_label+1)) - percolating_domains
    for n in non_percolating_domains:
        geo[domains==n] = 0
    # ========================================================================== Remove small isolated fluid clusters
    # -------------------------------------------------------------------------- Find boundary nodes
    fluid_solid_boundary = set(find_boundary_nodes_set_3d(cc_boundary, geo, (1, 2), (0,)))
    fluid1_fluid2_boundary = set(find_boundary_nodes_set_3d(cc_boundary, geo, (1,), (2,)))
    fluid2_fluid1_boundary = set(find_boundary_nodes_set_3d(cc_boundary, geo, (2,), (1,)))
    fluid_fluid_boundary = fluid1_fluid2_boundary | fluid2_fluid1_boundary
    all_boundary = fluid_solid_boundary | fluid_fluid_boundary
    # -------------------------------------------------------------------------- Find isolated fluid clusters
    domains1, max_label1 = bwlabel3D_set(cc, geo, (1,))
    domains2, max_label2 = bwlabel3D_set(cc, geo, (2,))
    domains = domains1 + domains2 + (domains2 > 0)*max_label1
    del domains1, domains2
    max_label = max_label1 + max_label2
    # -------------------------------------------------------------------------- Remove cluster less than a given size
    for ij in all_boundary:
        geo[ij] = -geo[ij]
    for n in range(1,max_label+1):
        if np.sum(geo[domains==n] > 0) < cluster_cutoff_size:
            geo[domains==n] = 0
    for ij in all_boundary:
        geo[ij] = -geo[ij]

    return geo

@njit
def gradient_linear_extended(sd):
    """
    Calculate the gradient based on the signed distance function. Extends sd 
    with an extra outer periodic layer that is a reflection around the system 
    boundaries to simplify the calculation of gradients

    Parameters
    ----------
    sd : numpy.ndarray[float]
        3D scalar array 

    Return
    -----
    : numpy.ndarray[float]
        gradient of the input array 
    """
    I, J, K = sd.shape
    I += 2
    J += 2
    K += 2
    a = np.zeros((I, J, K), dtype=sd.dtype)
    a[1:-1, 1:-1, 1:-1] = sd[:,:,:]
    # -------------------------------------------------------------------------- Periodic boundaries
    # sides 
    a[0, :, :]  = 2*a[1, :, :] - a[2, :, :]
    a[:, 0, :]  = 2*a[:, 1, :] - a[:, 2, :]
    a[:, :, 0]  = 2*a[:, :, 1] - a[:, :, 2]
    a[-1, :, :]  = 2*a[-2, :, :] - a[-3, :, :]
    a[:, -1, :]  = 2*a[:, -2, :] - a[:, -3, :]
    a[:, :, -1]  = 2*a[:, :, -2] - a[:, :, -3]

    grad = np.zeros((3,) + sd.shape)
    grad[0] = 0.5*(a[2:, 1:-1, 1:-1] - a[:-2, 1:-1, 1:-1])
    grad[1] = 0.5*(a[1:-1, 2:, 1:-1] - a[1:-1, :-2, 1:-1])
    grad[2] = 0.5*(a[1:-1, 1:-1, 2:] - a[1:-1, 1:-1, :-2])

    return grad

def calculate_normals(sd):
    """
    Calculate the normals based on a signed distance function.
    The boundaries are treated with linear extension. The normals
    are set to zero if the gradient of the sdf is zero

        Parameters
    ----------
    sd : numpy.ndarray[float]
        3D scalar array 

    Return
    -----
    : numpy.ndarray[float]
        array of normals 
    """
    grad = gradient_linear_extended(sd)
    length = np.sqrt(np.sum(grad**2, axis=0))
    length[length == 0] = 1
    normal = grad/length
    return normal

@njit
def extend_geometry_periodic(geo, flow_direction):
    """
    Extends the geometry with an extra outer periodic layer to simplify the boundary 
    recognition algorithm. Pressure boundaries are added for the inlet
    and outlet in the flow direction, according to the following label 
    scheme:
        SOLID        : 0
        PHASE 1      : 1
        PHASE 2      : 2
        PRESSURE BND : 3

    It is assumed that that system boundaries, other that the inlet and 
    outlet, are periodic. 

    Parameters
    ----------
    geo : numpy.ndarray[int]
        It is assumed that geo matrix that that contains the markers:
        SOLID : 0
        FLUID PHASES : 1 or 2

    flow_direction: int
        indicates the flow direction    

    Return
    ------
    : numpy.ndarray[int]
        Extended marked geometry with size = geo.size + 2        
    """
    # ========================================================================== Extend geometry
    I, J, K = geo.shape
    I += 2
    J += 2
    K += 2
    a = np.zeros((I, J, K), dtype=geo.dtype)
    a[1:-1, 1:-1, 1:-1] = geo[:,:,:]
    # -------------------------------------------------------------------------- Periodic boundaries
    # sides 
    for il, ir in [(0, -2), (-1, 1)]:
        a[il, :, :] = a[ir, :, :]
        a[:, il, :] = a[:, ir, :]
        a[:, :, il] = a [:, :, ir]
        # edges
        for jl, jr in [(0, -1), (-2, 1)]:
            a[il, jl, :] = a[ir, jr, :]
            a[il, :, jl] = a[ir, :, jr]
            a[:, il, jr] = a[:, il, jr]
            # corners
            for kl, kr in [(0, -1), (-2, 1)]:
                a[il, jl, kl] = a[ir, jr, kr]
    # ========================================================================== Set pressure boundaries
    if flow_direction == 0:
        a[0, :, :] = 3
        a[-1, :, :] = 3
    elif flow_direction == 1:
        a[:, 0, :] = 3
        a[:, -1, :] = 3
    elif flow_direction == 2:
        a[:, :, 0] = 3
        a[:, :,-1] = 3
    return a

@njit
def extend_geometry_solid(geo, flow_direction):
    """
    Extends the geometry with an extra outer solid layer to simplify the boundary 
    recognition algorithm. Pressure boundaries are added for the inlet
    and outlet in the flow direction, according to the following label 
    scheme:
        SOLID        : 0
        PHASE 1      : 1
        PHASE 2      : 2
        PRESSURE BND : 3

    It is assumed that that system boundaries, other that the inlet and 
    outlet, are periodic. 

    Parameters
    ----------
    geo : numpy.ndarray[int]
        It is assumed that geo matrix that that contains the markers:
        SOLID : 0
        FLUID PHASES : 1 or 2

    flow_direction: int
        indicates the flow direction    

    Return
    ------
    : numpy.ndarray[int]
        Extended marked geometry with size = geo.size + 2        
    """
    # ========================================================================== Extend geometry
    I, J, K = geo.shape
    I += 2
    J += 2
    K += 2
    geo_extended = np.zeros((I, J, K), dtype=geo.dtype) # Sets all nodes to solid
    geo_extended[1:-1, 1:-1, 1:-1] = geo[:,:,:] # Fills in the interior of the extended geo

    # ========================================================================== Set pressure boundaries
    if flow_direction == 0:
        geo_extended[0, :, :] = 3
        geo_extended[-1, :, :] = 3
    elif flow_direction == 1:
        geo_extended[:, 0, :] = 3
        geo_extended[:, -1, :] = 3
    elif flow_direction == 2:
        geo_extended[:, :, 0] = 3
        geo_extended[:, :,-1] = 3
    return geo_extended

@njit
def mark_boundaries(cc, geo_extended):
    """
    Marks the boundaries in the interior of 'geo__extended'. 
    Assumes that 'geo__extended' has the markers:
         0: solid
         1: fluid phase 1
         2: fluid phase 2 
         3: pressure boundary

    The marked geometry has the additional markers
         4: fluid-fluid interface
         8: fluid-solid interface
        16: pressure boundary

    Parameters
    ----------
    cc : numpy.ndarray[int]
        List of lattice directions with shape (Q, D). Defines 
        the connectivity of the clusters

    geo_extended : numpy.ndarray[int]
        It is assumed that geo_extended matrix contains the markers:
        SOLID : 0
        FLUID PHASES : 1 or 2
        PRESSURE BOUNDARY : 3

    Return
    ------
    : numpy.ndarray[int]
        Boundary marked geometry
    """ 
    I, J, K = geo_extended.shape    
    geo_marked = np.copy(geo_extended)
    for i in range(1, I-1):
        for j in range(1, J-1):
            for k in range(1, K-1):
                val = geo_extended[i, j, k]
                fluidbnd = 0
                solidbnd = 0
                pressurebnd = 0
                if val in [1, 2]:
                    for c in cc:
                        ni = i + c[0]
                        nj = j + c[1]
                        nk = k + c[2]
                        nval = geo_extended[ni, nj, nk]
                        # fluid-fluid (4)
                        if (nval in [1, 2]) and (nval != val):
                            fluidbnd = 4
                        # fluid-solid (8)
                        elif nval == 0:
                            solidbnd = 8
                        # pressure (16)
                        elif nval == 3:
                            pressurebnd = 16
                geo_marked[i, j, k] = val + fluidbnd + solidbnd + pressurebnd
    return geo_marked

def mark_boundaries_solid(cc, geo, flow_direction):
    """
    Marks the boundaries in 'geo'. Assumes that 'geo' has the markers:
         0: solid
         1: fluid phase 1
         2: fluid phase 2 

    Note the geometry is extended by a solid rim in all directions.

    The marked geometry has the additional markers
         3: pressure boundary' ghost node (treated as a solid node)
         4: fluid-fluid interface
         8: fluid-solid interface
        16: pressure boundary

    The pressure boundary is given by the inlet and outlet in the
    geometry. This are set based on the flow_direction, that 
    defines the surface normal to the inlet and outlet.

    Parameters
    ----------
    cc : numpy.ndarray[int]
        List of lattice directions with shape (Q, D). Defines 
        the connectivity of the clusters

    geo : numpy.ndarray[int]
        It is assumed that geo matrix that that contains the markers:
        SOLID : 0
        FLUID PHASES : 1 or 2

    flow_direction: int [0, 1, 2]
        indicates the flow direction    

    Return
    ------
    : numpy.ndarray[int]
        Boundary marked geometry. The size of the marked geometry
        is geo.shape + 2. 

    """ 
    geo_extended = extend_geometry_solid(geo, flow_direction)
    geo_marked = mark_boundaries(cc, geo_extended)
    return geo_marked

def find_fluid_phase_clusters(cc, geo_marked, flow_direction):
    """
    Find fluid phase clusters and marks as:
      1 : Phase clusters percolating in the flow direction.

    Parameters
    ----------
    cc : numpy.ndarray[int]
        List of lattice directions with shape (Q, D). Defines 
        the connectivity of the clusters

    geo_marked : numpy.ndarray[int]
        Assumed to be an array returned by `mark_boundaries_solid`.

    flow_direction: int [0, 1, 2]
        indicates the flow direction    

    Return
    ------
    : numpy.ndarray[int], int
        Cluster marked geometry. Same size as geo_marked and 
        maximum label marker  
    """
    # Find all combinations of boundary marks that can be added to
    # the fluid phase markers 1 and 2
    bvals = tuple((sum((x*y for x, y in zip(xx, (4, 8, 16)))) for xx in product([0, 1], repeat=3)))
    set_of_label1 = tuple(1 + x for x in bvals)
    set_of_label2 = tuple(2 + x for x in bvals)
    domains1, max_label1 = bwlabel3D_set(cc, geo_marked, set_of_label1)
    domains2, max_label2 = bwlabel3D_set(cc, geo_marked, set_of_label2)
    domains = domains1 + domains2 + (domains2 > 0)*max_label1
    max_label = max_label1 + max_label2    
    return domains, max_label

if __name__ == "__main__":
    # ########################################################################## Load data
    # ========================================================================== Data folders
    data_path = "/home/AD.NORCERESEARCH.NO/esje/Programs/Python/LB/Data/"
    output_path = "/home/AD.NORCERESEARCH.NO/esje/Programs/Python/CSSR/RelPerm/Data/"
    # -------------------------------------------------------------------------- Solid
    # pore >= 0 : SOLID
    # pore  < 0 : FLUID
    pore =  np.load(data_path + 'LV60A_PoreSolid_200x200x300_SDF.npy')
    # -------------------------------------------------------------------------- Fluid
    # fluid >= 0 : FLUID-PHASE 1
    # fluid  < 0 : FLUID-PHASE 2
    fluids = np.load(data_path + "LV60A_NWP_Sw062_200x200x300_SDF.npy") 
    # -------------------------------------------------------------------------- geo
    # SOLID         : 0
    # FLUID PHASE 1 : 1
    # FLUID PHASE 2 : 2
    geo = np.ones(pore.shape, dtype=np.int32) 
    geo[fluids < 0] = 2 # mark fluid phase 2
    geo[pore >= 0] = 0 # mark solid
    # ########################################################################## CLEAN GEOMETRY
    # ========================================================================== Define lattice
    # -------------------------------------------------------------------------- D3Q7
    cc = np.ones((7, 3), dtype=np.int64)
    cc[:4, :] = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [0, 0, 0]
        ]
    cc[:-4:-1, :] = -cc[:3, :]
    cc_d3q7 = np.copy(cc)
    # -------------------------------------------------------------------------- D3Q19
    cc = np.ones((19, 3), dtype=np.int64)
    cc[:10, :] = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 0],
        [1, 0, 1],
        [0, 1, 1],
        [-1, 1, 0],
        [-1, 0, 1],
        [0, -1, 1],
        [0, 0, 0]
        ]
    cc[:-10:-1, :] = -cc[:9, :]
    cc_d3q19 = np.copy(cc)
    del cc
    # ========================================================================== Test clean_geometry
    from time import time_ns

    print("TEST clean geometry:")

    t0 = time_ns()
    
    while True:
        print("inside clean geo")
        geo_run = clean_geometry(cc_d3q7, geo)
        if np.all(geo_run == geo):
            break
        geo[:] = geo_run[:]
    t1 = time_ns()
    print("    Run time for clean_geometry: ", (t1-t0)*1.0e-9)

    t0 = time_ns()
    geo_marked = mark_boundaries(cc_d3q19, geo_run, 2)
    t1 = time_ns()
    print("    Run time for mark_boundaries", (t1-t0)*1.0e-9)

    # ========================================================================== Plot clean_geometry
    import pyvista as pv

    # Plot domains
    grid = pv.ImageData()
    grid.dimensions = geo_marked.shape
    grid.origin = (0, 0, 0)
    grid.spacing = (1, 1, 1)
    grid.point_data["geo"] = geo_marked.flatten(order="F")
    grid.plot(show_edges=False)


