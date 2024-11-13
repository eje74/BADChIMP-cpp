import numpy as np
from matplotlib import pyplot as plt
from numba import njit
import vtklb as vlb
from itertools import product

# ################################################################################
#                             FUNCTIONS                                          #
# ################################################################################
# ================================================================================ Mask polygon
@njit
def jit_in_polygon(point, line_segments, line_segments_right_shift):
    """
    Will return true if the point is within or on the polygon defined 
    by the line_segments and false otherwise.
    Assumes the line_segments contains the points that defines the 
    polygon, and that they are order in the counter clockwise
    direction. 
      The shape of line segments is assumed to be on the 
    following form:
        numpy.array([[x0, y0], [x1, y1], ..., [xn, yn]])
    
    """
    res = 0
    for right, left in zip(line_segments, line_segments_right_shift):
        a = right - point
        b = left - point
        cross = a[0]*b[1] - a[1]*b[0]
        dot = np.sum(a*b) 
        if cross == 0 and dot <= 0:
            return 2*np.pi
        a_norm = np.linalg.norm(a)
        b_norm = np.linalg.norm(b)
        if a_norm == 0 or b_norm == 0:
            print(a_norm, b_norm)
        cross /= np.linalg.norm(a)*np.linalg.norm(b)
        cross = min(cross, 1)
        cross = max(cross, -1)
        theta0 = np.arcsin(cross)
        theta = (( 2*(theta0>0) - 1)*np.pi - 2*theta0)*(dot < 0) + theta0
        res += theta
    return res

def in_polygon(point, line_segments):
    line_segments_right_shift = np.roll(line_segments, -1, axis=0)
    return jit_in_polygon(point, line_segments, line_segments_right_shift)

@njit
def jit_mask_polygon(bounding_box, II, JJ, line_segments, line_segments_right_shift):
    mask = np.ones(II.shape, dtype=np.int32)
    for i, j in bounding_box:
        x = II[i, j]
        y = JJ[i, j]
        val = jit_in_polygon(np.array([x, y]), line_segments, line_segments_right_shift) 
        mask[i, j] = val < 3   
    return mask
# -------------------------------------------------------------------------------- mask_polygon

def mask_polygon(II:np.ndarray, JJ:np.ndarray, line_segments:np.ndarray):
    """
    Will generate a mask matrix where nodes within the polygon defined 
    by line_segments are given the value 0 while nodes outside is given 
    the value 1.

    Parameters
    ----------
    II : numpy array
        Array of component of the 0'th position component, i.e. x-component, 
        such that II[i+1, j] = II[i, j] + dx
    JJ : numpy array
        Array of component of the 1'st position component, i.e. y-component,
        such that JJ[i, j+1] = JJ[i, j] + dy
    line_segments : numpy array
        Assumes the line_segments contains the points that defines the 
        polygon, and that they are order in the counter clockwise direction. 
        The shape of line segments is assumed to be on the following form:
        numpy.array([[x0, y0], [x1, y1], ..., [xn, yn]])

    Returns
    -------
    numpy array:
        Array with the same shape as 'II', containing '1's and '0's.
    """
    # Find the bounding box for the polygon
    i_min, j_min = np.min(line_segments, axis=0)
    i_max, j_max = np.max(line_segments, axis=0)
    indics_within_bounding_box = np.argwhere((i_min <= II) & (II <= i_max) & (j_min <= JJ) & (JJ <= j_max))
    line_segments_right_shift = np.roll(line_segments, -1, axis=0)
    return jit_mask_polygon(indics_within_bounding_box, II, JJ, line_segments, line_segments_right_shift)

@njit
def markBoundaryNodes(geo:np.ndarray, basis):
    """
    Marks all boundary nodes

    Parameters
    ----------
    geo: array
        geo > 0 fluid and geo == 0 solid
    basis: array
        define neighborhood. List of di and dj pairs.

    Returns
    -------
    np.array
        Boundary nodes = 1, other nodes = 0 
    """
    I, J = geo.shape
    bnd = np.zeros_like(geo)
    for i in range(1, I-1):
        for j in range(1, J-1):
            if geo[i, j] > 0:
                for di, dj in basis:
                    if (geo[i+di, j+dj] == 0):
                        bnd[i, j] = 1
    return bnd

@njit
def segment_intersection(sega, segb):
    """
    The function finds the point of intersection between 
    the line segments sega and segb. The position (x, y) is given 
    by the parameters s, t so that

        x = sega[0][0]*(1-s) + sega[1][0]*s
        y = sega[0][1]*(1-s) + sega[1][1]*s

                    or

        x = segb[1][0]*(1-t) + segb[0][0]*t
        y = segb[1][1]*(1-t) + segb[0][1]*t
                    
    
    if s,t > 1 or s,t < 0 the segments do not intersect.    

    We assume that seg = np.array[[x1, y1],[x2, y2]]

    Parameters
    ----------
    sega : array
        A line segment defined as np.array[[x1, y1],[x2, y2]]
    segb : array
        A line segment defined as np.array[[x1, y1],[x2, y2]]

    Returns
    -------
    tuple
        (s, t) if s,t > 1 or s,t < 0 the segments do not intersect.
    """
    a = sega[1] - sega[0]
    b = segb[1] - segb[0]

    c = segb[1] - sega[0] # right hand side
    detA = a[0]*b[1] - a[1]*b[0]

    s = -1
    t = -1

    if np.abs(detA) > 0:
        detA1 = c[0]*b[1] - c[1]*b[0]
        detA2 = a[0]*c[1] - a[1]*c[0]
        s = detA1/detA
        t = detA2/detA

    return s, t

@njit
def surface_length(r, n, gamma, pnts):
    t = np.array([-n[1], n[0]])
    a0 = np.array(2)
    a1 = np.array(2)
    cnt = 0
    for i in range(4):
        ip = i % 4
        p2 = pnts[ip]
        p1 = pnts[i]
        dp = p2 - p1
        detA = -dp[0]*t[1] + dp[1]*t[0]
        if (gamma < 1e-8) and (np.abs(detA) < 1e-8):
            return 1.0 
        b = r - p1 - gamma*n
        detAb = -b[0]*t[1] + b[1]*t[0]
        beta = detAb/detA
        if (0 <= beta) and (beta < 1):
            if cnt == 0:
                a0[:] = beta*dp + p1
                cnt += 1
            elif cnt == 1:
                a1[:] = beta*dp + p1
                return np.sqrt( (a1-a0)**2)
            else:
                print(f"surface_length, error, n={n}")
                exit(0)
    print(f"surface_length, error could not find a length")
    exit(1)

@njit
def link_parameters(seg, line_segments):
    """
    Returns parameters for the link defined seg:
    - gamma: relative distance from fluid node to boundary
    - nx : x-component of the outward pointing normal vector
    - ny : y-component of the outward pointing normal vector

    Parameters
    ----------
    seg: numpy array
        A line segment defined as np.array[[x1, y1],[x2, y2]]    
    line_segments : numpy array
        Assumes the line_segments contains the points that defines the 
        polygon, and that they are order in the counter clockwise direction. 
        The shape of line segments is assumed to be on the following form:
        numpy.array([[x0, y0], [x1, y1], ..., [xn, yn]])     
    
    Returns
    -------
    list
        gamma, nx, ny
    """
    N, D = line_segments.shape
    sb = np.zeros((2, 2), dtype=np.float64)
    test = False
    for nl in range(N):
        nr = (nl + 1) % N
        sb[0, :] = line_segments[nl, :]
        sb[1, :] = line_segments[nr, :]
        s, t = segment_intersection(seg, sb)
        if 0 <= s and s <= 1 and 0 <= t and t <= 1:
            tangent = sb[1]-sb[0]
            nx = tangent[1]
            ny = -tangent[0]
            tmp = np.sqrt(nx**2 + ny**2)
            nx = nx/tmp
            ny = ny/tmp
            test = True
            return s, nx, ny
    return -1, -1, -1

@njit 
def interpolation_weights(x, y, interpolation_points):
    """
    Find the weights for interpolation/extrapolation for the point (x, y).
    It is assumed that the ```interpolation_points``` array is on the form:
              [[x0, y0], [x1, y1], [x2, y2], [x3, y3]],
    and that the relative positioning of the points are give as

                 (0, 1) 3---------2 (1, 1)
                        |         | 
                        |         |
                        |         |
                 (0, 0) 0---------1 (1, 0)

    Parameters
    ----------
    x, y : floats
        x- and y-coordinate of the point you want interpolate
    interpolation_points: array
        position of the four interpolation points

    Returns
    -------
    array:
        array of weights, w, for each of the interpolation_points so that
        
            f(x, y) = \sum_i w[i]*f(xi, yi)

    """
    w = np.zeros(4)

    x0, y0 = interpolation_points[0]
    x1, y1 = interpolation_points[1]
    x2, y2 = interpolation_points[2]
    x3, y3 = interpolation_points[3]

    w[0] = (x1 - x)*(y3 - y)
    w[1] = (x - x0)*(y2 - y)
    w[2] = (x - x3)*(y - y1)
    w[3] = (x2 - x)*(y - y0)

    return w

@njit 
def interpolation_weights_2nd_order(xp, yp, interpolation_points):
    """
    Find the weights for interpolation/extrapolation for the point (x, y).
    It is assumed that the ```interpolation_points``` array is on the form:
              [[x0, y0], [x1, y1], [x2, y2], [x3, y3], ..., [x8, y8]],
    and that the relative positioning of the points are give as

            (0,  2) 6---------7---------8 (2,  2)
                    |       (1, 2)      |
                    |         |         |
                    |         |         |
            (0,  1) 5---------4---------3 (2, 1)
                    |       (1, 1)      |
                    |         |         |
                    |         |         |
            (0, 0)  0---------1---------2 (2, 0)
                            (1, 0)
    Parameters
    ----------
    xp, yp : floats
        x- and y-coordinate of the point you want interpolate
    interpolation_points: array
        position of the four interpolation points

    Returns
    -------
    array:
        array of weights, w, for each of the interpolation_points so that
        
            f(x, y) = \sum_i w[i]*f(xi, yi)

    """
    w = np.zeros(9)
    x0, y0 = interpolation_points[0]
    x = 0.5*(xp-x0)
    y = 0.5*(yp-y0)

    w[0] = (1-x)*(1-2*x)*(1-y)*(1-2*y)
    w[1] = 4*x*(1-x)*(1-y)*(1-2*y)
    w[2] = x*(2*x-1)*(1-y)*(1-2*y)
    w[3] = 4*x*(2*x-1)*y*(1-y)
    w[4] = 16*x*(1-x)*y*(1-y)
    w[5] = 4*(1-x)*(1-2*x)*y*(1-y)
    w[6] = (1-x)*(1-2*x)*y*(2*y-1)
    w[7] = 4*x*(1-x)*y*(2*y-1)
    w[8] = x*(2*x-1)*y*(2*y-1)
    
    return w

@njit
def jit_write_boundary_file(boundary_nodes, interpolation_points_center, interpolation_points, line_segments, basis):
    N, D = boundary_nodes.shape
    Q, D = basis.shape
    boundary_normal = np.zeros(2)
    boundary_data = np.zeros((Q*N, 18))
    cnt = 0
    for n, pos in enumerate(boundary_nodes):
        for q, c in enumerate(basis):
            #------------------------------------- position
            boundary_data[cnt, 0] = pos[0] # x 
            boundary_data[cnt, 1] = pos[1] # y 
            #------------------------------------- unknown direction
            boundary_data[cnt, 2] = q
            #------------------------------------- intersection
            neig_pos = pos - c
            seg = np.array([[pos[0], pos[1]], [neig_pos[0], neig_pos[1]]])
            gamma, nx, ny = link_parameters(seg ,line_segments)
            boundary_data[cnt, 3] = gamma
            #------------------------------------- interpolation point
            boundary_normal[0] = nx
            boundary_normal[1] = ny
            pos_bnd = pos - gamma*c
            if gamma < 0 :
                ind = -1
            else:
                ind = find_nearest_interpolation_point(pos_bnd, boundary_normal, interpolation_points_center)
                points = interpolation_points[ind]
                boundary_data[cnt, 4:6] = points[0]
                boundary_data[cnt, 6:8] = points[1]
                boundary_data[cnt, 8:10] = points[2]
                boundary_data[cnt, 10:12] = points[3]
                #------------------------------------- interpolation weights
                w = interpolation_weights(pos_bnd[0], pos_bnd[1], points)
                boundary_data[cnt, 12:16] = w
                #------------------------------------- boundary normal
                boundary_data[cnt, 16] = nx
                boundary_data[cnt, 17] = ny
            cnt += 1 
    return boundary_data
    
def write_boundary_file(file_name, boundary_nodes, interpolation_points, line_segments, basis):
    """
    Parameters
    ----------
    file_name: string
        File name 
    boundary_node: array
        Array of boundary positions [[i0, j0], [i1, j1], ...]
    interpolation_points: array
        Output from ```find_interpolation_points```        
    line_segments : numpy array
        Assumes the line_segments contains the points that defines the 
        polygon, and that they are order in the counter clockwise direction. 
        The shape of line segments is assumed to be on the following form:
        numpy.array([[x0, y0], [x1, y1], ..., [xn, yn]])     
    basis: array
        Array defining the local neighborhood on the form 
        [[c0x, c0y], [c1x, c1y], ...] 

    Returns
    -------
    np.array
        One entry for each boundary nodes' basis direction: 
        position               |  x, y    |  0,  1
        unknown direction      |     q    |  2
        distance from boundary |  gamma   |  3
        point0                 | (x0, y0) |  4,  5
        point1                 | (x1, y1) |  6,  7
        point2                 | (x2, y2) |  8,  9
        point3                 | (x3, y3) | 10, 11
        weight p0              |    w0    | 12
        weight p1              |    w1    | 13
        weight p2              |    w2    | 14
        weight p3              |    w3    | 15
        wall normal            | nx, ny   | 16, 17
    """
    interpolation_points_center = np.mean(interpolation_points, axis=1)
    data = jit_write_boundary_file(boundary_nodes, interpolation_points_center, np.array(interpolation_points, dtype=np.int64), line_segments, basis)
    np.savetxt(file_name, data, fmt="%d %d %d %.15f" + " %d"*8 + " %.15f"*4 + " %.15f"*2)
    return data

@njit
def jit_write_boundary_file_node(boundary_nodes, interpolation_points_center, interpolation_points, line_segments, basis):
    N, D = boundary_nodes.shape
    Q, D = basis.shape
    boundary_normal = np.zeros(2)
    boundary_data = np.zeros((N, 26))
    cnorm = np.zeros(Q)
    for q, c in enumerate(basis):
        cnorm[q] = np.sqrt(np.sum(c**2))
    for bndno, pos in enumerate(boundary_nodes):
        min_q = -1
        min_dist = 1e10
        for q, c in enumerate(basis):
            #------------------------------------- link intersections
            neig_pos = pos - c
            seg = np.array([[pos[0], pos[1]], [neig_pos[0], neig_pos[1]]])
            gamma, nx, ny = link_parameters(seg ,line_segments)
            if gamma > -1:
                dist =  gamma*np.sqrt(np.sum(c**2)) #   np.linalg.norm(gamma*c)
                if dist < min_dist:
                    boundary_normal[0] = nx
                    boundary_normal[1] = ny
                    min_q = q
        if min_q == -1:
            print("ERROR (jit_write_boundary_file_node): Boundary do not cross the boundary!")
        else:
            #------------------------------------- node intersection
            neig_pos = pos - boundary_normal*cnorm[min_q]
            seg = np.zeros((2, 2))
            seg[0, 0] = pos[0]
            seg[0, 1] = pos[1]
            seg[1, 0] = neig_pos[0]
            seg[1, 1] = neig_pos[1]
            gamma, nx, ny = link_parameters(seg ,line_segments)
            if gamma > -1:
                boundary_normal[0] = nx
                boundary_normal[1] = ny
            else:
                print("ERROR (jit_write_boundary_file_node): Boundary do not cross the boundary in its normal direction!")
            gamma = gamma*cnorm[min_q]
            boundary_data[bndno, 0] = pos[0]
            boundary_data[bndno, 1] = pos[1]
            boundary_data[bndno, 2] = gamma
            #------------------------------------- find interpolation points
            ipsind = find_nearest_interpolation_point(pos, boundary_normal, interpolation_points_center)
            ippos = interpolation_points[ipsind]
            for n, i in enumerate((3, 5, 7, 9)):
                boundary_data[bndno, i:i+2] = ippos[n]
            ippos_mean = 0.25*np.sum(ippos, axis=0)
            ips_gamma = np.dot(ippos_mean-pos, boundary_normal)
            boundary_data[bndno, 11] = ips_gamma
            #------------------------------------- find interpolation weights boundary
            x = pos[0] - gamma*boundary_normal[0]
            y = pos[1] - gamma*boundary_normal[1]
            w = interpolation_weights(x, y,ippos)
            boundary_data[bndno, 12:16] = w
            #------------------------------------- find interpolation weights interpolation points
            x = pos[0] + ips_gamma*boundary_normal[0]
            y = pos[1] + ips_gamma*boundary_normal[1]
            w = interpolation_weights(x, y,ippos)
            boundary_data[bndno, 16:20] = w
            #------------------------------------- find interpolation weights node
            x = pos[0]
            y = pos[1]
            w = interpolation_weights(x, y,ippos)
            boundary_data[bndno, 20:24] = w
            #------------------------------------- boundary normal
            boundary_data[bndno, 24] = boundary_normal[0]
            boundary_data[bndno, 25] = boundary_normal[1]
            #------------------------------------- surface length
            print(surface_length(pos, boundary_normal, gamma, ippos))

    return boundary_data

def write_boundary_file_node(file_name, boundary_nodes, interpolation_points, line_segments, basis):
    """
    Parameters
    ----------
    file_name: string
        File name 
    boundary_node: array
        Array of boundary positions [[i0, j0], [i1, j1], ...]
    interpolation_points: array
        Output from ```find_interpolation_points```        
    line_segments : numpy array
        Assumes the line_segments contains the points that defines the 
        polygon, and that they are order in the counter clockwise direction. 
        The shape of line segments is assumed to be on the following form:
        numpy.array([[x0, y0], [x1, y1], ..., [xn, yn]])     
    basis: array
        Array defining the local neighborhood on the form 
        [[c0x, c0y], [c1x, c1y], ...] 

    Returns
    -------
    np.array
        One entry for each boundary node
        position               |  x, y    |  0,  1
        distance from boundary |  gamma   |  2
        point0                 | (x0, y0) |  3,  4
        point1                 | (x1, y1) |  5,  6
        point2                 | (x2, y2) |  7,  8
        point3                 | (x3, y3) |  9, 10
        distance to cm interpol|  gamma2  | 11
        weight bnd intersection|  w_a[:4] | 12, 13, 14, 15
        weight interp. points  |  w_b[:4] | 16, 17, 18, 19
        weight boundary node   |  w_c[:4] | 20, 21, 22, 23
        wall normal            | nx, ny   | 24, 25

    Notes
    -----
    Explanatory comments to the return values
    - position: is the x and y coordinate to the boundary node. Is also
       used for identification.
    - distance from boundary: -gamma*normal is the nearest point
       on the surface relative to 'position'.
    - point0-3: gives the coordinates to the four points used for the 
       bi-linear interpolation
    - distance to cm interpol: gamma2*normal is the point nearest to the
       center of the four interpolation points along the normal centered 
       at 'position'.
    - weight bnd intersection: interpolation weights such that
             sum w_a[n]*point_n = (x, y) - gamma*normal
    - weight interp. points: interpolation weights such that
             sum w_b[n]*point_n = (x, y) + gamma2*normal
    - weight boundary node: interpolation weights such that
             sum w_c[n]*point_n = (x, y)
    - wall normal: is the outward pointing wall normal
    """
    interpolation_points_center = np.mean(interpolation_points, axis=1)
    data = jit_write_boundary_file_node(boundary_nodes, interpolation_points_center, np.array(interpolation_points, dtype=np.int64), line_segments, basis)
    np.savetxt(file_name, data, fmt="%d %d" + " %.15f" + " %d"*8 + " %.15f"*(1 + 4 + 4 + 4) + " %.15f"*2)
    return data

@njit
def jit_write_boundary_file_node_2nd_order(boundary_nodes, interpolation_points_center, interpolation_points, line_segments, basis):
    N, D = boundary_nodes.shape
    Q, D = basis.shape
    boundary_normal = np.zeros(2)
    boundary_data = np.zeros((N, 26))
    cnorm = np.zeros(Q)
    for q, c in enumerate(basis):
        cnorm[q] = np.sqrt(np.sum(c**2))
    for bndno, pos in enumerate(boundary_nodes):
        min_q = -1
        min_dist = 1e10
        for q, c in enumerate(basis):
            #------------------------------------- link intersections
            neig_pos = pos - c
            seg = np.array([[pos[0], pos[1]], [neig_pos[0], neig_pos[1]]])
            gamma, nx, ny = link_parameters(seg ,line_segments)
            if gamma > -1:
                dist =  gamma*np.sqrt(np.sum(c**2)) #   np.linalg.norm(gamma*c)
                if dist < min_dist:
                    boundary_normal[0] = nx
                    boundary_normal[1] = ny
                    min_q = q
        if min_q == -1:
            print("ERROR (jit_write_boundary_file_node): Boundary do not cross the boundary!")
        else:
            #------------------------------------- node intersection
            neig_pos = pos - boundary_normal*cnorm[min_q]
            seg = np.zeros((2, 2))
            seg[0, 0] = pos[0]
            seg[0, 1] = pos[1]
            seg[1, 0] = neig_pos[0]
            seg[1, 1] = neig_pos[1]
            gamma, nx, ny = link_parameters(seg ,line_segments)
            if gamma > -1:
                boundary_normal[0] = nx
                boundary_normal[1] = ny
            else:
                print("ERROR (jit_write_boundary_file_node): Boundary do not cross the boundary in its normal direction!")
            gamma = gamma*cnorm[min_q]
            boundary_data[bndno, 0] = pos[0]
            boundary_data[bndno, 1] = pos[1]
            boundary_data[bndno, 2] = gamma
            #------------------------------------- find interpolation points
            ipsind = find_nearest_interpolation_point(pos, boundary_normal, interpolation_points_center)
            
            ippos = interpolation_points[ipsind]
            for n, i in enumerate((3, 5, 7, 9, 11, 13, 15, 17, 19)):
                boundary_data[bndno, i:i+2] = ippos[n]
            ippos_mean = np.sum(ippos, axis=0)/9
            ips_gamma = np.dot(ippos_mean-pos, boundary_normal)
            boundary_data[bndno, 21] = ips_gamma
            # #------------------------------------- find interpolation weights boundary
            # x = pos[0] - gamma*boundary_normal[0]
            # y = pos[1] - gamma*boundary_normal[1]
            # w = interpolation_weights(x, y,ippos)
            # boundary_data[bndno, 12:16] = w
            # #------------------------------------- find interpolation weights interpolation points
            # x = pos[0] + ips_gamma*boundary_normal[0]
            # y = pos[1] + ips_gamma*boundary_normal[1]
            # w = interpolation_weights(x, y,ippos)
            # boundary_data[bndno, 16:20] = w
            # #------------------------------------- find interpolation weights node
            # x = pos[0]
            # y = pos[1]
            # w = interpolation_weights(x, y,ippos)
            # boundary_data[bndno, 20:24] = w
            # #------------------------------------- boundary normal
            # boundary_data[bndno, 24] = boundary_normal[0]
            # boundary_data[bndno, 25] = boundary_normal[1]
    return boundary_data

def write_boundary_file_node_2nd_order(file_name, boundary_nodes, interpolation_points, line_segments, basis):
    """
    Parameters
    ----------
    file_name: string
        File name 
    boundary_node: array
        Array of boundary positions [[i0, j0], [i1, j1], ...]
    interpolation_points: array
        Output from ```find_interpolation_points```        
    line_segments : numpy array
        Assumes the line_segments contains the points that defines the 
        polygon, and that they are order in the counter clockwise direction. 
        The shape of line segments is assumed to be on the following form:
        numpy.array([[x0, y0], [x1, y1], ..., [xn, yn]])     
    basis: array
        Array defining the local neighborhood on the form 
        [[c0x, c0y], [c1x, c1y], ...] 

    Returns
    -------
    np.array
        One entry for each boundary node
        position               |  x, y    |  0,  1
        distance from boundary |  gamma   |  2
        point0                 | (x0, y0) |  3,  4
        point1                 | (x1, y1) |  5,  6
        point2                 | (x2, y2) |  7,  8
        point3                 | (x3, y3) |  9, 10
        point4                 | (x4, y4) | 11, 12
        point5                 | (x5, y5) | 13, 14
        point6                 | (x6, y6) | 15, 16
        point7                 | (x7, y7) | 17, 18
        point8                 | (x8, y8) | 19, 20
        distance to cm interpol|  gamma2  | 21
        weight bnd intersection|  w_a[:4] | 12, 13, 14, 15
        weight interp. points  |  w_b[:4] | 16, 17, 18, 19
        weight boundary node   |  w_c[:4] | 20, 21, 22, 23
        wall normal            | nx, ny   | 24, 25

    Notes
    -----
    Explanatory comments to the return values
    - position: is the x and y coordinate to the boundary node. Is also
       used for identification.
    - distance from boundary: -gamma*normal is the nearest point
       on the surface relative to 'position'.
    - point0-3: gives the coordinates to the four points used for the 
       bi-linear interpolation
    - distance to cm interpol: gamma2*normal is the point nearest to the
       center of the four interpolation points along the normal centered 
       at 'position'.
    - weight bnd intersection: interpolation weights such that
             sum w_a[n]*point_n = (x, y) - gamma*normal
    - weight interp. points: interpolation weights such that
             sum w_b[n]*point_n = (x, y) + gamma2*normal
    - weight boundary node: interpolation weights such that
             sum w_c[n]*point_n = (x, y)
    - wall normal: is the outward pointing wall normal
    """
    interpolation_points_center = np.mean(interpolation_points, axis=1)
    data = jit_write_boundary_file_node_2nd_order(boundary_nodes, interpolation_points_center, np.array(interpolation_points, dtype=np.int64), line_segments, basis)
    # np.savetxt(file_name, data, fmt="%d %d" + " %.15f" + " %d"*8 + " %.15f"*(1 + 4 + 4 + 4) + " %.15f"*2)
    return data


@njit    
def jit_find_interpolation_points(boundary_nodes, geo_org):
    ret = []
    geo = np.copy(geo_org)
    # Remove the boundary nodes as possible interpolation points
    for ib, jb in boundary_nodes:
        geo[ib, jb] = 0    
    # Find the interpolation points
    for ib, jb in boundary_nodes:
        for di in (-1, 0, 1, 2):
            for dj in (-1, 0, 1, 2):
                sum = 0
                for isq, jsq in ((0, 0), (1, 0), (1, 1), (0, 1)):
                    i = ib - di + isq
                    j = jb - dj + jsq
                    sum += geo[i, j] > 0
                if sum == 4: # Check that all nodes are fluid nodes
                    i = ib - di
                    j = jb - dj
                    box = [[i+x, j+y] for x, y in ((0, 0), (1, 0), (1, 1), (0, 1))]
                    ret.append(box)
    return ret

def find_interpolation_points(boundary_nodes, geo):
    """
    Find tuples of four points that can be used to perform bi-linear
    interpolations.
    Ordering of the interpolation points:
             (0, 1) 3---------2 (1, 1)
                    |         | 
                    |         |
                    |         |
             (0, 0) 0---------1 (1, 0)
                          
    Parameters
    ----------
    boundary_nodes: array
        Array of boundary positions [[i0, j0], [i1, j1], ...]
    geo : array
        2D array where geo = 0 marks solid nodes and 
        geo > 0 marks fluid nodes

    Returns
    -------
    list
        list of points that can be used for interpolations:
          (((x0, y0), (x1, y1), (x2, y2), (x3, y3)), ...)
        where x and y are coordinates and the tuple of four points 
        define a square of neighboring fluid nodes

    Notes
    -----
    The list-set-tuple construction, in the return statement, removes repeated squares. Need 
    use tuples to crate an immutable set of objects that can be fed
    into the set that removes equal entries. NB the entires of ret
    are the  ((x0, y0), (x1, y1), (x2, y2), (x3, y3)) tuples, so all
    node positions must equal if a entry is to be removed. (They also 
    need to be in the same order, which should be fixed)
    """
    ret = jit_find_interpolation_points(boundary_nodes, geo)
    return list(set(tuple(tuple(tuple(y) for y in x) for x in ret)))

@njit    
def jit_find_interpolation_2nd_order_points(boundary_nodes, geo_org):
    ret = []
    geo = np.copy(geo_org)
    # Remove the boundary nodes as possible interpolation points
    for ib, jb in boundary_nodes:
        geo[ib, jb] = 0    
    # Find the interpolation points
    node_positions = (-1, -1), (0, -1), (1, -1), (1, 0), (0, 0), (-1, 0), (-1, 1), (0, 1), (1, 1)
    for ib, jb in boundary_nodes:
        for di in (-2, -1, 0, 1, 2):
            for dj in (-2, -1, 0, 1, 2):
                sum = 0
                for isq, jsq in node_positions:
                    i = ib - di + isq
                    j = jb - dj + jsq
                    sum += geo[i, j] > 0
                if sum == 9: # Check that all nodes are fluid nodes
                    i = ib - di
                    j = jb - dj
                    box = [[i+x, j+y] for x, y in node_positions]
                    ret.append(box)
    return ret

def find_interpolation_2nd_order_points(boundary_nodes, geo):
    """
    Find tuples of nine points that can be used to perform bi-quadratic
    interpolations.
    Ordering of the interpolation points:
                            
           (-1,  1) 6---------7---------8 (1,  1)
                    |      (0,  1)      |
                    |         |         |
                    |         |         |
           (-1,  0) 5---------4---------3 (1,  0)
                    |      (0,  0)      |
                    |         |         |
                    |         |         |
           (-1, -1) 0---------1---------2 (1, -1)
                           (0, -1)
                          
    Parameters
    ----------
    boundary_nodes: array
        Array of boundary positions [[i0, j0], [i1, j1], ...]
    geo : array
        2D array where geo = 0 marks solid nodes and 
        geo > 0 marks fluid nodes

    Returns
    -------
    list
        list of points that can be used for interpolations:
          (((x0, y0), (x1, y1), (x2, y2), (x3, y3)), ...)
        where x and y are coordinates and the tuple of four points 
        define a square of neighboring fluid nodes

    Notes
    -----
    The list-set-tuple construction, in the return statement, removes repeated squares. Need 
    use tuples to crate an immutable set of objects that can be fed
    into the set that removes equal entries. NB the entires of ret
    are the  ((x0, y0), (x1, y1), (x2, y2), ..., (x8, y8)) tuples, so all
    node positions must equal if a entry is to be removed. (They also 
    need to be in the same order)
    """
    ret = jit_find_interpolation_2nd_order_points(boundary_nodes, geo)
    return list(set(tuple(tuple(tuple(y) for y in x) for x in ret)))

@njit
def find_nearest_interpolation_point(pos, boundary_normal, interpolation_points_center, cut_off=3.5):
    """
    Parameters
    ----------
    pos: array
        Position on the form [x, y]
    interpolation_points_center: array
        Center position of the interpolation points from ```find_interpolation_points```.
    cut_off: double (default=3.5)
    
    Returns
    -------
    int:
        index of interpolation points
    """
    # ----------------------------------------------- Mean pos
    #square_centers = np.mean(interpolation_points, axis=1)   
    dr = interpolation_points_center - pos
    dr_sq = np.sum(dr**2, axis=1)
    dr_dot_n = np.dot(dr, boundary_normal)
    val = np.sqrt(dr_sq - dr_dot_n**2)
    val[dr_sq > cut_off**2] = 2*cut_off
    val[dr_dot_n < 0] = 2*cut_off
    ret = np.argmin(val)
    if dr_sq[ret] > cut_off**2:
        print(ret, boundary_normal)

    return ret
    



if __name__ == "__main__":
    II, JJ = np.meshgrid(np.arange(100), np.arange(100), indexing='ij')

    # wing profile
    theta = np.linspace(0, 2*np.pi, 200)[:-1]
    r = lambda x: 0.3*np.sin(3*x) + 0.7
    lin_seg = np.array(
        [[40*r(t)*np.cos(t) + 49.5, 40*r(t)*np.sin(t) + 49.5] for t in theta]
    )

    # # simple box
    # lin_seg = np.array(
    #     [[39.5, 39.5], [59.5, 39.5], [59.5, 59.5], [39.5, 59.5]]
    # )
    # TEST
    #seg = np.array([[1, 1],[2, 2]])

    geo = mask_polygon(II, JJ, lin_seg)
    #geo[:50, :] = 2*geo[:50,:]

    vtk = vlb.vtklb(geo, "D2Q9", "xy", path="/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/input/mpi/")
    vtk.append_data_set("ii", II)
    vtk.append_data_set('jj', JJ)

    basis = vtk.get_basis("D2Q9")
    bnd = markBoundaryNodes(geo, basis)
    destination_folder = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/input/mpi/"

    plt.figure()
    plt.pcolormesh(II, JJ, bnd)



    for rank in range(1, 2):
        ij = np.argwhere((bnd == 1) & (geo == rank))
        interpolation_points = find_interpolation_2nd_order_points(ij, geo)

        data_node = write_boundary_file_node_2nd_order(
            destination_folder + f"boundary{rank-1}.txt", 
            ij, 
            interpolation_points, 
            lin_seg, 
            basis
            )

        for pts in interpolation_points:
            tmp = list(pts) + [pts[0]] 
            x = [p[0] for p in tmp]
            y = [p[1] for p in tmp]
            plt.plot(x, y, '-')
        interpolation_points_mean = np.mean(interpolation_points, axis=1)
        for x, y in interpolation_points_mean:
            plt.plot(x, y, 'cx')
        # for px, py, gamma, p0x, p0y, p1x, p1y, p2x, p2y, p3x, p3y, ip_gamma, *tmp, nx, ny in data_node:
        #     w0, w1, w2, w3 = tmp[8:12]
        #     x = [px, px - gamma*nx]       
        #     y = [py, py - gamma*ny]
        #     plt.plot(x, y, 'k-')       
        #     x = [px, px + ip_gamma*nx]       
        #     y = [py, py + ip_gamma*ny]
        #     plt.plot(x, y, 'c-')
        #     x = w0*p0x + w1*p1x + w2*p2x + w3*p3x
        #     y = w0*p0y + w1*p1y + w2*p2y + w3*p3y
        #     plt.plot(x, y, 'wo')       

    plt.axis("equal")
    plt.show()
