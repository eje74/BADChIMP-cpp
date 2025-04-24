import numpy as np
import os
from vtklb import vtklb
from clean_geometry import bwlabel3D_set


# ======================================================================== Remove non-perculating domains
def remove_non_perculating_domains(geo, marker=(1,)):
    """
    Remove the non percolating clusters, in the z-direction,
    from the geometry.

    Paramters
    ---------
    geo : nd.array(int)
        3d matrix defining the geometry
    
    marker: tuple
        List of lables in `geo` that make up the domains

    Return
    ------
    numpy.array(int)
        Cleand geometry fluid=1, solid=0 
    """
    # ------------------------------------------------------------------------ neighborhood
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
    # ======================================================================== Domains in geo
    domains, max_label = bwlabel3D_set(cc, geo, marker)
    # ========================================================================== Remove non-percolating pore space
    percolating_domains = set(domains[:, :, 0].flatten()) & set(domains[:, :, -1].flatten())
    non_percolating_domains = set(range(1,max_label+1)) - percolating_domains
    for n in non_percolating_domains:
        domains[domains==n] = 0
    domains[domains>0] = 1
    return domains

# ======================================================================== Remove non-perculating domains
def domain_decomposition_for_mpi(geo_mpi, num_proc=(2,)*3):
    """
    Domain decomposition of geometries using a cubical grid
    defined by `num_proc` that gives the linear number of 
    bins for each direction.

    Paramaters
    ----------
    geo_mpi: np.ndarray[int]
        Binary geometry images where it is assumed that
        0=solid and 1=fluid

    num_proc: list-like, default: (2, 2, 2)
        Number of bins in each spatial direction so that
        `prod(num_proc) is the maximum number of domains

    Return
    ------
    np.ndarray[int]:
        Domain decomposed geometry where the domains are
        enuerated from 1 up to the total number of domains

    int:
        Total number of domains

    """
    nx, ny, nz = geo_mpi.shape
    xbe = tuple((x[0], x[-1]+1) for x in np.array_split(np.arange(nx), num_proc[0]))
    ybe = tuple((x[0], x[-1]+1) for x in np.array_split(np.arange(ny), num_proc[1]))
    zbe = tuple((x[0], x[-1]+1) for x in np.array_split(np.arange(nz), num_proc[2]))

    nproc = 0
    for x0, x1 in xbe:
        for y0, y1 in ybe:
            for z0, z1 in zbe:
                if np.any(geo_mpi[x0:x1, y0:y1, z0:z1]>0):
                    nproc += 1
                    geo_mpi[x0:x1, y0:y1, z0:z1] *= nproc
    return geo_mpi, nproc

def generate_geometry_mpi(geo, num_proc=(3,)*3):
    """
    Remove the non percolating clusters, in the z-direction,
    from the geometry.

    Paramters
    ---------
    geo : nd.array(int)
        3d matrix defining the geometry
    
    num_proc: list-like, default: (3, 3, 3)
        Number of bins in each spatial direction so that
        `prod(num_proc) is the maximum number of domains

    Return
    ------
    int :
        0   : no percolating cluster (no geometry files are written)
        > 0 : number of 
    """
    # Default return (No percolating structure)
    numproc = 0
    # Remove the non perculating clusters
    geo_org = remove_non_perculating_domains(geo)
    # Check if there are any percolating structures
    if np.any(geo_org>0):
        # Add a rim to treat the lb-boundary conditions
        new_shape = tuple(x+2 for x in geo_org.shape)
        geo = np.zeros(new_shape, dtype=np.int64)
        geo[1:-1, 1:-1, 1:-1] = geo_org
        del geo_org 
        # Tags for inlet and outlet
        geo_tag = np.zeros_like(geo)
        geo_tag[1:-1, 1:-1, 0] = 1
        geo_tag[1:-1, 1:-1, -1] = 1 
        # Decomposition of the geometry for parallell runs
        geo, numproc = domain_decomposition_for_mpi(geo)
        # Write the geometry for the lb-run
        vtk = vtklb(geo, "D3Q19", "", path="../../../input/mpi/")
        # and add the tag
        vtk.append_data_set("geo_tag", geo_tag)
    return numproc

