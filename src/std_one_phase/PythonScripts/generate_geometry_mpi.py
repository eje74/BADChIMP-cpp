import numpy as np
import os
import subprocess
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

# ======================================================================== Domain decomposition
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

# ======================================================================== Generate geometry
def generate_geometry_mpi(geo, num_proc=(3,)*3, inputpath = r"../../../input/mpi/"):
    """
    Takes the geometry array, removes all non percolating 
    clusters and generates the lb-geometry inputfiles if
    there exist at least on percolating cluster.

    Returns the number of processors needed for the 
    simulations. If the number of processors are zero
    the there were no percolating clusters.

    Paramters
    ---------
    geo : ndarray(int)
        3d matrix defining the geometry
    
    num_proc: list-like, default: (3, 3, 3)
        Number of bins in each spatial direction so that
        `prod(num_proc) is the maximum number of domains

    Return
    ------
    int :
        0   : no percolating cluster (no geometry files are written)
        > 0 : number of processors needed for the simulation.
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
        geo, numproc = domain_decomposition_for_mpi(geo, num_proc)
        # Write the geometry for the lb-run
        vtk = vtklb(geo, "D3Q19", "", path=inputpath)
        # and add the tag
        vtk.append_data_set("geo_tag", geo_tag)
    return numproc

# ======================================================================== Generate geometry
def write_input_file(pathlb,
                     max_iterations,
                     write_interval,
                     tau,
                     bodyforce_z,
                     filebasename):
    """
    Writes the input file used by BADChIMP

    Paramters
    ---------
    pathlb : string
        Path to main/source badchimp folder

    max_iterations: int
        Maximum number of iterations
    
    write_interval: int
        Number of iterations between file write

    tau: float
        LB relaxation time

    bodyforce_z: float
        Body force component in the z-spatial direction

    filebasename: string
        Name used for output files.
    """
    with open(pathlb + r"input/input.dat", "w") as file:
        file.write("# ----------------------------\n")
        file.write("# input for relperm run test\n")
        file.write("# ----------------------------\n")
        file.write("<iterations>\n")
        file.write(f"  max   {max_iterations}         # stop simulation after\n")
        file.write(f"  write {write_interval}           # write interval in steps\n")
        file.write("<end>\n")
        file.write("# ----------------------------\n")
        file.write("# fluid input:\n")
        file.write("# ----------------------------\n")
        file.write("<fluid>\n")
        file.write(f"  tau {tau}\n")
        file.write(f"  bodyforce 0.0 0.0 {bodyforce_z}            # Body force\n")
        file.write("<end>\n")
        file.write("<filenames>\n")
        file.write(fr"  basename {filebasename}" + "\n")
        file.write("<end>\n")
        file.close()   

        
# ======================================================================== New run
def new_run(geo,
            pathlb,
            filebase,
            num_proc=(3, 3, 3),
            key_phrase="ITERATION",
            max_iterations=50000,
            write_interval=500,
            tau=0.8,
            bodyforce_z=1.0e-6
            ):
    """
    Starts a new badchimp run. Takes the geo and generate the lbvtk input files.
    Chekcs if there exists any percolating structures. 
    If yes:
        runs the badchimp code.
        Genreates the filbase.dat file contianing the sum of velocities, and the
        filebase.vtk file for paraview plotting
    else:
        Just generates the filbase.dat file with zero velocity
      
    Parameters
    ----------
    geo : ndarray(int)
        3d matrix defining the geometry
    
    pathlb : string
        Path to main/source badchimp folder

    filebase: string
        Name used for output files namse.

    num_proc: list-like, default: (3, 3, 3)
        Number of bins in each spatial direction so that
        `prod(num_proc) is the maximum number of domains

    max_iterations: int
        Maximum number of iterations
    
    write_interval: int
        Number of iterations between file write

    tau: float
        LB relaxation time

    bodyforce_z: float
        Body force component in the z-spatial direction

    Return
    ------
    int:
        0   : no percolating cluster (no geometry files are written)
        > 0 : number of processors needed for the simulation.

    None/subprocess-object:
        None if no percolating cluster else subprocess-object of the 
        badchimp run. 
    """
    # Defualt value of proc is set to None. Will change if nproc > 0
    proc = None

    # Removes non percolating clusters and partition the system for
    # parallell runs
    nproc = generate_geometry_mpi(geo, num_proc, pathlb + r"input/mpi/")

    if nproc > 0:
        # write the input file
        write_input_file(pathlb,
                        max_iterations,
                        write_interval,
                        tau,
                        bodyforce_z,
                        filebase)
        # mpirun command with correct number of processes
        command = ["nohup", r"mpirun", r"-np", str(nproc), r"bin/bdchmp", "&"]
        # set the working dir to build so that we do not need to
        # changed the paths in the badchimp code
        workingdir = pathlb + r"build/"
        proc = subprocess.Popen(" ".join(command), 
                                # We must use the join-command to accomodate 
                                # the shell=True choice
                                cwd=workingdir, 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE, 
                                text=True,
                                shell=True)
                                # We use shell=True so that BADChIMP can be 
                                # run in the background (ie. nohup and &)
        # When badchimp outputs "ITERATION" we know that it has 
        # read the geometry files
        for line in proc.stdout:
            if key_phrase in line:
                break
        #... and we are ready to begin another run
    else:
        # Write zero flux to file
        with open(pathlb + "output/" + filebase + ".dat", "w") as file:
            file.write("0, 0.0, 0.0, 0.0\n")
            file.close()

    return nproc, proc
