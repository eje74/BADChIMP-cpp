import numpy as np
import os
import subprocess
import time
from generate_geometry_mpi import new_run

# ======================================================================== Automate runs
running_processes = []
key_phrase = "ITERATION"
max_lim = 2
# ======================================================================== Path to BADChIMP
pathlb = os.path.expanduser(r"~/GitHub/BADChIMP-cpp/")

# ======================================================================== Data format
# Solid geometry:
#   element > 0: solid
# Fluid geometry:
#   element > 0: wetting

# ======================================================================== Path to data
# ------------------------------------------------------------------------ Clone onedrive
os.system("rclone-mount.sh")
time.sleep(2)
# ------------------------------------------------------------------------ root path 
pathinput = os.path.expanduser(r"~/OneDrive/NORCE/CSSR/RelPerm LB LS/") 
# ------------------------------------------------------------------------ system
pathinput += r"Castlegate_Tow20_LVC_Oil/Sat_control/PrimaryDrainage/"
# ------------------------------------------------------------------------ fluid phases
filenames = [x[:-4] for x in os.listdir(pathinput) if r"CG_NWP" in x]

# ======================================================================== Rock data
filebase = r"CG_PoreSolid_200x200x200_SDF_PD"
# ------------------------------------------------------------------------ load data
pore = np.load(pathinput + filebase + r".npy")
# ------------------------------------------------------------------------ generate geometry
geo = np.ones(pore.shape, dtype=np.int32)
geo[pore>0] = 0
# ------------------------------------------------------------------------ new proc rock
nproc, proc = new_run(geo, pathlb, filebase)
if proc:
    running_processes.append(proc)

for filebase in filenames:
  while len(running_processes) == max_lim:
    # Take some time between checks
    time.sleep(10)
    for proc in running_processes:
      # If process has terminated
      if proc.poll() is not None:
        # Propely halt the process
        proc.communicate()
        # Remove it from the list of
        # running processes
        running_processes.remove(proc)
        # Break the for-loop
        break
        # We know that the number of running processes
        # are less then max_lim so we will also exit
        # the while-loop
  # We add a new process to the running processes
  fluid = np.load(pathinput + filebase +r".npy")
  geo[:] = 1
  # Wetting phase
  geo[fluid<=0] = 0
  geo[pore>0] = 0
  nproc, proc = new_run(geo, pathlb, filebase + r"_W")
  if proc:
    running_processes.append(proc)
  while len(running_processes) == max_lim:
    time.sleep(10)
    for proc in running_processes:
      if proc.poll() is not None:
        proc.communicate()
        running_processes.remove(proc)
        break
  # Non-wetting phase
  geo[:] = 1
  geo[fluid>0] = 0
  geo[pore>0] = 0
  nproc, proc = new_run(geo, pathlb, filebase + r"_NW")
  if proc:
    running_processes.append(proc)
# End the rest of the processes
while len(running_processes) > 0:
    for proc in running_processes:
        if proc.poll() is not None:
            proc.communicate()
            running_processes.remove(proc)

#nproc = generate_geometry_mpi(geo, (2,)*3, pathlb + r"input/mpi/")


# nproc = generate_geometry_mpi(geo, (2, 1, 1), pathlb + r"input/mpi/")

# if nproc > 0:
#     # write the input file
#     write_input_file(pathlb,
#                      5,
#                      1,
#                      0.8,
#                      1.0e-6,
#                      filebase)
#     # mpirun command with correct number of processes
#     command = ["nohup", r"mpirun", r"-np", str(nproc), r"bin/bdchmp", "&"]
#     # set the working dir to build so that we do not need to
#     # changed the paths in the badchimp code
#     workingdir = pathlb + r"build/"
#     proc = subprocess.Popen(" ".join(command), 
#                             # We must use the join-command to accomodate 
#                             # the shell=True choice
#                             cwd=workingdir, 
#                             stdout=subprocess.PIPE, 
#                             stderr=subprocess.PIPE, 
#                             text=True,
#                             shell=True)
#                             # We use shell=True to be able to run
#                             # BADChIMP in the background (ie using
#                             # nohup and &)
#     # When badchimp outputs "ITERATION" we know that it has 
#     # read the geometry files
#     for line in proc.stdout:
#         if key_phrase in line:
#             break
#     #... and we are ready to begin another run
# else:
#     # Write zero flux to file
#     with open(pathlb + "output/" + filebase + ".dat", "w") as file:
#         file.write("0, 0.0, 0.0, 0.0\n")
#         file.close()
    
# print("number of procs = ", nproc, flush=True)

# # ======================================================================== Fluid phase data
# # ------------------------------------------------------------------------ load data
# fluid = np.load(pathinput + r"CG_NWP_Sw0746_200x200x200_SDF_PD.npy")
# # ------------------------------------------------------------------------ wetting generate geometry
# wphase = np.ones_like(geo)
# wphase[fluid<=0] = 0
# wphase[geo==0] = 0
# # ------------------------------------------------------------------------ non-wetting generate geometry
# nwphase = np.ones_like(geo)
# nwphase[fluid>0] = 0
# nwphase[geo==0] = 0



# if False:

#     # ------------------------------------------------------------------------ Plot rock
#     grid = pv.ImageData()
#     grid.dimensions = geo_tag.shape
#     grid.origin = (0, 0, 0)
#     grid.spacing = (1, 1, 1)
#     grid.point_data["geo"] = geo_tag.flatten(order="F")
#     #grid.plot(show_edges=False, opacity="linear")
#     grid.plot(volume=True, opacity=[0.0, 0.9], shade=False)

# if False:
#     grid = pv.ImageData()
#     grid.dimensions = geo.shape
#     grid.origin = (0, 0, 0)
#     grid.spacing = (1, 1, 1)
#     grid.point_data["geo"] = nwphase.flatten(order="F")
#     #grid.plot(show_edges=False, opacity="linear")
#     grid.plot(volume=True, opacity=[0, 0.5], shade=False)

#     grid = pv.ImageData()
#     grid.dimensions = geo.shape
#     grid.origin = (0, 0, 0)
#     grid.spacing = (1, 1, 1)
#     grid.point_data["geo"] = wphase.flatten(order="F")
#     #grid.plot(show_edges=False, opacity="linear")
#     grid.plot(volume=True, opacity=[0, 0.5], shade=False)
