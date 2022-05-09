#!/usr/bin/env python3
from vtklb import vtklb
#import numpy as np
#import matplotlib.pyplot as plt
import sys
from pathlib import Path
from pyvista import read as vtkread, UnstructuredGrid

dim = 3
chimpdir = Path.home()/"github/BADChIMP-cpp"
geofile = Path(sys.argv[1])
mesh = vtkread(geofile)

UnstructuredGrid.cell

# generate geometry input file(s)
vtk = vtklb(geo, 'D2Q9', 'x', 'tmp', chimpdir/'input'/'mpi') 

# Set initial density
# - get the Cartesian coordinates
X, Y = np.mgrid[0:40, 0:32]
# - set rho
rho = 1.0 + 0.1*np.sin(2*np.pi*X/40)

# Add rho to the geometry files
vtk.append_data_set("init_rho", rho)

plt.figure(1)
rho[geo==0] = np.nan
ax = plt.pcolormesh(rho.transpose())
plt.axis('equal')
plt.axis('off')
plt.colorbar()
#plt.savefig(path_badchimp + "std_init_rho.pdf")
plt.show()

