#!/usr/bin/env python3

import pyvista as pv
import numpy as np
from scipy.spatial import KDTree
from pathlib import Path
import sys

path = Path(sys.argv[1])
dx = float(sys.argv[2])

# Read surface mesh
mesh = pv.read(path)
mesh_cells = mesh.cell_centers().points.astype(np.double)

# Create uniform grid
xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
size = np.array((xmax-xmin, ymax-ymin, zmax-zmin))
dim = np.ceil(size/dx).astype(int)
grid = pv.UniformGrid()
grid.dimensions = dim + 1 # Cell data
grid.origin = (xmin, ymin, zmin)
grid.spacing = (dx, dx, dx)
grid_cells = grid.cell_centers().points.astype(np.double)
print(dim)

# Shortest distance from grid-cells and mesh surface 
# with index of mesh-cell
dist_idx = np.array(KDTree(mesh_cells).query(grid_cells))

# Inside or outside surface?
idx = dist_idx[1].astype(int) 
norm_dot_vec = mesh.cell_normals[idx] * (mesh_cells[idx]-grid_cells)
inside = np.sum(norm_dot_vec, axis=1) > 0
sign = -np.ones(inside.shape)
sign[inside] = 1

# Save distance from grid node to surface, positive values inside surface 
grid.cell_data['distance'] = dist_idx[0]*sign
map = path.parent/'distance_map.vtk'
grid.save(map)
print(f'  Saved {map}')

#p = pv.Plotter()
# p.add_mesh(mesh, show_edges=True, line_width=1)
# p.add_mesh(mcen, color="r", point_size=8.0, render_points_as_spheres=True)
# p.show()

