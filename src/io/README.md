IO Module Guide
===============

Overview
--------
This folder provides the input parser and the VTK output helpers used by the
simulation code. Most users only interact with Input (for reading config files)
and Output (for writing VTK files).

Files
-----
- Input.h: Hierarchical input parser. Reads tagged blocks and scalar/vector
  values into a tree of Block objects.
- Output.h: Convenience wrapper around VTK::Output for writing fields.
- VTK.h: Low-level VTK writer for VTU/PVTU files, including data marshaling.

Input.h quick start
-------------------
Example file:
  <fluid>
    viscosity 0.001
    <bodyforce>
      0.0 0.0 0.0
    <end>
  <end>

Example usage:
  Input input("input.dat");
  double nu = input["fluid"]["viscosity"];
  std::vector<double> bf = input["fluid"]["bodyforce"];

Notes:
- Use <tag> ... <end> to define nested blocks.
- Values default to double; "int" and "char" can be specified on the tag line.
- Simple math expressions are supported unless disabled via flags.

Output.h quick start
-------------------
Example usage:
  Output<LT> out(grid, nodes, outputDir, rank, nprocs);
  out.add_file("lb_run");
  out.add_scalar_variables({"rho"}, {rho});
  out.add_vector_variables({"vel"}, {vel});
  out.write(step);

Float32 binary output:
  Output<LT, float> out(grid, nodes, outputDir, rank, nprocs);
  // All variables added to this instance are written as Float32.

Notes:
- Output owns a VTK::Output instance configured by template parameters:
  T controls the output precision (float or double) and FMT controls ASCII vs
  binary (default is binary).
- For 2D lattices, VTK::voxel is automatically mapped to VTK::pixel via the
  Cell<voxel,2> specialization in VTK.h.

VTK.h quick start
----------------
VTK::Output is the low-level writer. It supports:
- Multiple output files per instance (add_file).
- Variables as std::vector or std::valarray.
- Automatic handling of contiguous and indexed (non-contiguous) data.

The cast wrappers in VTK::Output allow writing Float32 even when input data is
stored as double.
