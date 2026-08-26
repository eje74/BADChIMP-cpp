#### GITHUB
* Add a standard file structure folder 
#### Coding
* Add valarray option for the python script generating the DxQy.h files [23/08 - 19]
* the lattice function should take std::vectors as input. ie cDot et c.

#### Data structure
* Add global dimensions to the grid class [5/7-2022]
* Add a nodes class. Stores data for:  
  * [OK]position -> added to the grid class.
  * [OK]node type:
    * bulkSolid
    * bulkNode
    * boundarySolid
    * boundaryNode
    * boundaryMpi?
  * rank ?  
* [ok]Use std::vector
* [ok]BndMpi.communciateScalarField() fungerer ikke med mer enn ett skalarfelt
* Bør ha muligheten å kommunisere vektor- og tensor-felt.
* Add regular grids
* [ok]apply boundary for multiple fields
##### MPI
* Add neighborhood lists extraction to python scripts, so that just neighborlists are sent to BADChIMP
* write mpi comunicated for multiple fields.
* write mpi comunicated for vector fields

#### Input / Output
* Format of VTK-files: structured and/or unstructured
* Should VTK be our base method to store data?

   If yes, then we need to have a python-lib that can read them.
* Use BADChIMP read functions for the input.dat files
* What should we do about the  *.mpi* files. [Made new geometry file format]

   We should have a method that is not dependent on defining the lattice-basis in the python file.

#### Capabilities
* Take gradients of multispecies scalar fields
* Take div of gradient fields 
$\nabla^2$

---

## Run script for rel.perm calculations  

Which scripts were we using for the old setup (eg. jupiter4): 
* Switch to `relperm_ls`-branch   
* Path to script (relative to local repo path): `/src/std_one_phase/PythonScripts`  
* `lbrelperm.py` Script used to run multiple simulations and to setup geometry  
* `lbrelperm.py` depends on `generate_geometry_mpi.py`
* `generate_geometry_mpi.py` depends on `vtklb.py` and `clean_geometry.py`
* dependency: `numba`, `numpy`, (`itertools`)

### Basic run workflow 
Definitions: 
* pathlb : Path to the root folder for the BADChIMP code  
  eg `r"~/Programs/GitHub/BADChIMP-cpp/"`  
* pathinput : Path to the root of the LS data. The rest of the paths will be specified in `run`.  
  eg `r"~/OneDrive/NORCE/CSSR/RelPerm LB LS/"`  

Script logic:  
* load rock data (file name contains something like `PoreSolid`)  
  `pore = np.load(filenamebase + r".npy")`  
* make a geometry file consisting of only 1's and 0's (solid and void):  
  `geo = np.ones(pore.shape, dtype=np.int32)`  
  `geo[pore>0] = 0`
