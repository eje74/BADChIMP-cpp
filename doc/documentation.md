# LB Documentation

## Python programs

### Geometry files
The geometry is defined by a multi dimensional array initiated as
```python
	geo = np.zeros([nz, ny, nx])
```
where ```nx```, ```ny```, ```nz``` is the number of nodes the x, y and z direction.  
The ```geo``` matrix only contains ones and zeros  

- __0__= __Fluid__ node
- __1__= __Solid__ node 

The geometry is written to file with the function [```write_geo_file```](#write-geo-file)

### Make the mpi grid files  
First we need to read the file made by the ```write_geo_file``` by using
```python
geo_input =  readGeoFile(file_name):
```
where ```file_name``` is the name of the geometry file.  
```geo_input``` has the same content as the geometry file, but ones are set to zeros and zeros to ones.  

- __0__= __Solid__ node
- __1__= __Fluid__ node  

The dimensions are also reversed, so that ```geo_input``` has the structure ```[nx, ny, nz]```.


### Structure of a one phase fluid solver

#### Preamble
First we will set the lattice type with:
```cpp
// SET THE LATTICE TYPE
#define LT LTYPE
```
where ```LTYPE``` is one of the predefined [lattice types](#lattice-types) (e. g. ```D3Q19```) 

#### main 
**Setup MPI** with
```cpp
    // SETUP MPI
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
```
_nProcs_ is the total number of processes  
_myRank_ is the rank of the current process  

**Setup input directories**  
```cpp
    // SETUP THE INPUT AND OUTPUT PATHS
    std::string chimpDir = "./";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
```
_chimpDir_ is the root directory  
_mpiDir_ is the directory with the files describing the partitioning of the system.  
_inputDir_ is the directory that contains the _input file_  

**Read the files with geometrical information**  
These file contains the information about the geometry, boundary conditions and mpi information
```cpp
    // READ BOUNDARY FILES WITH MPI INFORMATION
    MpiFile<LT> rankFile(mpiDir + "rank_" + std::to_string(myRank) + "_rank.mpi");
    MpiFile<LT> localFile(mpiDir + "rank_" + std::to_string(myRank) + "_local_labels.mpi");
    MpiFile<LT> globalFile(mpiDir + "rank_" + std::to_string(myRank) + "_global_labels.mpi");
```
All these files are [MpiFile](#mpifile) objects.

**Setup the grid and node structures**  
```cpp
    // SETUP GRID
    std::cout << "grid" << std::endl;
    auto grid  = Grid<LT>::makeObject(localFile, rankFile, myRank);
    std::cout << "grid.size = " << grid.size() << std::endl;

    // SETUP NODE
    std::cout << "nodes" << std::endl;
    auto nodes = Nodes<LT>::makeObject(localFile, rankFile, myRank, grid);
```
The [_grid_](#grid) object contains the information about the neighbors of a given node. 

## APPENDIX

### Python function
#### write geo file
```python
	def write_geo_file(filename, geo, res=1.0):
```
where ```filename``` is the output filename (including path) for a geometry file and ```geo```  is the integer multidimensional matrix
### Class descriptions
#### MpiFile 
We should improve this part of the code.


#### Grid 
The Grid calls holds all information about the neighborhood of a given node in the system.  
Grid objects are called using templates  
```cpp
Grid<LTYPE> grid;
```
where ```LTYPE``` is a [lattice type](#lattice-types).  
Functions

* _neighbor_: Overloaded function. Returns the node number of a node in a given direction, _qNo_, or a list of all neighbors for a given node, _nodeNo_.  
```int neighbor(const int qNo, const int nodeNo)```  
Returns the node number of _nodeNo_ neighbor in the _qNo_ direction  
```std::vector<int> neighbor(const int nodeNo)```  
Returns a list of all neighbors for _nodeNo_  
* _pos_: Overloaded function. Returns the Cartesian coordinates of a nodes position.  
```std::vector<int> pos(const int nodeNo)```  
returns a vector containing the Cartesian coordinates of _nodeNo_.  
```int& pos(const int nodeNo, const int index)```  
returns the value of the _nodeNo_'s coordinate with index _index_. 


#### Fields
There are three type of fields in the code. Scalars, Vectors and LB fields:
```cpp
	// Scalar
	ScalarField s_field(num _phases, num_nodes);
	// Vector
	VectorField<LTYPE> v_field(num _phases, num_nodes);
	// LB
	LbField<LTYPE> f_field(num _phases, num_nodes);
```
where ```LTYPE``` is a [lattice type](#lattice-types), ```num_phases``` is the number of different field (e.g. number of phases or components), and ```num_nodes``` is the number of nodes in the system, usually given by ```Grid.size()```.

### Lattice types
```D3Q19```