# LB Documentation
### Structure of a one phase fluid solver
#### Preamble
First we will set the lattice type with:
```cpp
// SET THE LATTICE TYPE
#define LT <LTYPE>
```
where ```<LTYPE>``` is one of the predefined [lattice types](#lattice-types) (eq. ```D3Q19```) 

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


### Lattice types
```D3Q19```