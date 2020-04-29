# LB Documentation

## Recipe: Add a geometry

1. Use the python script ```./PythonScripts/writeGeoFile.py``` as a guid for writing a geo file that is to be read by ```./PythonScripts/mpiGrid.py```
2. In  ```mpiGrid.py``` set  
   ```write_dir``` to the directory where the geo-file was written  
   ```file_name``` is the name of the geo-file  
3. Run the  ```mpiGrid.py``` script
4. The .mpi files are written to ```write_dir```. Move these files to the directory read by _BADChIMP_ which is given in ```mpiDir``` in the main-file.

  

## Python programs

### Geometry files
The geometry is defined by a multi dimensional array initiated as
```python
	geo = np.zeros([nz, ny, nx])
```
where ```nx```, ```ny```, ```nz``` is the number of nodes the x, y and z direction.  
The ```geo``` matrix only contains non negative integer values  

- __0__= __Solid__ node
- __0__<  __Fluid__ node 


The integer value of a fluid nodes give the processor rank it belongs to.
The geometry is written to file with the function [```write_geo_file```](#write-geo-file)

### Make the mpi grid files  
First we need to read the file made by the ```write_geo_file``` by using
```python
geo_input =  readGeoFile(file_name):
```
where ```file_name``` is the name of the geometry file.  
```geo_input``` has the same content as the geometry file.  

- __0__= __Solid__ node
- __0__ < __Fluid__ node  


The _domain decomposition_  is given in ```geo_input``` so that the fluid node _n_ belongs to the process with rank ```geo_input[n]-1```.

Next we need to find the size of the rim that we need to include, that is, the size of the neighborhood given by ```c```,
```python
rim_width = getRimWidth(c)
```
and make a new geometry matrix, ```geo```, including the rim
```python
# -- add rim
geo = addRim(geo_input, rim_width)
``` 
Here it is assumed that the rim values are periodic (marked with -1). Other rim values can be given following the line of comments 

```python
# SETUP RIM VALUES
## Here we need to set aditional values for the rim
# set as periodic -1 (default value)
# set as solid = -2
# set as fluid = -<#rank> - 3
```

Afther the comment 
```python
# -- END user input
```
there should be no more need for user input. The rest of the script is described in the [appendix](#mpigrid)

#### Files generated  

The files are written into the folder in ```write_dir```, as an example
```python
write_dir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/PythonScripts/"
```
For each processor we generate three files with ```.mpi```-extensions. The file names for processor 0 are ```rank_0_global_labels.mpi```, ```rank_0_local_labels.mpi ``` and ```rank_0_rank.mpi```.  
These files should be moved to the directory read by _BADChIMP_ which is given in ```mpiDir``` (see box below)
```cpp
    // SETUP THE INPUT AND OUTPUT PATHS
    std::string chimpDir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/";
    std::string mpiDir = chimpDir + "input/mpi/";
```



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
where ```filename``` is the output filename (including path) for a geometry file and ```geo```  is an integer multidimensional matrix, with ```geo.shape = (nz, ny, nx)```
The file has the structure
```
dim nx ny nz
res <double>
<string>
geo(0,0,0)geo(0,0,1)...geo(0,0,nx-1)
geo(0,1,0)geo(0,1,1)...geo(0,1,nx-1)
...
geo(0,ny-1,0)geo(0,ny-1,1)...geo(0,ny-1,nx-1)
geo(1,0,0)geo(1,0,1)...geo(1,0,nx-1)

```

### Python scripts
#### mpiGrid
Firstly the geometry file is read by 
```python
geo_input = readGeoFile(file_name) 
```
that has the shape ```geo_input.shape = (nz,ny,nx)``` and 0 is a _solid_ node and a non zero value is a _fluid_ node.  
Example of a system that is partitioned between two processors with rank 0 and 1:

$$
\begin{matrix}0 & 0 & 0 & 0 & 0\\2 & 2 & 2 & 2 & 2\\2 & 2 & 2 & 2 & 2\\2 & 2 & 2 & 2 & 2\\1 & 1 & 1 & 1 & 1\\1 & 1 & 1 & 1 & 1\\0 & 0 & 0 & 0 & 0\end{matrix}
$$
and then we added a default rim 
```python
# -- add rim
geo = addRim(geo_input, rim_width)
```
$$
\begin{matrix}-1 & -1 & -1 & -1 & -1 & -1 & -1\\-1 & 0 & 0 & 0 & 0 & 0 & -1\\-1 & 2 & 2 & 2 & 2 & 2 & -1\\-1 & 2 & 2 & 2 & 2 & 2 & -1\\-1 & 2 & 2 & 2 & 2 & 2 & -1\\-1 & 1 & 1 & 1 & 1 & 1 & -1\\-1 & 1 & 1 & 1 & 1 & 1 & -1\\-1 & 0 & 0 & 0 & 0 & 0 & -1\\-1 & -1 & -1 & -1 & -1 & -1 & -1\end{matrix}
$$
Other alternatives for boundary nodes are
```python
# set as periodic -1 (default value)
# set as solid = -2
# set as fluid = -<#rank> - 3
```


Then we make the ```node_labels``` array
```python
# Make node labels
node_labels = setNodeLabels(geo, num_proc)
node_labels[ind_solid] = 0
```
$$
\begin{matrix}0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 1 & 2 & 3 & 4 & 5 & 0\\0 & 6 & 7 & 8 & 9 & 10 & 0\\0 & 11 & 12 & 13 & 14 & 15 & 0\\0 & 1 & 2 & 3 & 4 & 5 & 0\\0 & 6 & 7 & 8 & 9 & 10 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\end{matrix}
$$

The periodic boundaries are added for ```geo``` with 
```python
geo = addPeriodicBoundary(ind_periodic, geo, rim_width)
```
changing the ```geo``` matrix to
$$ 
\begin{matrix}0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\\2 & 2 & 2 & 2 & 2 & 2 & 2\\2 & 2 & 2 & 2 & 2 & 2 & 2\\2 & 2 & 2 & 2 & 2 & 2 & 2\\1 & 1 & 1 & 1 & 1 & 1 & 1\\1 & 1 & 1 & 1 & 1 & 1 & 1\\0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\end{matrix}
$$
and similarly we add the periodic boundary conditions for the for the ```node_labels```
$$
\begin{matrix}0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\\5 & 1 & 2 & 3 & 4 & 5 & 1\\10 & 6 & 7 & 8 & 9 & 10 & 6\\15 & 11 & 12 & 13 & 14 & 15 & 11\\5 & 1 & 2 & 3 & 4 & 5 & 1\\10 & 6 & 7 & 8 & 9 & 10 & 6\\0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\end{matrix}
$$

Now we write the processes specific labeling for process ```my_rank```
```python
node_labels_local = setNodeLabelsLocal(geo, node_labels, ind_periodic, my_rank, c, rim_width)
```
For ```my_rank = 2``` it is
$$
\begin{matrix}20 & 16 & 17 & 18 & 19 & 20 & 16\\0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\\15 & 11 & 12 & 13 & 14 & 15 & 11\\5 & 1 & 2 & 3 & 4 & 5 & 1\\10 & 6 & 7 & 8 & 9 & 10 & 6\\20 & 16 & 17 & 18 & 19 & 20 & 16\\0 & 0 & 0 & 0 & 0 & 0 & 0\end{matrix}
$$
and for ```my_rank = 2``` it is
$$
\begin{matrix}0 & 0 & 0 & 0 & 0 & 0 & 0\\20 & 16 & 17 & 18 & 19 & 20 & 16\\5 & 1 & 2 & 3 & 4 & 5 & 1\\10 & 6 & 7 & 8 & 9 & 10 & 6\\15 & 11 & 12 & 13 & 14 & 15 & 11\\25 & 21 & 22 & 23 & 24 & 25 & 21\\0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0\\20 & 16 & 17 & 18 & 19 & 20 & 16\end{matrix}
$$
We note that the solid boundary nodes also are given a node label.

### Mpi files
#### _rank.mpi


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