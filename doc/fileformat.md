# VTK like Geometry file format


## Overview of the file format (same as vtk)
```
# BADChIMP GeoFile Version 0.1   <header>
one line descriptioin            <title> 
ASCII | BINARY                   <data type>
DATASET type                     <geometry/topology>
...
POINT_DATA n                     <dataset attributes>
...
```

### Geometry/Topology
```
DATASET UNSTRUCTURED_LB_GRID
NUM_DIMENSIONS nd               <nd: int number of spatial dimension>
USE_ZERO_GHOST_NODE             <use node 0 as the default ghost node>
POINTS n dataType               <n: number of points>
p1_xp1_y...                     <point/node spatial position>
p2_xp2_y...
...
LATTICE nq dataType             <nq: number of basis vectors>
c0_xc0_y...                     <basis vecotor elements>
c1_xc1_y...
...
NEIGHBORS dataType              <list of neighbor nodes for each node>
i1_0i1_1...i1_(nq-1)            <node number of neighbor nodes>
i2_0i2_1...i2_(nq-1)
...
in_0in_1...in_(nq-1)
PARALLEL_COMPUTING rank         <rank: rank of the current processor>
PROCESSOR n rank                <n: number of point, rank: neighbor processor>
i0j0                            <i: node this rank, j: node neighbor rank>
i1j1
...
i(n-1)j(n-1)
```
#### Comments to entries

- ```NUM_DIMENSIONS```: (Optional, Default: nd=3)  Specify number of spatial dimensions. In vtk it seems that all vectors and positions are given with three components. 

- ```USE_ZERO_GHOST_NODE```: (Optional) Treat node 0 as a placeholder for a node that is not in use (i.e solid nodes). This means that we begin numbering points from 1, _not_ 0. If the key word is written then this feature is enabled.

- ```POINT```: Same as for the vtk-format. NB if ```USE_ZERO_GHOST_NODE``` is enabled then the the first entry has node number 1, the next has 2 and so on. If ```USE_ZERO_GHOST_NODE``` is disabled  the first entry has node number 0.

- ```LATTICE```: Description of the set of basis vector. We can check these with the ones used in the BADChIMP code and map one set to the other if that is needed.

- ```NEIGHHBORS```:  For each node there is a list of ```nq``` nodes. The node number corresponds to the position in the list under the _points_ keyword and the position in the list of neighbors corresponds to the basic vector in the list under the _lattice_ keyword.

- ```PARALLEL_COMPUTING```: Begins the block describing the processor-processor communication. _rank_ is the rank of the current processor.

- ```PROCESSOR```:  Information about neighboring processor. _n_ is the number of nodes at the current processor that is used to represent the nodes at a neighboring process. _rank_ is the rank of the neighboring node. The list under the key-word the first entry, _i_, is the node number in the geometry in the current processor that represent the node with number _j_ in the neighboring rank, the second entry on the line.

### Dataset attributes
Here we can use the same formalism as the vtk-file format. That is use the key-words:

- ```SCALARS dataName dataType```

- ```VECTORS dataName dataType```

- ```TENSORS dataName dataType```
