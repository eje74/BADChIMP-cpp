#### Coding
* Add valarray option for the python script generating the DxQy.h files [23/08 - 19]

#### Data structure
* Add a nodes class. Stores data for:  
  * position
  * node type:
    * bulkSolid
    * bulkNode
    * boundarySolid
    * boundaryNode
    * boundaryMpi?
  * rank ?  
* Use std::vector
* Which structures should we hide?
#### Input / Output
* Format of VTK-files: structured and/or unstructured
* Should VTK be our base method to store data?

   If yes, then we need to have a python-lib that can read them.
* Use BADChIMP read functions for the input.dat files
* What should we do about the  *.mpi* files.

   We should have a method that is not dependent on defining the lattice-basis in the python file.
