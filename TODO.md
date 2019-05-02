#### Data structure
* Add a node class. Stores data for:  
  * position
  * node type:
    * bulkSolid
    * bulkNode
    * boundarySolid
    * boundaryNode
    * ...
  * rank ?  
* Use std::vector
* What structures should we hide?
#### Input / Output
* Format of VTK-files: structured and unstructured
* Should VTK be our base method to store data?

   If yes, then we need to have a python-lib that can read them.
* Use BADChIMPS read functions for the input.dat files
* What should we do about the  *.mpi* files.

   We should have a method that is not dependent on defining the lattice-basis in the python file. 
