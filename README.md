# BADChIMP-cpp
## Build and compile 
**Linux:** Run `./make.sh <name_of_folder_with_main_file>`, in root folder, `std_case` is built if no argument to make.sh is given.

**Windows:** Make sure that open [MPI is installed](https://docs.microsoft.com/en-us/archive/blogs/windowshpc/how-to-compile-and-run-a-simple-ms-mpi-program). Download and run `msmpisetup.exe` and `msmpisdk.msi`.  Install [cmake for Windows](https://cmake.org/). Run cmake from root directory to generate Visual Studio C++ project. Or simply use VSCode.

## A c++ port of the C BADChIMP code

### *Deadlock* prevention in BADChIMP
In the LB code each process has a list of neighboring processes which it
sends and receives data from. The list has the following properties:
1) if process A is in process B's list of neighbors, then process B will be in
process A's list of neighbors.
2) Process A is not part of its own list of neighbors.
3) A list of neighbors is sorted in ascending order, based on rank.  

The send/receive code follows the structure given below:
```cpp
for (auto neigRank: listOfNeighbors) {
  if (myRank < neigRank) {
    // SEND DATA to neigRank
    // RECEIVE DATA from neigRank
  } else {
    // RECEIVE DATA from neigRank
    // SEND DATA to neigRank
  }
}
```
We argue, based  on the
[Coffman conditions](https://en.wikipedia.org/wiki/Deadlock),
that this structure is enough to avoid *deadlock*. Here, we will show that the
assumption of a *circular wait condition* will lead to a contraction.  
We assume that there are *n* processes, *p*, waiting for each other, so that
*p<sub>0</sub>* is waiting for *p<sub>1</sub>*, *p<sub>1</sub>* is waiting for *p<sub>2</sub>*, and so on until *p<sub>n-1</sub>*
is waiting for *p<sub>0</sub>*. We call the set of processes that are part of the circular
wait loop for a circular wait-set (CWS).  
Let *p<sub>k</sub>* be the process with the lowest rank in a CWS. Since it has the lowest
rank, it must be waiting to send to a process, *p<sub>k+1</sub>*, with a higher rank, as
given by the code structure above, and since *p<sub>k+1</sub>* is in *p<sub>k</sub>*'s list of
neighbors, then *p<sub>k</sub>* is in *p<sub>k+1</sub>*'s list of neighbors.  
*p<sub>k+1</sub>* must either be waiting to send to, or receive from, *p<sub>k+2</sub>*. If
*p<sub>k+2</sub>* is equal to *p<sub>k</sub>* we know that *p<sub>k</sub>*'s rank is lower than *p<sub>k+1</sub>*'s
and *p<sub>k+1</sub>* should be waiting to receive data from *p<sub>k</sub>*, but this makes the
*circular wait condition* stated above invalid. Hence, *p<sub>k+2</sub>* must be a
different process that *p<sub>k</sub>*. We then know, by assumption, that the rank of
*p<sub>k+2</sub>* is higher than *p<sub>k</sub>*'s, so that it is in a later position in *p<sub>k+1</sub>*'s
list of neighbors (since it was sorted in ascending order). Hence, *p<sub>k+1</sub>*
should then already have been waiting to receive data from *p<sub>k</sub>*, as the list of
neighbors is traversed from lowest to highest values. But since *p<sub>k</sub>* is already
waiting to send to *p<sub>k+1</sub>*, this will again make the *circular wait condition*
false, as *p<sub>k</sub>* would no longer wait to send to *p<sub>k+1</sub>*.  
Thus, the assumption of a *circular wait* condition leads to a contraction,
which proves, by *reductio ad absurdum*, that the proposed parallel
communication protocol cannot lead to a *deadlock* situation.


### Pseudo code

```
//   input: - konfigurasjons data. (tau ++)
//          - fil med geo.
//          - (evnt. fil med makroskopiske variable)
//
//   lattice boltzmann algorithm:
//         - INIT:
//             f_\alpha(x_i, t) = f^eq(rho, vel) + f^neq(rho, vel, grad\rho, grad\vel, tau , +++)
//
//         Vi antar at effektiviteten til loopene miker med økende type nummer. Slik at vi vil
//         bruke type 1 hvis det er mulig
//         - LOOP TYPE 1:
//            for all bulk nodes do {
//            * Macrosopics
//                 rho(x, t) = func(momenter av f-feltet inngår, ++)
//            * Collision:
//                f*_\alpha(x, t) = f_\alpha(x, t) + \Omega(f_\alpha, tau, rho, vel)
//            * Propagation:
//                f_\alpha(x + c_\alpha, t + 1) = f*_\alpha(x, t)
//            } end all bulk nodes
//
//            for all boundary nodes of type 1 do {
//                * Macroscopic
//                * 'Collision'
//                * Propagation
//            }
//
//
//            for all boundary nodes of type 2 do {
//            }
//
//         - LOOP TYPE 2: (Used for color-gradients)
//            for all bulk nodes do {
//            * Macrosopics
//                 rho(x, t) = func(momenter av f-feltet inngår, ++)
//            } end for all bulk node Macroscopic
//
//            for all boundary nodes of type 1 do {
//                * Macroscopic
//            } end for all boundary node Macroscopic
//
//            for all bulk nodes do {
//            * Collision:
//                f*_\alpha(x, t) = f_\alpha(x, t) + \Omega(f_\alpha, tau, rho, vel)
//            * Propagation:
//                f_\alpha(x + c_\alpha, t + 1) = f*_\alpha(x, t)
//            } end all bulk nodes
//
//            for all boundary nodes of type 1 do {
//                * 'Collision'
//                * Propagation
//            }
//
//
//         - LOOP TYPE 3:
//            for all bulk nodes do {
//            * Macrosopics
//                 rho(x, t) = func(momenter av f-feltet inngår, ++)
//            } end for all bulk node Macroscopic
//
//            for all boundary nodes of type 1 do {
//                * Macroscopic
//            } end for all boundary node Macroscopic
//
//            for all bulk nodes do {
//            * Collision:
//                f*_\alpha(x, t) = f_\alpha(x, t) + \Omega(f_\alpha, tau, rho, vel)
//            } end for all bulk node Collision
//
//            for all boundary nodes of type 1 do {
//                * Collision
//            } end for all boundary node Collision
//
//            * Propagation:
//            for all bulk nodes do {
//                f_\alpha(x + c_\alpha, t + 1) = f*_\alpha(x, t)
//            } end all bulk nodes
//
//            for all boundary nodes of type 1 do {
//                * Propagation
//            } end all boundary nodes
```

### Data structure
#### Bulk Nodes
Collision based:  
f[x][&alpha;] = `f0(x1)` `f1(x1)` `f2(x1)`...`fq-1(x1)`...  
Propagation based:  
f[&alpha;][x] = `f0(x1)` `f0(x2)`...`f0(xN)` `f1(x1)` `f1(x2)`...  
All fields node based:  
`f0(x1)` `f1(x1)` `f2(x1)`...`fq-1(x1)` `rho(x1)` `ux(x1)` `uy(x1)`...

```
// GRID STRUCTURE
// element 0                       to nBulk                                     [ bulk noder              ]
// element nBulk                   to nBulk + nBoundaryNodes1                   [ Boundary nodes of type 1]
// element nBulk + nBoundaryNodes1 to nBulk + nBoundaryNodes1 + nBoundaryNodes2 [ Boundary nodes of type 1]
//  or
//
//  listBulkNodes[] = [(x1,y1,z1), (x2,y2,z2) ,,,]
//  listBoundaryNodes1[] = [...]
// Last is simpler to update dynamically
```


#### Boundary nodes
```
// Different linke types
// beta  - betaHat  : one link crosses boundary. beta: unknow, and betaHat known.
// gamma - gammaHat : no links cross boundary. both known
// delta - deltaHat : two links cross boundary. both unknown
// example ||1,6(beta)     |3(gamma)     |4(delta)     ||
//         || 5,2(betaHat) | 7(gammaHat) | 0 (deltaHat)||
// linkList_ = [1, 6, 3, 4]
// nTypes_   = [2, 1, 1]
//
// Boundary. listOfAllWallBoundaryNodes = [2, 6, 8, 100]
// Boundary. linkList [1,4,5,6,   2,4,7,1,   2,4,7,1,   2,4,7,1,]
//
```  
**General**  
In general we will have *nDirections* / 2 pair of links if we have no rest particle and  
(*nDirections* - 1) / 2 if we have a rest particle.
A pair of directions is a set of two directions where the basis velocities are pointing in opposite directions (in a sense, this can include the rest particle as it is its own reverse directions).  

**Periodic nodes**  
Why copy information from a ghost node to a ghost node.  
Why not copy them directly from a ghost node to a periodic node?  
Need then to know how the boundary nodes are stored in the boundary object.  
We will assume that the periodic nodes will pull their values.  
Assume that we have a periodic node at $(x, y)$ then the links "intersecting" the periodic boundary will either be a $\beta$ or $\gamma$ type of link. Then we know that the outgoing distribution should have gone to eg. $(x + c_{\beta x}, y + c_{\beta y})$ this is then the boundary node $(x', y') = (x + c_{\beta x} \mod nX, y + c_{\beta y} \mod nY)$ (*nX* and *nY* are the size of the system without boundary nodes), so that we should the copy the value to $(x' - c_{\beta x}, y -  c_{\beta y})$.  $x, x', y$, and $y'$ are positions relative to the real geometry (that is, without ghost nodes).

#### Project structure
Use [Hiltmon](https://hiltmon.com/blog/2013/07/03/a-simple-c-plus-plus-project-structure/) suggestion for project structure.

