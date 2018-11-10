# BADChIMP-cpp
## A c++ port of the C BADChIMP code

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

**Multiple makefiles**  
Here, I will need to use multiple *makefiles*. `make` can be run with the *-f* option to supply a makefile name.  
`make -f <my_makefile>`.
