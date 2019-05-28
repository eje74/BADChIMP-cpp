#ifndef LBGEOMETRY_H
#define LBGEOMETRY_H

#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBgrid.h"
#include "LBnodes.h"
#include "LBboundary.h"
#include "LBhalfwaybb.h"
#include "LBbulk.h"

template<typename DXQY>
std::vector<int> findBulkNodes(const int &myRank, const Grid<DXQY> &grid)
// makeBulkNodes : make a list of bulk node labels. Here we assume that all
//  fluid nodes are also bulk nodes.
{
    std::vector<int> bulkNodes;
    for (int n = 1; n < grid.size(); ++n)
        if (grid.getRank(n) == myRank)
            bulkNodes.push_back(n);

    return bulkNodes;
}


template<typename DXQY>
std::vector<int> findSolidBndNodes(const int &myRank, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
{
    std::vector<int> ret; // List of node numbers of all solid boundary nodes for myRank process
    for (int n = 1; n < grid.size(); n++) { // Loop over all grid nodes excpet the default node (node number = 0)
        if (nodes.getType(n) == 0) { // The node is a solid node
            bool hasFluidNeig = false;
            // Check if the node has a solid neighbor
            for (int q = 0; q < DXQY::nQNonZero_; ++q) {
                int neigNode = grid.neighbor(q, n);
                if (nodes.getRank(neigNode) == myRank) // fluid node
                    hasFluidNeig = true;
            }
            if (hasFluidNeig)  ret.push_back(n);
        }
    }
    return ret;
}




template<typename DXQY>
std::vector<int> findFluidBndNodes(const int &myRank, const Grid<DXQY> &grid)
{
    std::vector<int> ret; // List of node numbers to all fluid boundary nodes for myRank process
    for (int n = 1; n < grid.size(); n++) { // Loop over all grid nodes excpet the default node (node number = 0)
        if (grid.getRank(n) == myRank) { // The node is a fluid node
            bool hasSolidNeig = false;
            // Check if the node has a solid neighbor
            for (int q = 0; q < DXQY::nQNonZero_; ++q) {
                int neigNode = grid.neighbor(q, n);
                if (grid.getType(neigNode) == 0) // Solid Node
                    hasSolidNeig = true;
            }
            if (hasSolidNeig)  ret.push_back(n);
        }
    }
    return ret;
}

// template <typename BndType, typename DXQY>

// template <typename DXQY>
template <template <class> class T,  typename DXQY>
T<DXQY> makeFluidBoundary(const int &myRank, const Grid<DXQY> &grid)
/* Sets up Bounce back bounray Boundary object by classifying and adding boundary nodes
 *
 * myRank    : rank of the current process
 * grid      : object of the Grid class
 */

{
    std::vector<int> bndNodes = findFluidBndNodes(myRank, grid);  // Find list of boundary nodes
    T<DXQY> bnd(bndNodes.size());  // Set size of the boundary node object

    for (auto nodeNo: bndNodes) {  // For all boundary nodes
        int nBeta = 0, beta[DXQY::nDirPairs_];
        int nGamma = 0, gamma[DXQY::nDirPairs_];
        int nDelta = 0, delta[DXQY::nDirPairs_];

        for (int q = 0; q < DXQY::nDirPairs_; ++q) {
            int neig_q = grid.neighbor(q, nodeNo);
            int neig_q_rev = grid.neighbor(bnd.dirRev(q), nodeNo);

            if (grid.getType(neig_q) != 0) { // FLUID
                if (grid.getType(neig_q_rev) != 0) { // FLUID :GAMMA
                    gamma[nGamma] = q;
                    nGamma += 1;
                }
                else { // SOLID :BETA (q) and BETA_HAT (qRev)
                    beta[nBeta] = q;
                    nBeta += 1;
                }
            }
            else { // SOLID
                if (grid.getType(neig_q_rev) != 0) { // FLUID :BETA (qDir) and BETA_HAT (q)
                    beta[nBeta] = bnd.dirRev(q);
                    nBeta += 1;
                }
                else { // SOLID: DELTA
                    delta[nDelta] = q;
                    nDelta += 1;
                }
            }
        } // END FOR DIR PAIRS
        bnd.addNode(nodeNo, nBeta, beta, nGamma, gamma, nDelta, delta);
    } // For all boundary nodes
    return bnd;
}



inline int my_mod(int p, int n)
{
    return (p + n) % n;
}


// GEOMETRY FUNCTIONS
// -- defined in LBgeometry.cpp
void newGeometry(const int nX, const int nY, int** &geo);
void newGeometry(const int nX, const int nY, const int nZ, int*** &geo); // 3D
void deleteGeometry(const int nX, const int nY, int** &geo);
void deleteGeometry(const int nX, const int nY, const int nZ, int*** &geo); // 3D
void inputGeometry(int nX, int nY, int** &geo);
void inputGeometry(int nX, int nY, int nZ, int*** &geo); // 3D

// -- template functions
template<typename DXQY>
void analyseGeometry(const int nX, const int nY, int** &geo)
/* Defines nodes as fluid bulk, fluid boundary, solid bulk, solid boundary for a 2d geometry.
 * Here it is assumend that the geometry is periodic.
 *
 *  fluid bulk     : a fluid node with only fluid nodes in the lattice neighborhood
 *  fluid boundary : a fluid node with at least one solid node in the lattice neighborhood
 *  fluid unknown  : not yet checked fluid node
 *  solid bulk     : a solid node with only solid nodes in the lattice neighborhood
 *  solid boundary : a solid node with at least one fluid node in the lattice neighborhood
 *  solid unknown  : not yet checked solid node
 *
 *  integer labels :
 *      FLUID < 3
 *      bulk fluid     : 0
 *      boundary fluid : 1
 *      fluid unknown  : 2
 *
 *      SOLID > 2
 *      boundary solid : 3
 *      bulk solid     : 4
 *      solid unknown  : 5
 *
 *
 * usage : analyseGeometry<D2Q9>(nX, nY, geo)
 *
 * The geo-matrix should now only contain {0, 1} for fluid or {3, 4} for solid
 *
 */
{
    // In the beginning all values are unkonwn
    // fluid nodes are then set to 2 and solid nodes are set to 5
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            int trans[] = {2, 5};  // trans[0(fluid)] = 2 and trans[1(solid)] = 5
            geo[y][x] = trans[geo[y][x]];
        }

    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            /* Calculate the number of fluid neighbors */
            int nFluidNeig = 0;
            for (int q = 0; q < DXQY::nQNonZero_; ++q) {
                int nx, ny;
                nx = my_mod(x + DXQY::c(q, 0), nX);
                ny = my_mod(y + DXQY::c(q, 1), nY);
                nFluidNeig += geo[ny][nx] < 3;
            }

            if (geo[y][x] < 3) { // node is fluid
                if (nFluidNeig < DXQY::nQNonZero_)  geo[y][x] = 1; // (boundary fluid) neighborhood contains solid nodes
                else geo[y][x] = 0;  // (bulk fluid)
            }
            else {  // node is solid
                if (nFluidNeig > 0)  geo[y][x] = 3; // (boundary solid) neighborhood contains fluid nodes
                else  geo[y][x] = 4; // (bulk solid)
            }
        }
}

// -- template functions
template<typename DXQY>
void analyseGeometry(const int nX, const int nY, const int nZ, int*** &geo) // 3D
/* Defines nodes as fluid bulk, fluid boundary, solid bulk, solid boundary for a 3d geometry.
 * Here it is assumend that the geometry is periodic.
 *
 *  fluid bulk     : a fluid node with only fluid nodes in the lattice neighborhood
 *  fluid boundary : a fluid node with at least one solid node in the lattice neighborhood
 *  fluid unknown  : not yet checked fluid node
 *  solid bulk     : a solid node with only solid nodes in the lattice neighborhood
 *  solid boundary : a solid node with at least one fluid node in the lattice neighborhood
 *  solid unknown  : not yet checked solid node
 *
 *  integer labels :
 *      FLUID < 3
 *      bulk fluid     : 0
 *      boundary fluid : 1
 *      fluid unknown  : 2
 *
 *      SOLID > 2
 *      boundary solid : 3
 *      bulk solid     : 4
 *      solid unknown  : 5
 *
 *
 * usage : analyseGeometry<DXQY>(nX, nY, nZ, geo)
 *
 * The geo-matrix should now only contain {0, 1} for fluid or {3, 4} for solid
 *
 */
{
    // In the beginning all values are unkonwn
    // fluid nodes are then set to 2 and solid nodes are set to 5
    for (int z = 0; z < nZ; ++z)
        for (int y = 0; y < nY; ++y)
            for (int x = 0; x < nX; ++x)
            {
                int trans[] = {2, 5};  // trans[0(fluid)] = 2 and trans[1(solid)] = 5
                geo[z][y][x] = trans[geo[z][y][x]];
            }
    for (int z = 0; z < nZ; ++z)
        for (int y = 0; y < nY; ++y)
            for (int x = 0; x < nX; ++x)
            {
                /* Calculate the number of fluid neighbors */
                int nFluidNeig = 0;
                for (int q = 0; q < DXQY::nQNonZero_; ++q) {
                    int nx, ny, nz;
                    nx = my_mod(x + DXQY::c(q, 0), nX);
                    ny = my_mod(y + DXQY::c(q, 1), nY);
                    nz = my_mod(z + DXQY::c(q, 2), nZ);
                    nFluidNeig += geo[nz][ny][nx] < 3;
                }

                if (geo[z][y][x] < 3) { // node is fluid
                    if (nFluidNeig < DXQY::nQNonZero_)  geo[z][y][x] = 1; // (boundary fluid) neighborhood contains solid nodes
                    else geo[z][y][x] = 0;  // (bulk fluid)
                }
                else {  // node is solid
                    if (nFluidNeig > 0)  geo[z][y][x] = 3; // (boundary solid) neighborhood contains fluid nodes
                    else  geo[z][y][x] = 4; // (bulk solid)
                }
            }
}


// LABEL FUNCTIONS
// -- defined in LBgeometry.cpp
void newNodeLabel(int nX, int nY, int** &nodeLabel);
void newNodeLabel(int nX, int nY, int nZ, int*** &nodeLabel); // 3D
void deleteNodeLabel(int nX, int nY, int** &nodeLabel);
void deleteNodeLabel(int nX, int nY, int nZ, int*** &nodeLabel); // 3D
int setBulkLabel(int nX, int nY, int** geo, int** &nodeLabel);
int setBulkLabel(int nX, int nY, int nZ, int*** geo, int*** &nodeLabel); // 3D
int setNonBulkLabel(int previousLabel, int nX, int nY, int** geo, int** &nodeLabel);
int setNonBulkLabel(int previousLabel, int nX, int nY, int nZ, int*** geo, int*** &nodeLabel); // 3D

// -- Template functions
template<typename DXQY>
void setupGrid(int nX,  int nY, int ** nodeLabel,  Grid<DXQY> &grid)
/*Sets up Grid object by adding active nodes and the respective neighbor nodes to the object
 *
 * nX        : number of grid points in the Cartesian x-direction
 * nY        : number of grid points in the Cartesian y-direction
 * nodeLabel : pointer to the nodeLabel matrix
 * grid      : object of the Grid class that is to be set up 
 */
{
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (nodeLabel[y][x] > 0) {
                int nodeNo = nodeLabel[y][x];
                grid.addNodePos(x, y, nodeNo); //LBgrid
                for (int q = 0; q < DXQY::nQ; ++q) {
                    int nx, ny;
                    nx = my_mod(x + DXQY::c(q, 0), nX);
                    ny = my_mod(y + DXQY::c(q, 1), nY);
                    grid.addNeigNode(q, nodeNo, nodeLabel[ny][nx]); //LBgrid
                }
            }
        }
}

template<typename DXQY>
void setupGrid(int nX,  int nY,  int nZ, int *** nodeLabel,  Grid<DXQY> &grid) // 3D
/*Sets up Grid object by adding active nodes and the respective neighbor nodes to the object
 *
 * nX        : number of grid points in the Cartesian x-direction
 * nY        : number of grid points in the Cartesian y-direction
 * nZ        : number of grid points in the Cartesian z-direction
 * nodeLabel : pointer to the nodeLabel matrix
 * grid      : object of the Grid class that is to be set up 
 */
{
  for (int z = 0; z < nZ; ++z)
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (nodeLabel[z][y][x] > 0) {
                int nodeNo = nodeLabel[z][y][x];
                grid.addNodePos(x, y, z, nodeNo); //LBgrid
                //grid.addNodePos({x, y, z}, nodeNo); //LBgrid
                for (int q = 0; q < DXQY::nQ; ++q) {
                    int nx, ny, nz;
                    nx = my_mod(x + DXQY::c(q, 0), nX);
                    ny = my_mod(y + DXQY::c(q, 1), nY);
                    nz = my_mod(z + DXQY::c(q, 2), nZ);
                    grid.addNeigNode(q, nodeNo, nodeLabel[nz][ny][nx]); //LBgrid
                }
            }
        }
}


// -- Inline functions
inline int nBoundaryNodes(int bndLabel, int nX, int nY, int** &geo)
{
    int nNodes = 0;
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (geo[y][x] == bndLabel)
                nNodes += 1;
        }
    return nNodes;
}

inline int nBoundaryNodes(int bndLabel, int nX, int nY, int nZ, int*** &geo) // 3D
{
    int nNodes = 0;
    for (int z = 0; z < nZ; ++z)
      for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (geo[z][y][x] == bndLabel)
                nNodes += 1;
        }
    return nNodes;
}

inline void setupBulk(int nX, int nY, int** &geo, int** &nodeLabel, Bulk &bulk)
/* Sets up Bulk object by adding active bulk nodes to the list of bulk nodes
 *
 * nX        : number of grid points in the Cartesian x-direction
 * nY        : number of grid points in the Cartesian y-direction
 * geo       : pointer to the geo matrix
 * bulk      : object of the bulk class that is to be set up 
 */
{
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (geo[y][x] < 2)
                bulk.addBulkNode(nodeLabel[y][x]);
        }

}

inline void setupBulk(int nX, int nY, int nZ, int*** &geo, int*** &nodeLabel, Bulk &bulk) // 3D
/* Sets up Bulk object by adding active bulk nodes to the list of bulk nodes
 *
 * nX        : number of grid points in the Cartesian x-direction
 * nY        : number of grid points in the Cartesian y-direction
 * nZ        : number of grid points in the Cartesian z-direction
 * geo       : pointer to the geo matrix
 * bulk      : object of the bulk class that is to be set up 
 */
{
   for (int z = 0; z < nZ; ++z)
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (geo[z][y][x] < 2)
                bulk.addBulkNode(nodeLabel[z][y][x]);
        }

}




#endif // LBGEOMETRY_H
