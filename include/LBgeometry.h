#ifndef LBGEOMETRY_H
#define LBGEOMETRY_H

#include "LBglobal.h"
#include "LBgrid.h"
#include "LBboundary.h"
#include "LBbulk.h"

inline int my_mod(int p, int n)
{
    return (p + n) % n;
}


// GEOMETRY FUNCTIONS
// -- defined in LBgeometry.cpp
void newGeometry(const int nX, const int nY, int** &geo);
void deleteGeometry(const int nX, const int nY, int** &geo);
void inputGeometry(int nX, int nY, int** &geo);
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


// LABEL FUNCTIONS
// -- defined in LBgeometry.cpp
void newNodeLabel(int nX, int nY, int** &nodeLabel);
void deleteNodeLabel(int nX, int nY, int** &nodeLabel);
int setBulkLabel(int nX, int nY, int** geo, int** &nodeLabel);
int setNonBulkLabel(int previousLabel, int nX, int nY, int** geo, int** &nodeLabel);

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

template <typename DXQY>
void setupBoundary(int bndLabel, int nX, int nY, int** &geo, int** &nodelLabel, Grid<DXQY> &grid, Boundary<DXQY> &bnd)
/* Sets up Boundary object by classifying and adding boundary nodes
 *
 * nX        : number of grid points in the Cartesian x-direction
 * nY        : number of grid points in the Cartesian y-direction
 * geo       : pointer to the geo matrix
 * nodeLabel : pointer to the nodeLabel matrix
 * grid      : object of the Grid class
 * bnd       : object of the Boundary class that is to be set up
 */
{
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (geo[y][x] == bndLabel) {
                int node = nodelLabel[y][x];
                int nBeta = 0, beta[DXQY::nDirPairs_];
                int nGamma = 0, gamma[DXQY::nDirPairs_];
                int nDelta = 0, delta[DXQY::nDirPairs_];

                for (int q = 0; q < DXQY::nDirPairs_; ++q) {
                    int* qPos = grid.pos( grid.neighbor(q, node) );
                    int* qRevPos = grid.pos( grid.neighbor(bnd.dirRev(q), node) );

                    if (geo[qPos[1]][qPos[0]] < 3) { // FLUID
                        if (geo[qRevPos[1]][qRevPos[0]] < 3) { // FLUID :GAMMA
                            gamma[nGamma] = q;
                            nGamma += 1;
                        }
                        else { // SOLID :BETA (q) and BETA_HAT (qRev)
                            beta[nBeta] = q;
                            nBeta += 1;
                        }
                    }
                    else { // SOLID
                        if (geo[qRevPos[1]][qRevPos[0]] < 3) { // FLUID :BETA (qDir) and BETA_HAT (q)
                            beta[nBeta] = bnd.dirRev(q);
                            nBeta += 1;
                        }
                        else { // SOLID: DELTA
                            delta[nDelta] = q;
                            nDelta += 1;
                        }
                    }
                } // END FOR DIR PAIRS
                bnd.addNode(node, nBeta, beta, nGamma, gamma, nDelta, delta);
            } // END FOR BOUNDARY NODE

        } // END FOR x TO nX
}


#endif // LBGEOMETRY_H
