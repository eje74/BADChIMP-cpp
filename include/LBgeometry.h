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



template<typename DXQY>
void setupGrid(int nX,  int nY, int ** nodeLabel,  Grid &grid)
{
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (nodeLabel[y][x] > 0) {
                int nodeNo = nodeLabel[y][x];
                grid.addNodePos(x, y, nodeNo);
                for (int q = 0; q < DXQY::nQ; ++q) {
                    int nx, ny;
                    nx = my_mod(x + DXQY::c(q, 0), nX);
                    ny = my_mod(y + DXQY::c(q, 1), nY);
                    grid.addNeigNode(q, nodeNo, nodeLabel[ny][nx]);
                }
            }
        }
}

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
{
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (geo[y][x] < 2)
                bulk.addBulkNode(nodeLabel[y][x]);
        }

}

template <typename DXQY>
void setupBoundary(int bndLabel, int nX, int nY, int** &geo, int** &nodelLabel, Grid &grid, Boundary<DXQY> &bnd)
{
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {
            if (geo[y][x] == bndLabel) {
                int node = nodelLabel[y][x];
                int nBeta = 0, beta[D2Q9::nDirPairs_];
                int nGamma = 0, gamma[D2Q9::nDirPairs_];
                int nDelta = 0, delta[D2Q9::nDirPairs_];

                for (int q = 0; q < D2Q9::nDirPairs_; ++q) {
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

// Here we have assumed that all boundaries are periodic.
template<typename DXQY>
void analyseGeometry(const int nX, const int nY, int** &geo)
{
    // BULK FLUID     : 0
    // BOUNDARY FLUID : 1
    // FLUID UNKNOWN  : 2

    // BOUNDARY SOLID : 3
    // BULK SOLID     : 4
    // BULK UNKNOWN   : 5

    // In the begning all values are unkonwn
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            int trans[] = {2, 5};
            geo[y][x] = trans[geo[y][x]];
        }

    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            /* Calculate the number of fluid neighbors */
            int nFluidNeig = 0;
            for (int q = 0; q < DXQY::nQ-1; ++q) {
                int nx, ny;
                nx = my_mod(x + DXQY::c(q, 0), nX);
                ny = my_mod(y + DXQY::c(q, 1), nY);
                nFluidNeig += geo[ny][nx] < 3;
            }
            if (geo[y][x] < 3) { // FLUID
                if (nFluidNeig < DXQY::nQ-1)  geo[y][x] = 1;
                else geo[y][x] = 0;
            }
            else {  // SOLID
                if (nFluidNeig > 0)  geo[y][x] = 3;
                else  geo[y][x] = 4;
            }
        }
}


void makeGeometry(const int nX, const int nY, int** &geo);
void deleteGeometry(const int nX, const int nY, int** &geo);
void printGeoScr(const int nX, const int nY, int** geo);

void makeNodeLabel(int nX, int nY, int** &nodeLabel);
void deleteNodeLabel(int nX, int nY, int** &nodeLabel);

int setBulkLabel(int nX, int nY, int** geo, int** &nodeLabel);
int setNonBulkLabel(int firstLabel, int nX, int nY, int** geo, int** &nodeLabel);




#endif // LBGEOMETRY_H
