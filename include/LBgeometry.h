#ifndef LBGEOMETRY_H
#define LBGEOMETRY_H

#include "LBglobal.h"
#include "LBgrid.h"

static int my_mod(int p, int n)
{
    return (p + n) % n;
}

template<typename DXQY>
void setupGridNeig(int nX,  int nY, int ** geo, int ** nodeLabel,  Grid &grid)
{
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x) {

        }

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


void makeGeometryX(const int nX, const int nY, int** &geo);
void deleteGeometryX(const int nX, const int nY, int** &geo);
void printGeoScr(const int nX, const int nY, int** geo);



void makeNodeLabel(int nX, int nY, int** &nodeLabel);
void deleteNodeLabel(int nX, int nY, int** &nodeLabel);

int setBulkLabel(int nX, int nY, int** geo, int** &nodeLabel);
int setNonBulkLabel(int firstLabel, int nX, int nY, int** geo, int** &nodeLabel);



#endif // LBGEOMETRY_H
