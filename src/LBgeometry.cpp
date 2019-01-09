#include "LBgeometry.h"
#include <iostream>
#include <iomanip>
#include <assert.h>


/* We need to specify for which nodes that we need to allocate memory, and
   for which nodes we need to get the link information.
  */
void setupGridNeig(const int nLinks, const int* links,  Grid &grid)
{
    for (int n = 0; n < nLinks; ++n)
    {
        int num = 3 * n;
        //               direction    from node       to node
        grid.addElement(links[num], links[num + 1], links[num + 2]);
    }
}


void makeNodeLabel(int nX, int nY, int** &nodeLabel)
{
    nodeLabel = new int* [nY];
    for (int y = 0; y < nY; ++y)
    {
        nodeLabel[y] = new int [nX];
        for (int x = 0; x < nX; ++x)
        {
            nodeLabel[y][x] = 0;
        }
    }
}

void deleteNodeLabel(int nX, int nY, int** &nodeLabel)
{
    for (int y = 0; y < nY; ++y)
        delete []  nodeLabel[y];
    delete [] nodeLabel;
    nodeLabel = nullptr;

}


/* We will number bulk labels from 1 to "#of bulk labels" */
int setBulkLabel(int nX, int nY, int** geo, int** &nodeLabel)
{
    int label = 0;
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            if (geo[y][x] < 3) {
                label += 1;
                nodeLabel[y][x] = label;
            }
        }
    return label;
}


int setNonBulkLabel(int firstLabel, int nX, int nY, int** geo, int** &nodeLabel)
{
    int label = firstLabel;
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            if (geo[y][x] == 3) {
                assert(nodeLabel[y][x] == 0); // Label should not be previously set
                label += 1;
                nodeLabel[y][x] = label;
            }
        }
    return label;
}



void makeGeometryX(const int nX, const int nY, int** &geo)
{
    geo = new int* [nY];
    for (int y = 0; y < nY; ++y)
        geo[y] = new int [nX];

    // MAKE GEOMETRY
    // FLUID = 0
    // SOLID = 1
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            if (y < 3)
                geo[y][x] = 1;
            else
                geo[y][x] = 0;
        }
}


void deleteGeometryX(const int nX, const int nY, int** &geo)
{
    for (int y = 0; y < nY; ++y)
        delete []  geo[y];
    delete [] geo;
    geo = nullptr;
}


void printGeoScr(const int nX, const int nY, int** geo)
{
    std::cout << " ";
    for (int x = 0; x < nX; ++x)
        std::cout << "-";
    std::cout << std::endl;
    for (int y = nY-1; y >= 0; --y) {
       std::cout << "|";
       for (int x = 0; x < nX; ++x) {
           std::cout << std::setw(3) << geo[y][x];
       }
       std::cout << "|" << std::endl;
    }
    std::cout << " ";
    for (int x = 0; x < nX; ++x)
        std::cout << "-";
    std::cout << std::endl;

}

