#include "LBgeometry.h"
#include <iostream>
#include <iomanip>
#include <assert.h>


/* ****************************************************
 * Functions for creating, filling and deleting a 2d
 * structured GEOMETRY matrix
 * ****************************************************/
void newGeometry(const int nX, const int nY, int** &geo)
/* geo : is a pointer variable
 *       the value of geo at point (x,y) is geo[y][x]
 * nX, and nY are the spatial dimensions */
{
    geo = new int* [nY];
    for (int y = 0; y < nY; ++y)
        geo[y] = new int [nX];
}


void deleteGeometry(const int nX, const int nY, int** &geo)
{
    for (int y = 0; y < nY; ++y)
        delete []  geo[y];
    delete [] geo;
    geo = nullptr;
}

void inputGeometry(int nX, int nY, int** &geo)
/*  Makes a geometry:
 *  The conventions is that:
 *    FLUID = 0
 *    SOLID = 1
 * This function could be replaced by a geo file
 */
{
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            if (y == 0)
                geo[y][x] = 1;
            else
                geo[y][x] = 0;
        }
}
/*======================================================*/


/*******************************************************
 * Node LABEL functions. Each node in the computational
 * domains should be given its unique label-tag (int number)
 *
 *  Example:
 *  geo-matrix:     label-matrix:
 *    4 4 4            0  0  0
 *    3 3 3           10 11 12
 *    1 1 1            1  2  3
 *    0 0 0            4  5  6
 *    1 1 1            7  8  9
 *    3 3 3           13 14 15
 *
 *  Here we have made the choices that
 *  - 0 is a default bulk solid tag. We will allocate
 *    memory for it, but we assume that it contains
 *    only dummy values
 *  - We also want the bulk node values to be numbered
 *    consequently from 1 (number one). In this example
 *    both fluid bulk and fluid boundary are treated as
 *    bulk nodes in the LB iterations for macroscopic
 *    varaible evaluation, collision and propagation.
 *
 *******************************************************/
void newNodeLabel(int nX, int nY, int** &nodeLabel)
/* nodeLabel is the tag-matrix.
 *  the tag for the node at point (x,y) is nodeLabel[y][x]
 * nX and nY is the matrix size
 *
 * The matrix is initiated to zero.
 */
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


// EJE START HER

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

