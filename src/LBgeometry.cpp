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
 * nX and nY is the matrix size, number of columns and
 *  number of rows, respectively
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


int setBulkLabel(int nX, int nY, int** geo, int** &nodeLabel)
/* Tags the bulk nodes in the nodeLabel matrix. Here bulk means
 *  fluid nodes means nodes that uses the standard LB update
 *  algorithm
 *
 * Input arguments:
 *  nX, nY: system size
 *  geo: geometry matrix analyzed by function 'analyseGeometry'
 *  nodeLabel: matrix crated by function newNodeLabel
 *
 * Description:
 *  - The bulk labels are number from 1 to the total number of
 *    bulk nodes.
 *  - 0 is used as a default node label, assumed to be a dummy
 *    variable
 *  - returns an integer that is the highest bulk tag.
 *  - Assumes that geo is initialized with zeros so that the
 *    bulk nodes are taged with unique labels ranging from 1 to the
 *    total number of bulk nodes, ie. the highest bulk tag.
 */
{
    int label = 0;
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            if (geo[y][x] < 3) { // Here we have set all fluid nodes to bulk nodes
                label += 1;
                nodeLabel[y][x] = label;
            }
        }
    return label;
}


int setNonBulkLabel(int previousLabel, int nX, int nY, int** geo, int** &nodeLabel)
/* Tags the non-bulk nodes in the nodeLabel matrix that is
 *  part of the computational domain. Here non bulk means nodes
 *  that do not use the standard LB update algorithm.
 *
 * Input arguments:
 *  previousLabel: last label set
 *  nX, nY: system size
 *  geo: geometry matrix analyzed by function 'analyseGeometry'
 *  nodeLabel: matrix crated by function newNodeLabel
 *
 * Description:
 * - Nodes are taged with labels beginning with [previoiusLabel + 1]
 * - returns the integer value of the last tag assigned.
 */
{
    int label = previousLabel;
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
        {
            if (geo[y][x] == 3) {
                assert(nodeLabel[y][x] == 0); // Label should not have been previously set
                label += 1;
                nodeLabel[y][x] = label;
            }
        }
    return label;
}
/*======================================================*/


