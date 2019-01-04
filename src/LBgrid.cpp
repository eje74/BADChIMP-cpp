#include "LBgrid.h"

Grid::Grid(const int maxNodeNo, const int stride) :maxNodeNo_(maxNodeNo), stride_(stride)
{
    neigList_ = new int [maxNodeNo_ * stride_]; // Data structure [nodeNo*stride_ + direction]
}

Grid::~Grid()
{
    delete [] neigList_;
}

void Grid::addElement(const int qDirection, const int nodeNo, const int nodeNeigNo)
{
    neigList_[nodeNo * stride_ + qDirection] = nodeNeigNo;
}
