#include "LBgrid.h"

Grid::Grid(const int maxNodeNo, const int stride) :maxNodeNo_(maxNodeNo), stride_(stride)
{
    // Needs to add one to include a node with label = maxNodeNo_
    neigList_ = new int [(maxNodeNo_ + 1) * stride_];
}

Grid::~Grid()
{
    delete [] neigList_;
}

void Grid::addElement(const int qDirection, const int nodeNo, const int nodeNeigNo)
{
    neigList_[nodeNo * stride_ + qDirection] = nodeNeigNo;
}
