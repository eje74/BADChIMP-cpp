#include "LBgrid.h"

Grid::Grid(const int maxNodeNo, const int stride) :maxNodeNo_(maxNodeNo), stride_(stride)
{
    // Needs to add one to include a node with label = maxNodeNo_
    neigList_ = new int [(maxNodeNo_ + 1) * stride_];
    pos_ = new double [(maxNodeNo_ + 1) * 2];
}

Grid::~Grid()
{
    delete [] neigList_;
}

void Grid::addNeigNode(const int qDirection, const int nodeNo, const int nodeNeigNo)
{
    neigList_[nodeNo * stride_ + qDirection] = nodeNeigNo;
}

void Grid::addNodePos(const double x, const double y, const int nodeNo)
{
    pos_[2*nodeNo] = x;
    pos_[2*nodeNo + 1] = y;
}

