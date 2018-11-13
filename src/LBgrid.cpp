#include "LBgrid.h"

// Constructor and destructor
GridRegular::GridRegular(const int nX, const int nY, Lattice *lattice): nX_(nX + 2), nY_(nY + 2), nElements_(nX_ * nY_)  // Input can also be one or many files ...
{
    lattice_ = lattice;
}
GridRegular::~GridRegular()
{}
