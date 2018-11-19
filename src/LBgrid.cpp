#include "LBgrid.h"

// Constructor and destructor
GridRegular::GridRegular(const int nX, const int nY, Lattice *lattice): nX_(nX + 2), nY_(nY + 2), nElements_(nX_ * nY_), lattice_(lattice)  // Input can also be one or many files ...
{
    //lattice_ = lattice;
    // Setup neighbor list
    for (int q = 0; q < (*lattice_).nQ(); ++q)
        neighborStride[q] = (*lattice_).c(q, 0) + nX_ * (*lattice_).c(q, 1);
}
GridRegular::~GridRegular()
{}
