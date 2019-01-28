#include "LBfield.h"


ScalarField::ScalarField(const int nFields, const int nNodes): nFields_(nFields), nNodes_(nNodes)
{
    data_ = new lbBase_t [nFields_ * nNodes_];
}

ScalarField::~ScalarField()
{
    delete [] data_;
}



VectorField::VectorField(const int nFields, const int nDimensions, const int nNodes)
    :nFields_(nFields), nDim_(nDimensions), elementSize_(nFields_ * nDim_), nNodes_(nNodes)
{
    data_ = new lbBase_t [elementSize_ * nNodes_];
}

VectorField::~VectorField()
{
    delete [] data_;
}



LbField::LbField(const int nFields, const int nDirections, const int nNodes)
    :nFields_(nFields), nDir_(nDirections), elementSize_(nFields_ * nDir_), nNodes_(nNodes)
{
    data_ = new lbBase_t [elementSize_ * nNodes_];
}

LbField::~LbField()
{
    delete [] data_;
}

