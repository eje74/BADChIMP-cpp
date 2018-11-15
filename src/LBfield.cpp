#include "LBfield.h"


ScalarField::ScalarField(const int nFields, const int nElements): nFields_(nFields), nElements_(nElements)
{
    data_ = new lbbase_t [nFields_ * nElements_];
}

ScalarField::~ScalarField()
{
    delete [] data_;
}



VectorField::VectorField(const int nFields, const int nDimensions, const int nElements)
    :nFields_(nFields), nDimensions_(nDimensions), elementSize_(nFields_ * nDimensions_), nElements_(nElements)
{
    data_ = new lbbase_t [elementSize_ * nElements_];
}

VectorField::~VectorField()
{
    delete [] data_;
}



LbField::LbField(const int nFields, const int nDirections, const int nElements)
    :nFields_(nFields), nDirections_(nDirections), elementSize_(nFields_ * nDirections_), nElements_(nElements)
{
    data_ = new lbbase_t [elementSize_ * nElements_];
}

LbField::~LbField()
{
    delete [] data_;
}



GeoField::GeoField(const int nFields, const int nElements): nFields_(nFields), nElements_(nElements)
{
    data_ = new int [nFields_ * nElements_];
}

GeoField::~GeoField()
{
    delete [] data_;
}

