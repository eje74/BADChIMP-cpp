#include "LBfield.h"


ScalarField::ScalarField(const int nFields, const int nNodes): nFields_(nFields), nNodes_(nNodes)
{
    data_ = new lbBase_t [nFields_ * nNodes_];
}

ScalarField::~ScalarField()
{
    delete [] data_;
}
