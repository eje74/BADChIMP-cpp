#ifndef LBFIELD_H
#define LBFIELD_H

#include "LBglobal.h"
#include "LBlattice.h"

// SCALARFIELD

class ScalarField
{
public:
    ScalarField(const int nFields, const int nElements);
    ~ScalarField();

    lbbase_t& operator () (const int fieldNo, const int elementNo) const;

private:
    const int nFields_;
    int nElements_;
    lbbase_t* data_;
};


inline lbbase_t& ScalarField::operator () (const int fieldNo, const int elementNo) const
{
    return data_[nFields_ * elementNo + fieldNo];
}
// END SCALARFIELD



// VECTORFIELD

class VectorField
{
public:
    VectorField(const int nFields, const int nDimensions, const int nElements);
    ~VectorField();

    lbbase_t& operator () (const int fieldNo, const int dimNo, const int elementNo) const; // Returns element
    lbbase_t* operator () (const int fieldNo, const int elementNo) const; // Returns pointer to beginning of a vector

private:
    const int nFields_;
    const int nDimensions_;
    const int elementSize_;
    int nElements_;
    lbbase_t* data_;
};


inline lbbase_t& VectorField::operator () (const int fieldNo, const int dimNo, const int elementNo) const
{
    return data_[elementSize_ * elementNo + nDimensions_ * fieldNo + dimNo];
}


inline lbbase_t* VectorField::operator () (const int fieldNo, const int elementNo) const
{
    return &data_[elementSize_ * elementNo + nDimensions_ * fieldNo];
}
// END VECTORFIELD


// LBFIELD

class LbField
{
public:
    LbField(const int nFields, const int nDirections_, const int nElements);
    ~LbField();

    lbbase_t& operator () (const int fieldNo, const int dirNo, const int elementNo) const; // Returns element
    lbbase_t* operator () (const int fieldNo, const int elementNo) const; // Returns pointer to beginning of a vector

    inline void swapData(LbField& field);

private:
    const int nFields_;
    const int nDirections_;
    const int elementSize_;
    int nElements_;
    lbbase_t* data_;
};



inline lbbase_t& LbField::operator () (const int fieldNo, const int dirNo, const int elementNo) const // Returns element
{
    return data_[elementSize_ * elementNo + nDirections_ * fieldNo + dirNo];
}


inline lbbase_t* LbField::operator () (const int fieldNo, const int elementNo) const // Returns pointer to beginning of a vector
{
    return &data_[elementSize_ * elementNo + nDirections_ * fieldNo];
}


inline void LbField::swapData(LbField& field)
{
    lbbase_t* dataTmp;
    dataTmp = field.data_;
    field.data_ = data_;
    data_ = dataTmp;
}

// END LBFIELD


// GEOFIELD
class GeoField
{
public:
    GeoField(const int nFields, const int nElements);
    ~GeoField();

    int& operator () (const int fieldNo, const int elementNo) const;

private:
    const int nFields_;
    int nElements_;
    int* data_;

};


inline int& GeoField::operator () (const int fieldNo, const int elementNo) const
{
    return data_[nFields_ * elementNo + fieldNo];
}

// END GEOFIELD

#endif // LBFIELD_H