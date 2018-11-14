#ifndef LBFIELD_H
#define LBFIELD_H

#include "LBlattice.h"

// SCALARFIELD
template <typename BASETYPE>
class ScalarField
{
public:
    ScalarField(const int nFields, const int nElements);
    ~ScalarField();

    BASETYPE& operator () (const int fieldNo, const int elementNo) const;

private:
    const int nFields_;
    int nElements_;
    BASETYPE* data_;
};

template <typename BASETYPE>
ScalarField<BASETYPE>::ScalarField(const int nFields, const int nElements): nFields_(nFields), nElements_(nElements)
{
    data_ = new BASETYPE [nFields_ * nElements_];
}

template <typename BASETYPE>
ScalarField<BASETYPE>::~ScalarField()
{
    delete [] data_;
}

template <typename BASETYPE>
inline BASETYPE& ScalarField<BASETYPE>::operator () (const int fieldNo, const int elementNo) const
{
    return data_[nFields_ * elementNo + fieldNo];
}
// END SCALARFIELD



// VECTORFIELD
template <typename BASETYPE>
class VectorField
{
public:
    VectorField(const int nFields, const int nDimensions, const int nElements);
    ~VectorField();

    BASETYPE& operator () (const int fieldNo, const int dimNo, const int elementNo) const; // Returns element
    BASETYPE* operator () (const int fieldNo, const int elementNo) const; // Returns pointer to beginning of a vector

private:
    const int nFields_;
    const int nDimensions_;
    const int elementSize_;
    int nElements_;
    BASETYPE* data_;
};

template <typename BASETYPE>
VectorField<BASETYPE>::VectorField(const int nFields, const int nDimensions, const int nElements)
    :nFields_(nFields), nDimensions_(nDimensions), elementSize_(nFields_ * nDimensions_), nElements_(nElements)
{
    data_ = new BASETYPE [elementSize_ * nElements_];
}

template <typename BASETYPE>
VectorField<BASETYPE>::~VectorField()
{
    delete [] data_;
}

template <typename BASETYPE>
inline BASETYPE& VectorField<BASETYPE>::operator () (const int fieldNo, const int dimNo, const int elementNo) const
{
    return data_[elementSize_ * elementNo + nDimensions_ * fieldNo + dimNo];
}

template <typename BASETYPE>
inline BASETYPE* VectorField<BASETYPE>::operator () (const int fieldNo, const int elementNo) const
{
    return &data_[elementSize_ * elementNo + nDimensions_ * fieldNo];
}
// END VECTORFIELD


// LBFIELD
template <typename BASETYPE>
class LbField
{
public:
    LbField(const int nFields, const int nDirections_, const int nElements);
    ~LbField();

    BASETYPE& operator () (const int fieldNo, const int dirNo, const int elementNo) const; // Returns element
    BASETYPE* operator () (const int fieldNo, const int elementNo) const; // Returns pointer to beginning of a vector

    inline void swapData(LbField<BASETYPE>& field);

private:
    const int nFields_;
    const int nDirections_;
    const int elementSize_;
    int nElements_;
    BASETYPE* data_;
};

template <typename BASETYPE>
LbField<BASETYPE>::LbField(const int nFields, const int nDirections, const int nElements)
    :nFields_(nFields), nDirections_(nDirections), elementSize_(nFields * nDirections), nElements_(nElements)
{
    data_ = new BASETYPE [elementSize_ * nElements_];
}

template <typename BASETYPE>
LbField<BASETYPE>::~LbField()
{
    delete [] data_;
}

template <typename BASETYPE>
inline BASETYPE& LbField<BASETYPE>::operator () (const int fieldNo, const int dirNo, const int elementNo) const // Returns element
{
    return data_[elementSize_ * elementNo + nDirections_ * fieldNo + dirNo];
}

template <typename BASETYPE>
inline BASETYPE* LbField<BASETYPE>::operator () (const int fieldNo, const int elementNo) const // Returns pointer to beginning of a vector
{
    return &data_[elementSize_ * elementNo + nDirections_ * fieldNo];
}

template <typename BASETYPE>
inline void LbField<BASETYPE>::swapData(LbField<BASETYPE>& field)
{
    BASETYPE* dataTmp;
    dataTmp = field.data_;
    field.data_ = data_;
    data_ = dataTmp;
}

// END LBFIELD



#endif // LBFIELD_H
