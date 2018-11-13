#ifndef LBFIELD_H
#define LBFIELD_H

#include "LBlattice.h"

template <typename BASETYPE>
class Field
{
public:
    Field(int elementSize, int nElements);
    ~Field();

    BASETYPE& operator () (int posInElement, int elementNo);
    BASETYPE* operator () (int elementNo); // Returns a pointer to the first position in elementNo

    void swapData(Field& field);

private:
    int elementSize_;
    int nElements_;
    BASETYPE* data_;
};

// Seems that template class construcutors and destructors need to be
// included in the header-file.
template <typename BASETYPE>
Field<BASETYPE>::Field(int elementSize, int nElements)
{
    elementSize_ = elementSize;
    nElements_ = nElements;
    data_ = new BASETYPE [elementSize_ * nElements_];
}

template <typename BASETYPE>
Field<BASETYPE>::~Field()
{
    delete [] data_;
}

template <typename BASETYPE>
inline BASETYPE& Field<BASETYPE>::operator () (int posInElement, int elementNo)
{
    return data_[elementNo * elementSize_ + posInElement];  // Here we have choosen a collision based data structure
}

template <typename BASETYPE>
inline BASETYPE* Field<BASETYPE>::operator () (int elementNo) // Returns a pointer to the first position in elementNo
{
    return &data_[elementNo * elementSize_];
}

template <typename BASETYPE>
inline void Field<BASETYPE>::swapData(Field<BASETYPE>& field)
{
    BASETYPE* dataTmp;
    dataTmp = field.data_;
    field.data_ = data_;
    data_ = dataTmp;
}



template <typename BASETYPE>
class LbField: public Field<BASETYPE> // public Field<double>
{
public:
    LbField(int elementSize, int nElements, Lattice* lattice);

    void zerothMoment(int elementNo, double* returnScalar);
    void firstMoment(int elementNo, double* returnVector);

private:
    Lattice *lattice_;
};

template <typename BASETYPE>
LbField<BASETYPE>::LbField(int elementSize, int nElements, Lattice* lattice) : Field<BASETYPE>(elementSize, nElements)
{
    lattice_ = lattice;
}

template <typename BASETYPE>
inline void LbField<BASETYPE>::zerothMoment(int elementNo, double* returnScalar)
{
    (*returnScalar) = 0.0;
    for (int q = 0; q < (*lattice_).nQ(); q++) {
        (*returnScalar) += (*this)(q, elementNo);
    }
}

template <typename BASETYPE>
inline void LbField<BASETYPE>::firstMoment(int elementNo, double* returnVector)
{
    for (int d = 0; d < (*lattice_).nD(); d++)
        returnVector[d] = (*lattice_).innerProductQMajor(d, (*this)(elementNo));
}




#endif // LBFIELD_H
