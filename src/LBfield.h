#ifndef LBFIELD_H
#define LBFIELD_H

#include "LBglobal.h"
#include <iostream>
#include <algorithm>
#include <vector>

// SCALARFIELD

/*********************************************************
 * class SCALARFIELD: Represents a given number of scalar
 *  fields
 *
 *********************************************************/
class ScalarField
{
public:
    /* Constructor:
     * nFields : number of field
     * nNodes  : number of nodes per field
     */
    ScalarField(const int nFields, const int nNodes): nFields_(nFields), nNodes_(nNodes),
        data_(static_cast<std::size_t>(nFields * nNodes)) {}

    /* operator overloading of (). */
    inline const lbBase_t& operator () (const int fieldNo,const int nodeNo) const;
    inline lbBase_t& operator () (const int fieldNo,const int nodeNo);

    //JLV
    std::vector<std::vector<lbBase_t>::iterator> get_iterator(const int fieldNo, const std::vector<int> nodes)  {
      std::vector<std::vector<lbBase_t>::iterator> itr(nodes.size());
      for (unsigned n=0; n< nodes.size(); ++n)
        itr[n] =  data_.begin() + (nFields_ * nodes[n] + fieldNo);
      return itr;
    }
    //JLV

    /*
     * fieldNo : the current field
     * nodeNo : the current node (tag)
     *
     * Returns the scalar value for the given field at the given node
     * Usage:
     * rho(0, 32) = 3.5 // sets the value of field_0's node 32 to 3.5
     */

    int size() {return nNodes_;} // Getter for nNodes_
    int num_fields() const {return nFields_;}
private:
    const int nFields_;  // Number of fields
    int nNodes_;  // Number of nodes in each field
    std::vector<lbBase_t> data_; // Container for scalar data
};


// Two wersions to handle const corectly
inline const lbBase_t& ScalarField::operator () (const int fieldNo, const int nodeNo) const
{
    return data_[static_cast<std::size_t>(nFields_ * nodeNo + fieldNo)];
}

inline lbBase_t& ScalarField::operator () (const int fieldNo, const int nodeNo)
{
    return data_[static_cast<std::size_t>(nFields_ * nodeNo + fieldNo)];
}


// END SCALARFIELD



/*********************************************************
 * class VECTORFIELD: Represents a given number of vector
 *  fields
 *
 *********************************************************/
/* template <typename DXQY>
class VectorFieldRet
{
public:
    VectorFieldRet(lbBase_t *dstBegin): dstBegin_(dstBegin) {}
    template <typename T>
    void operator = (const std::vector<T> &rhs)
    {
        std::copy(rhs.begin(), rhs.begin() + DXQY::nD, dstBegin_);
    }

private:
    lbBase_t *dstBegin_;
}; */


template <typename DXQY>
class VectorField
{
public:
    /* Constructor
     * nFields     : number of vector fields
     * nDimensions : number of spatial dimensions
     * nNodes      : number of nodes
     */
    VectorField(const int nFields, const int nNodes):
        nFields_(nFields), elementSize_(nFields_ * DXQY::nD), nNodes_(nNodes), data_(nFields * nNodes * DXQY::nD){}


    /* operator overloading of () */
    inline const lbBase_t& operator () (const int fieldNo, const int dimNo, const int nodeNo) const // Returns element
    {
        return data_[elementSize_ * nodeNo + DXQY::nD * fieldNo + dimNo];
    }
    inline lbBase_t& operator () (const int fieldNo, const int dimNo, const int nodeNo) // Returns element
    {
        return data_[elementSize_ * nodeNo + DXQY::nD * fieldNo + dimNo];
    }
    /* Returns a reference to a vector component of a node.
     * Example:
     *  vectorField(0, 1, 29) returns the y-component of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * dimNo   : the spatial dimension [We follow the convenction that 0:x-direction, 1:y-direction, 2:z-direction
     * nodeNo  : the current node (tag)
     */

    /* operator overloading of () */
    inline const std::valarray<lbBase_t> operator () (const int fieldNo, const int nodeNo) const // Returns pointer to beginning of a vector
    {        
        return data_[std::slice(elementSize_ * nodeNo + DXQY::nD * fieldNo,  DXQY::nD, 1)];
    }
    inline std::valarray<lbBase_t> operator () (const int fieldNo, const int nodeNo) // Returns pointer to beginning of a vector
    {
            return data_[std::slice(elementSize_ * nodeNo + DXQY::nD * fieldNo,  DXQY::nD, 1)];
    }

/*    template <typename T>
    inline void set(const int fieldNo, const int nodeNo, const T &vec)
    {
        for (auto d = 0; d < DXQY::nD; ++d)
            (*this)(fieldNo, d, nodeNo) = vec[d];
    } */
    inline std::slice_array<lbBase_t> set(const int fieldNo, const int nodeNo)
    {
        return data_[std::slice(elementSize_ * nodeNo + DXQY::nD * fieldNo,  DXQY::nD, 1)];
    }

    //   VectorFieldRet<DXQY> operator () (const int fieldNo, const int nodeNo)
 //   {
 //       VectorFieldRet<DXQY> ret(data_.data() + elementSize_ * nodeNo + DXQY::nD * fieldNo);
 //       return ret;
 //   }
    /* returns a pointer to a vector at a given node for a given field number
     * Example :
     *  vectorField(0, 29) returns a pointer to the vecotor at the 29th node for vector field 0
     *  vectorField(0, 29)[1] gives the y-component of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * nodeNo  : the current node (tag)
     */

    int getNumNodes() {return nNodes_;} // Getter for nNodes_

private:
    const int nFields_;  // Number of fields
    const int elementSize_;  // Size of a memory block
    int nNodes_;  // Number of nodes per field
    std::valarray<lbBase_t> data_;  // Pointer to the vector data
};


// Two version to handle const correctly
/*template <typename DXQY>
inline const lbBase_t& VectorField<DXQY>::operator () (const int fieldNo, const int dimNo, const int nodeNo) const
{
    return data_[elementSize_ * nodeNo + DXQY::nD * fieldNo + dimNo];
}
template <typename DXQY>
inline lbBase_t& VectorField<DXQY>::operator () (const int fieldNo, const int dimNo, const int nodeNo)
{
    return data_[elementSize_ * nodeNo + DXQY::nD * fieldNo + dimNo];
} */


// END VECTORFIELD


/*********************************************************
 * class LBFIELD: Represents a given number of lattice
 *  boltzmann distribution fields
 *
 *********************************************************/
template <typename DXQY>
class LbField
{
public:
    /* Constructor */
    LbField(const int nFields, const int nNodes):
        nFields_(nFields), elementSize_(nFields_ * DXQY::nQ), nNodes_(nNodes), data_(elementSize_ * nNodes_) {}
    /* nFields : number of vector fields
     * nNodes  : number of nodes
     */


    /* operator overloading of () */
    inline const lbBase_t& operator () (const int fieldNo, const int dirNo, const int nodeNo) const; // Returns element
    inline lbBase_t& operator () (const int fieldNo, const int dirNo, const int nodeNo); // Returns element
    /* Returns a reference to a distribution component at a node.
     * Example:
     *  lbField(0, 5, 29) returns the distribution in the 5th velocoity direction of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * dirNo   : the distribution direction (q) \in [0, ..., nQ-1]
     * nodeNo  : the current node (tag)
     */

    /* operator overloading of () */
//    inline lbBase_t* operator () (const int fieldNo, const int nodeNo) const; // Returns pointer to beginning of a vector

/*    inline const std::vector<lbBase_t> operator () (const int fieldNo, const int nodeNo) const // Returns element
    {
        auto begin_pos = data_.begin() + elementSize_ * nodeNo + DXQY::nQ * fieldNo;
        return std::vector<lbBase_t> (begin_pos, begin_pos + DXQY::nQ);
    }
    inline std::vector<lbBase_t> operator () (const int fieldNo, const int nodeNo) // Returns element
    {
        auto begin_pos = data_.begin() + elementSize_ * nodeNo + DXQY::nQ * fieldNo;
        return std::vector<lbBase_t> (begin_pos, begin_pos + DXQY::nQ);
    } */
    inline const std::valarray<lbBase_t> operator () (const int fieldNo, const int nodeNo) const // Returns element
    {
        return data_[std::slice(elementSize_ * nodeNo + DXQY::nQ * fieldNo, DXQY::nQ, 1)];
    }
    inline std::valarray<lbBase_t> operator () (const int fieldNo, const int nodeNo) // Returns element
    {
        return data_[std::slice(elementSize_ * nodeNo + DXQY::nQ * fieldNo, DXQY::nQ, 1)];

    }

    /* Returns a pointer to a lb distribution at a given node for a given field number
     * Example:
     *  lbField(0, 29) returns a pointer to the lb distribution at the 29th node for vector field 0.
     *  vectorField(0, 29)[5] returns 5th velocoity direction distribution of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * nodeNo  : the current node (tag)
     */

    inline void swapData(LbField& field);
    /* swapData swaps the data_ pointer for the this object and the object, 'field', given
     *  as input. This is done straight after propagation.
     *
     * field :  the temporary field where the propagated distribution values are stored
     */

    int getNumNodes() {return nNodes_;} // Getter for nNodes_

private:
    const int nFields_;  // Number of fields
    const int elementSize_;  // Size of a memory block
    int nNodes_;  // number of nodes per field
    std::valarray<lbBase_t> data_;  // Container for the field
};


template <typename DXQY>
inline const lbBase_t& LbField<DXQY>::operator () (const int fieldNo, const int dirNo, const int nodeNo) const // Returns element
{
    return data_[elementSize_ * nodeNo + DXQY::nQ * fieldNo + dirNo];
}
template <typename DXQY>
inline lbBase_t& LbField<DXQY>::operator () (const int fieldNo, const int dirNo, const int nodeNo) // Returns element
{
    return data_[elementSize_ * nodeNo + DXQY::nQ * fieldNo + dirNo];
}


template <typename DXQY>
inline void LbField<DXQY>::swapData(LbField<DXQY>& field)
{
    data_.swap(field.data_);
}

// END LBFIELD


#endif // LBFIELD_H
