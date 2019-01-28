#ifndef LBFIELD_H
#define LBFIELD_H

#include "LBglobal.h"

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
    ScalarField(const int nFields, const int nNodes);
    /* Destructor
     */
    ~ScalarField();

    /* operator overloading of (). */
    lbBase_t& operator () (const int fieldNo, const int nodeNo) const;
    /*
     * fieldNo : the current field
     * nodeNo : the current node (tag)
     *
     * Returns the scalar value for the given field at the given node
     * Usage:
     * rho(0, 32) = 3.5 // sets the value of field_0's node 32 to 3.5
     */

private:
    const int nFields_;  // Number of fields
    int nNodes_;  // Number of nodes in each field
    lbBase_t* data_;  // Pointer to the scalar data
};


inline lbBase_t& ScalarField::operator () (const int fieldNo, const int nodeNo) const
{
    return data_[nFields_ * nodeNo + fieldNo];
}
// END SCALARFIELD



/*********************************************************
 * class VECTORFIELD: Represents a given number of vector
 *  fields
 *
 *********************************************************/
class VectorField
{
public:
    /* Constructor
     * nFields     : number of vector fields
     * nDimensions : number of spatial dimensions
     * nNodes      : number of nodes
     */
    VectorField(const int nFields, const int nDim, const int nNodes);

    /* Destructor
     */
    ~VectorField();

    /* operator overloading of () */
    lbBase_t& operator () (const int fieldNo, const int dimNo, const int nodeNo) const; // Returns element
    /* Returns a reference to a vector component of a node.
     * Example:
     *  vectorField(0, 1, 29) returns the y-component of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * dimNo   : the spatial dimension [We follow the convenction that 0:x-direction, 1:y-direction, 2:z-direction
     * nodeNo  : the current node (tag)
     */

    /* operator overloading of () */
    lbBase_t* operator () (const int fieldNo, const int nodeNo) const; // Returns pointer to beginning of a vector
    /* returns a pointer to a vector at a given node for a given field number
     * Example :
     *  vectorField(0, 29) returns a pointer to the vecotor at the 29th node for vector field 0
     *  vectorField(0, 29)[1] gives the y-component of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * nodeNo  : the current node (tag)
     */

private:
    const int nFields_;  // Number of fields
    const int nDim_;  // Numper of spatial dimensions
    const int elementSize_;  // Size of a memory block
    int nNodes_;  // Number of nodes per field
    lbBase_t* data_;  // Pointer to the vector data
};


inline lbBase_t& VectorField::operator () (const int fieldNo, const int dimNo, const int nodeNo) const
{
    return data_[elementSize_ * nodeNo + nDim_ * fieldNo + dimNo];
}


inline lbBase_t* VectorField::operator () (const int fieldNo, const int nodeNo) const
{
    return &data_[elementSize_ * nodeNo + nDim_ * fieldNo];
}
// END VECTORFIELD


/*********************************************************
 * class LBFIELD: Represents a given number of lattice
 *  boltzmann distribution fields
 *
 *********************************************************/
class LbField
{
public:
    /* Constructor */
    LbField(const int nFields, const int nDir, const int nNodes);
    /* nFields : number of vector fields
     * nDir    : number of lattice directions (q)
     * nNodes  : number of nodes
     */

    /* Destructor
     */
    ~LbField();

    /* operator overloading of () */
    lbBase_t& operator () (const int fieldNo, const int dirNo, const int nodeNo) const; // Returns element
    /* Returns a reference to a distribution component at a node.
     * Example:
     *  lbField(0, 5, 29) returns the distribution in the 5th velocoity direction of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * dirNo   : the distribution direction (q) \in [0, ..., nQ-1]
     * nodeNo  : the current node (tag)
     */

    /* operator overloading of () */
    lbBase_t* operator () (const int fieldNo, const int nodeNo) const; // Returns pointer to beginning of a vector
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

private:
    const int nFields_;  // Number of fields
    const int nDir_;  // Number of velocies in the lb distribtuion
    const int elementSize_;  // Size of a memory block
    int nNodes_;  // number of nodes per field
    lbBase_t* data_;  // Pointer to the distribution data
};



inline lbBase_t& LbField::operator () (const int fieldNo, const int dirNo, const int nodeNo) const // Returns element
{
    return data_[elementSize_ * nodeNo + nDir_ * fieldNo + dirNo];
}


inline lbBase_t* LbField::operator () (const int fieldNo, const int nodeNo) const // Returns pointer to beginning of a vector
{
    return &data_[elementSize_ * nodeNo + nDir_ * fieldNo];
}


inline void LbField::swapData(LbField& field)
{
    lbBase_t* dataTmp;
    dataTmp = field.data_;
    field.data_ = data_;
    data_ = dataTmp;
}

// END LBFIELD


#endif // LBFIELD_H
