#ifndef LBFIELD_H
#define LBFIELD_H

#include "LBglobal.h"
#include "LBgrid.h"
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
    // return reference to data_ vector
    const std::valarray<lbBase_t>& get_data() const {return data_;}

    // use operator() to calculate the index of a given field in the data_ vector
    std::vector<int> get_field_index(int fieldNo, std::vector<int>& nodes) const {
      std::vector<int> ind(nodes.size());
      for (auto n=0; n<nodes.size(); ++n) {
        ind[n] = &(*this)(fieldNo, nodes[n]) - &data_[0];
      }
      return ind;
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
    std::valarray<lbBase_t> data_; // Container for scalar data
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
    /* returns the local vector at a given node for a given field number, as a valarray
     * Example :
     *  vectorField(0, 29) returns the vector at the 29th node for vector field 0
     *  vectorField(0, 29)[1] gives the y-component of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * nodeNo  : the current node (tag)
     */
    inline const std::valarray<lbBase_t> operator () (const int fieldNo, const int nodeNo) const // Returns pointer to beginning of a vector
    {        
        return data_[std::slice(elementSize_ * nodeNo + DXQY::nD * fieldNo,  DXQY::nD, 1)];
    }
    inline std::valarray<lbBase_t> operator () (const int fieldNo, const int nodeNo) // Returns pointer to beginning of a vector
    {
            return data_[std::slice(elementSize_ * nodeNo + DXQY::nD * fieldNo,  DXQY::nD, 1)];
    }

    inline std::slice_array<lbBase_t> set(const int fieldNo, const int nodeNo)
    {
        return data_[std::slice(elementSize_ * nodeNo + DXQY::nD * fieldNo,  DXQY::nD, 1)];
    }

    int getNumNodes() {return nNodes_;} // Getter for nNodes_
    int getNumNodes() const {return nNodes_;}
    int size() {return nNodes_;} // laternative Getter for nNodes_    
    int size() const {return nNodes_;} // laternative Getter for nNodes_    
    int num_fields() const {return nFields_;} // Getter for nFields_
   //JLV
    const std::valarray<lbBase_t>& get_data() {return data_;}

    // use operator() to calculate the index of a given field in the data_ vector
    std::vector<int> get_field_index(int fieldNo, std::vector<int>& nodes) {
      std::vector<int> ind; ind.reserve(nodes.size()*DXQY::nD);
      for (auto n=0; n<nodes.size(); ++n) {
        for (auto d=0; d<DXQY::nD; ++d) {
          ind.push_back( &(*this)(fieldNo, d, nodes[n]) - &data_[0] );
        }
      }
      return ind;
    }
    //JLV


private:
    const int nFields_;  // Number of fields
    const int elementSize_;  // Size of a memory block
    int nNodes_;  // Number of nodes per field
    std::valarray<lbBase_t> data_;  // Pointer to the vector data
};


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
    inline const lbBase_t& operator () (const int fieldNo, const int dirNo, const int nodeNo) const // Returns element
    {
        return data_[elementSize_ * nodeNo + DXQY::nQ * fieldNo + dirNo];
    }
    inline lbBase_t& operator () (const int fieldNo, const int dirNo, const int nodeNo) // Returns element
    {
        return data_[elementSize_ * nodeNo + DXQY::nQ * fieldNo + dirNo];
    }
    /* Returns a reference to a distribution component at a node.
     * Example:
     *  lbField(0, 5, 29) returns the distribution in the 5th velocoity direction of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * dirNo   : the distribution direction (q) \in [0, ..., nQ-1]
     * nodeNo  : the current node (tag)
     */

    /* operator overloading of () */
    /* returns a valarray for the distribution at a given node for a given field number
     * Example :
     *  LbField(0, 29) returns the distributions for the 29th node for field 0
     *  LbField(0, 29)[1] gives the 1-velocity distribution of field 0's 29th node.
     *
     * fieldNo : the current field [remember 0 if only one field]
     * nodeNo  : the current node (tag)
     */
    inline const std::valarray<lbBase_t> operator () (const int fieldNo, const int nodeNo) const // Returns element
    {
        return data_[std::slice(elementSize_ * nodeNo + DXQY::nQ * fieldNo, DXQY::nQ, 1)];
    }
    inline std::valarray<lbBase_t> operator () (const int fieldNo, const int nodeNo) // Returns element
    {
        return data_[std::slice(elementSize_ * nodeNo + DXQY::nQ * fieldNo, DXQY::nQ, 1)];

    }

    inline std::slice_array<lbBase_t> set(const int fieldNo, const int nodeNo)
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

    /* Propagate values from the node with node number nodeNo and values f_omega to the neighboring distribution on
     * int this field.
     * */
    template <typename T> 
    inline void propagateTo(const int & fieldNo, const int & nodeNo, const T& f_omega, const Grid<DXQY> &grid)
    {
        int ind_shift = fieldNo * DXQY::nQ;
        for (int q = 0; q < DXQY::nQ; ++q) {
            data_[elementSize_ * grid.neighbor(q, nodeNo) + ind_shift + q] = f_omega[q];
        }
    }

    inline void swapData(LbField& field) { data_.swap(field.data_); }
    /* swapData swaps the data_ pointer for the this object and the object, 'field', given
     *  as input. This is done straight after propagation.
     *
     * field :  the temporary field where the propagated distribution values are stored
     */

    int getNumNodes() const {return nNodes_;} // Getter for nNodes_
    int num_fields() const {return nFields_;} // Getter for nFields_

private:
    const int nFields_;  // Number of fields
    const int elementSize_;  // Size of a memory block
    int nNodes_;  // number of nodes per field
    std::valarray<lbBase_t> data_;  // Container for the field
};

// END LBFIELD


#endif // LBFIELD_H
