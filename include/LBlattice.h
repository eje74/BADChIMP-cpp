#ifndef LBLATTICE_H
#define LBLATTICE_H

/*! LATTICE class
 *  Contains all information about the basis vectors, and operators for vectors and tensors
 *  like scalarproducs.
 */


class Lattice // Can we make static class types to increase speed. And add it with templates
{
public:
    // Constructure and destructor
    Lattice(const int nDirections, const int nDimensions);
    ~Lattice();

    // Getter for directions and dimensions
    int nQ() const;
    int nD() const;

    // Getter and setter function for lattice weights
    double w(const int qDirection) const;
    void setW(const int qDirection, const double value);

    // Getter function for basis vector components
    int c(const int qDirection, const int dimension) const;
    void setC(const int qDirection, const int dimension, const int value);

    // Returns the reverse direction of a lattice directions
    int reverseDirection(const int qDirection) const;

    // Standard scalar products for Cartesian vectors
    // dot
    template <typename BASETYPE>
    BASETYPE dot(const BASETYPE* leftVec, const BASETYPE* rightVec) const;
    // The inner product of a cartesian vector and a basis vector
    // Used for inner products with Cartesian vecotrs (\sum_i c_{\alpha i} u_i)
    // cDot
    template <typename BASETYPE>
    BASETYPE cDot(const int qDir, const BASETYPE* rightVec) const;

    // The projection of lattice botlzmann field and a basis vector
    // Used to calculate first moments (\sum_\alpha c_{\alpha i} f_\alpha)
    // qSum
    // qSumC
    // qSumCC
    template <typename BASETYPE>
    BASETYPE qSum(const BASETYPE* dist) const;
    template <typename BASETYPE>
    BASETYPE qSumC(const int dim, const BASETYPE* dist) const;
    // Different powers of the sound speed.
    const double c2Inv_ = 3.0;
    const double c2_ = 1.0 / c2Inv_;
    const double c4Inv_ = 9.0;
    const double c4Inv0_5_ = 4.5;
    const double c4_ = 1.0 / c4Inv_;

private:
    const int nDirections_;  // Size of velocity basis
    const int nDimensions_;  // Number of spacial dimensions
    double *w_; // Lattice weights
    int *cDMajor_;  // Velocity basis on the form c_\alpha,i = c[nD*alpha + i]
    int *cQMajor_;  // Velocity basis on the form c_\alpha,i = c[nQ*i + alpha]
    int reverseStep_;  // Used to reverse direction (assuming that \hat{\alpha} = \alpha + reverseStep_
    int nNonZeroDirections_;  // = nQ_ if no reste particle, = (nQ_ - 1) if rest particle
};

// Getter for directions and dimensions
inline int Lattice::nQ() const// Getter for the number of lattice directions
{
    return nDirections_;
}

inline int Lattice::nD() const  // Getter for the number of spatial dimensions
{
    return nDimensions_;
}


// Getter and setter function for lattice weights
inline double Lattice::w(const int qDirection) const
{
    return w_[qDirection];
}

inline void Lattice::setW(const int qDirection, const double value)
{
    w_[qDirection] = value;
}


// Getter function for basis vector components
inline int Lattice::c(const int qDirection, const int dimension) const
{
    return cDMajor_[nDimensions_*qDirection + dimension];
}

inline void Lattice::setC(const int qDirection, const int dimension, const int value)
{
    cDMajor_[qDirection * nDimensions_ + dimension] = value;
    cQMajor_[dimension * nDirections_ + qDirection] = value;
}


// Returns the reverse direction of a lattice directions
inline int Lattice::reverseDirection(const int qDirection) const // Returns the reverse direction
{
    /* Here we have assumed that the reverse directions are given as:
     *     alpha_reverse = alpha + <number of non zero directons>/2 ,
     * and that the rest particle velocity is that last element in the
     * basis vector.
     */
    return (qDirection + reverseStep_) % nNonZeroDirections_;
}


// The inner product of a cartesian vector and a basis vector
template <typename BASETYPE>
inline BASETYPE Lattice::dot(const BASETYPE* leftVec, const BASETYPE* rightVec) const
{
    BASETYPE ret = 0;
    for (int d = 0; d < nDimensions_; d++) {
        ret += leftVec[d]*rightVec[d];
    }
    return ret;
}

// The projection of lattice botlzmann field and a basis vector
template <typename BASETYPE>
inline BASETYPE Lattice::cDot(const int qDir, const BASETYPE* vec) const
{
    BASETYPE ret = 0;
    for (int d = 0; d < nDimensions_; d++) {
        ret += cDMajor_[qDir * nDimensions_ + d] * vec[d];
    }
    return ret;
}



template <typename BASETYPE>
inline BASETYPE Lattice::qSum(const BASETYPE* dist) const
{
    double ret = 0;
    for (int q = 0; q < nDirections_; ++q)
        ret += dist[q];
    return ret;
}

template <typename BASETYPE>
inline BASETYPE Lattice::qSumC(const int dim, const BASETYPE* dist) const
{
    double ret = 0;
    for (int q = 0; q < nDirections_; ++q)
        ret +=  cQMajor_[dim * nDirections_ + q] * dist[q];
    return ret;
}



#endif // LBLATTICE_H

