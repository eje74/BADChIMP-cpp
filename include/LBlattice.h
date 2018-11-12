#ifndef LBLATTICE_H
#define LBLATTICE_H

/*! LATTICE class
 *  Contains all information about the basis vectors, and operators for vectors and tensors
 *  like scalarproducs.
 */

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

    // The inner product of a cartesian vector and a basis vector
    // Used for inner products with Cartesian vecotrs (\sum_i c_{\alpha i} u_i)
    double innerProductDMajor(const int qDirection, const double* vec) const;

    // The projection of lattice botlzmann field and a basis vector
    // Used to calculate first moments (\sum_\alpha c_{\alpha i} f_\alpha)
    double innerProductQMajor(const int dimension, const double* vec) const;

    // Standard scalar products for Cartesian vectors
    double innerProduct(const double *leftVec, const double *rightVec) const;

    // Different powers of the sound speed.
    const double c2Inv_ = 3.0;
    const double c2_ = 1.0 / c2Inv_;
    const double c4Inv_ = 9.0;
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


#endif // LBLATTICE_H
