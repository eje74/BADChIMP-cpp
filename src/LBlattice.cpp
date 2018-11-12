#include "LBlattice.h"

Lattice::Lattice(const int nDirections, const int nDimensions) : nDirections_(nDirections), nDimensions_(nDimensions)  // Lattice constructor
{
    w_ = new double [nDirections_];                     /* List of lattice weights */
    cDMajor_ = new int [nDirections_ * nDimensions_];   /*  Basis velocities on the form c[alpha][dim]
                                                         *  So that c = [c(alpha=0)_x|c(alpha=0)_y, ...]
                                                         */
    cQMajor_ = new int [nDirections_ * nDimensions_];   /* Basis velocities on the form c[dim][alpha].
                                                         * So that c = [c(alpha=0)_x|c(alpha=1)_x, ...]
                                                         */
    reverseStep_ = nDirections_ / 2;                    /* Number of bounce back pairs */
    nNonZeroDirections_ = 2 * reverseStep_;             /* Number of non-rest particles velocities */
}

Lattice::~Lattice()  // Lattice destructor
{
    delete [] w_;
    delete [] cDMajor_;
    delete [] cQMajor_;
//        std::cout << "Destructed a Lattice object" << std::endl;
}


// Getter for directions and dimensions
int Lattice::nQ() const// Getter for the number of lattice directions
{
    return nDirections_;
}

int Lattice::nD() const  // Getter for the number of spatial dimensions
{
    return nDimensions_;
}


// Getter and setter function for lattice weights
double Lattice::w(const int qDirection) const
{
    return w_[qDirection];
}

void Lattice::setW(const int qDirection, const double value)
{
    w_[qDirection] = value;
}


// Getter function for basis vector components
int Lattice::c(const int qDirection, const int dimension) const
{
    return cDMajor_[nDimensions_*qDirection + dimension];
}

void Lattice::setC(const int qDirection, const int dimension, const int value)
{
    cDMajor_[qDirection * nDimensions_ + dimension] = value;
    cQMajor_[dimension * nDirections_ + qDirection] = value;
}


// Returns the reverse direction of a lattice directions
int Lattice::reverseDirection(const int qDirection) const // Returns the reverse direction
{
    /* Here we have assumed that the reverse directions are given as:
     *     alpha_reverse = alpha + <number of non zero directons>/2 ,
     * and that the rest particle velocity is that last element in the
     * basis vector.
     */
    return (qDirection + reverseStep_) % nNonZeroDirections_;
}


// The inner product of a cartesian vector and a basis vector
double Lattice::innerProductDMajor(const int qDirection, const double* vec) const // Used for inner products with Cartesian vecotrs (\sum_i c_{\alpha i} u_i)
{
    double ret = 0;
    for (int d = 0; d < nDimensions_; d++) {
        ret += vec[d]*cDMajor_[qDirection*nDimensions_ + d];
    }
    return ret;
}


// The projection of lattice botlzmann field and a basis vector
double Lattice::innerProductQMajor(const int dimension, const double* vec) const // Used to calculate first moments (\sum_\alpha c_{\alpha i} f_\alpha)
{
    double ret = 0;
    for (int q = 0; q < nDirections_; q++) {
        ret += vec[q]*cQMajor_[dimension*nDirections_ + q];
    }
    return ret;
}


// Standard scalar products for Cartesian vectors
double Lattice::innerProduct(const double *leftVec, const double *rightVec) const
{
    double ret = 0;
    for (int d = 0; d < nDimensions_; d++) {
        ret += leftVec[d]*rightVec[d];
    }
    return ret;
}
