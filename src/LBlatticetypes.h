#ifndef LBLATTICETYPES_H
#define LBLATTICETYPES_H

#include "LBglobal.h"
#include "LBd2q9.h"
#include "LBd3q19.h"

/**************************************************************
 * The different lattice types are defined in separate head files
 * with the naming convention 'LBd<num dim>q<num dir>.h'.
 *
 * Notes:
 *  In the current definitions we have assumend that the last
 *   position is reserved for the "rest-particle direction"
 *
 * Example format for a lattice structure, using LBd2q9.h:
 *
 * ******* CONSTANTS AND ARRAYS ********
 *
 * BASIC LATTICE STRUCTURE
static constexpr int nD = 2;         // Number of spatial dimensions
static constexpr int nQ = 9;         // Number of lattice directions
static constexpr int nDirPairs_ = 4; // Number of unsigned bounce back directions
static constexpr int nQNonZero_ = 8; // Number of non zeros basis directions

 * LATTICE SOUND VELOCITY ( = sqrt(C_2) )
 * 2nd order \sum_\alpha w_\alpha c_{\alpha_i}c_{\alpha_j} = C_2 \delta_{ij}
 * 4th order \sum_\alpha w_\alpha c_{\alpha_i}c_{\alpha_j}c_{\alpha_k}c_{\alpha_l} = C_4(\delta_{ij}\delta_{kl} + \delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})
static constexpr lbBase_t c2Inv = 3.0;             // = 1/C_2
static constexpr lbBase_t c2 = 1.0 / c2Inv;        // = C_2
static constexpr lbBase_t c4Inv = 9.0;             // = 1/C_4
static constexpr lbBase_t c4 = 1.0 / c4Inv;        // = C_4
static constexpr lbBase_t c4Inv0_5 = 0.5 * c4Inv;  // = 1/(2*C_4)

 * LATTICE WEIGHTS
static constexpr lbBase_t w0 = 4.0/9.0;  // Rest particle
static constexpr lbBase_t w1 = 1.0/9.0;  // Nearest neighbor
static constexpr lbBase_t w2 = 1.0/36.0; // Next nearest neighbor
 * Derived qunatetatis
static constexpr lbBase_t w0c2Inv = w0*c2Inv;
static constexpr lbBase_t w1c2Inv = w1*c2Inv;
static constexpr lbBase_t w2c2Inv = w2*c2Inv;
 * Weights for each lattice direction
static constexpr lbBase_t w[9] = {w1, w2, w1, w2, w1, w2, w1, w2, w0};

 * LATTICE BASIS VECTORS
 * Basis vector for each lattice direction on the form [vx[0], vy[0], vx[1], vy[1], ...]
static constexpr int cDMajor_[18] = {1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 0, 0};
 * Norm of each lattice vector
static constexpr lbBase_t cNorm[9] = {1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 0.0};

 * COLOR GRADIENT CONSTATNS
 * Weights used in the surface tension perturbation. Ref Reis, Phillips (2007)
static constexpr lbBase_t B0 = -4.0/27.0;  // Rest particle
static constexpr lbBase_t B1 = 2.0/27.0;   // Nearest Neighbor
static constexpr lbBase_t B2 = 5.0/108.0;  // Next nearest neighbor
 * Array of weights for each lattice direction.
 * NB. Remember to define the static arrays in the LBd2q9.cpp file as well.
static constexpr lbBase_t B[9] = {B1, B2, B1, B2, B1, B2, B1, B2, B0};


 * ******* FUNCTIONS ********
 * LATTICE BASIS VECTORS
inline static int c(const int qDirection, const int dimension)
 * returns a component of a basis vector.
 *  qDirection : lattice direction for the basis vector
 *  dimension  : component of the basis vector (0:x-component, 1:y-component, ...)

 * BASIS LATTICE STRUCTURE
inline static int reverseDirection(const int qDirection) {return (qDirection + nDirPairs_) % nQNonZero_;}
 * returns the reverse lattice direction of the input direction
 *  qDirection : lattice direction

 * VECTOR OPERATORS
 * Vector product
static lbBase_t dot(const lbBase_t* leftVec, const lbBase_t* rightVec);
 * returns the vector product of the arrays
 *  leftVec  : left hand vector
 *  rightVec : right hand vector
 *
template<typename T>
static T cDot(const int qDir, const T* rightVec);
* returns the vector product of a basis vector and an array
*  qDir : lattice direction
*  rightVec : vector to be doted with the basis vector
*
static void cDotAll(const lbBase_t* vec, lbBase_t* ret);
* calculates a list of the dot procuect of a vector with all basis vectors
*  vec : the vector to be used in the dot product
*  ret : the list of dot products to be returned [vec*c[0], vec*c[1], ...]
*
static void grad(const lbBase_t* rho, lbBase_t* ret);
* calculates the gradient of rho
*  rho : list of scalar values for each lattice direction
*  ret : the gradient to be returned.
*
static void qSum(const lbBase_t* dist, lbBase_t& ret);
* calculates the zeroth moment of distribution
*  dist: list of the elements of the distribution for each lattice direction
*  ret : the zeroth moment, to be returned.
*
static void qSumC(const lbBase_t* dist, lbBase_t* ret);
* calculates the first moment of distribution
*  dist: list of the elements of the distribution for each lattice direction
*  ret : the first moment, to be returned.
*
// Two phase
static void gradPush(const lbBase_t& scalarVal, const int* neighList, VectorField<D2Q9>& grad);

};

*/


#endif // LBLATTICETYPES_H
