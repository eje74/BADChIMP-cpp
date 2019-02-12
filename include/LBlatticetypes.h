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
 * Example format for a lattice structure:
 *

struct D2Q9 {  // Structure naming convention 'D<num dim>Q<num dir>'

static constexpr int nD = 2;  // Number of spatial dimensions
static constexpr int nQ = 9;  // Number of lattice directions
static constexpr int nDirPairs_ = 4; // Number of unsigned bounce back directions
static constexpr int nQNonZero_ = 8; // Number of non zeros basis directions

// * Constants assosiated with the lattice symmetries
// * 2nd order \sum_\alpha w_\alpha c_{\alpha_i}c_{\alpha_j} = C_2 \delta_{ij}
// * 4th order \sum_\alpha w_\alpha c_{\alpha_i}c_{\alpha_j}c_{\alpha_k}c_{\alpha_l} = C_4(\delta_{ij}\delta_{kl} + \delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})
static constexpr lbBase_t c2Inv = 3.0;  // = 1/C_2
static constexpr lbBase_t c2 = 1.0 / c2Inv;  // = C_2
static constexpr lbBase_t c4Inv = 9.0;  // = 1/C_4
static constexpr lbBase_t c4 = 1.0 / c4Inv; // = C_4
static constexpr lbBase_t c4Inv0_5 = 0.5 * c4Inv;  // = 1/(2*C_4)

// * list of lattice weights
static constexpr lbBase_t w[9] = {1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 4.0/9.0};

// * list of lattice vectors given on the form [c_{\alpha x}, c_{\alpha y}, ...]
static constexpr int cDMajor_[18] = {1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 0, 0};

// * The function c(\alpha, i) returns c_{\alpha i}
inline static int c(const int qDirection, const int dimension)  {return cDMajor_[nD*qDirection + dimension];}

// * The function reverseDirection(\alpha) returns the reverse direction of alpha so that
// *         c_{reverseDirection(\alpha), i} + c_{\alpha, i} = 0
inline static int reverseDirection(const int qDirection) {return (qDirection + nDirPairs_) % nQNonZero_;}

// * Function dot(vec1, vec2) is the standard scalar product, that is:
// *  dot(vec1, vec2) = \sum_{i=0}^{nD-1} vec1_i*vec2_i
static lbBase_t dot(const lbBase_t* leftVec, const lbBase_t* rightVec);

// * Function cDot(\alpha, vec) is the scalar product of lattice direction \alpha and vec, that is
// *   cDot(\alpha, vec) = \sum_{i=0}^{nD-1} c_{\alpha, i}*vec_i
static lbBase_t cDot(const int qDir, const lbBase_t* rightVec);

// * Function cDotAll(vec, ret) returns a list of vec's scalar product with all lattice direction
// * so that ret = {cDot(0, vec), cDot(1, vec), ..., cDot(nQ-1, vec)}, that is
// *   ret[\alpha] = \sum_{i=0}^{nD-1} c_{\alpha, i}*vec_i
static void cDotAll(const lbBase_t* vec, lbBase_t* ret);

// * Function qSum(f, ret) returns the zeroth moment of f, that is
// *  ret = \sum_\alpha f_\alpha
static void qSum(const lbBase_t* dist, lbBase_t& ret);

// * Function qSumC(f, ret) returns the first moment of f, that is
// *  ret[i] = \sum_\alpha c_{\alpha i}f_\alpha
static void qSumC(const lbBase_t* dist, lbBase_t* ret);


// * Function gradPush
static void gradPush(const lbBase_t& rho, const int* neighList, VectorField<D2Q9>& grad);

};


inline lbBase_t D2Q9::dot(const lbBase_t* leftVec, const lbBase_t* rightVec)
{
    return leftVec[0]*rightVec[0] + leftVec[1]*rightVec[1];
}

inline lbBase_t D2Q9::cDot(const int qDir, const lbBase_t* rightVec)
{
    return c(qDir, 0)*rightVec[0] + c(qDir, 1)*rightVec[1];
}

inline void D2Q9::cDotAll(const lbBase_t* vec, lbBase_t* ret)
{
    ret[0] = vec[0];
    ret[1] = vec[0] + vec[1];
    ret[2] = vec[1];
    ret[3] = vec[1] - vec[0];
    ret[4] = -vec[0];
    ret[5] = -vec[0] - vec[1];
    ret[6] = -vec[1];
    ret[7] = vec[0] - vec[1];
    ret[8] = 0;
}

inline void D2Q9::qSum(const lbBase_t* dist, lbBase_t& ret)
{
    ret = 0.0;
    for (int q = 0; q < nQ; ++q)
        ret += dist[q];
}

inline void D2Q9::qSumC(const lbBase_t* dist, lbBase_t* ret)
{
    ret[0] = dist[0] + dist[1]            - dist[3] - dist[4]  - dist[5]           + dist[7];
    ret[1] =           dist[1] + dist[2]  + dist[3]            - dist[5] - dist[6] - dist[7];
}
*/

#endif // LBLATTICETYPES_H
