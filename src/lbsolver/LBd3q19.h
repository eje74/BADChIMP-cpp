#ifndef LBD3Q19_H
#define LBD3Q19_H

#include "LBglobal.h"
#include <vector>

// See "LBlatticetypes.h" for description of the structure

// TO MAKE CHANGES TO THIS FILE, MAKE THE CHANGES IN "PythonScripts/writeLatticeFile.py",
// RUN THE SCRIPT AND PLACE RESULTING FILES IN "src/"

struct D3Q19{

static constexpr int nD = 3;
static constexpr int nQ = 19;
static constexpr int nDirPairs_ = 9;
static constexpr int nQNonZero_ = 18;

static constexpr lbBase_t c2Inv = 3.0;
static constexpr lbBase_t c4Inv = 9.0;
static constexpr lbBase_t c2 = 1.0 / c2Inv;
static constexpr lbBase_t c4 = 1.0 / c4Inv;
static constexpr lbBase_t c4Inv0_5 = 0.5 * c4Inv;

static constexpr lbBase_t w0 = 12.0/36.0;
static constexpr lbBase_t w0c2Inv = w0*c2Inv;
static constexpr lbBase_t w1 = 2.0/36.0;
static constexpr lbBase_t w1c2Inv = w1*c2Inv;
static constexpr lbBase_t w2 = 1.0/36.0;
static constexpr lbBase_t w2c2Inv = w2*c2Inv;

static constexpr lbBase_t w[19] = {w1, w1, w1, w2, w2, w2, w2, w2, w2, w1, w1, w1, w2, w2, w2, w2, w2, w2, w0};
static constexpr int cDMajor_[57] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, -1, 0, 1, 0, 1, 1, 0, -1, 0, 1, 1, 0, 1, -1, -1, 0, 0, 0, -1, 0, 0, 0, -1, -1, -1, 0, -1, 1, 0, -1, 0, -1, -1, 0, 1, 0, -1, -1, 0, -1, 1, 0, 0, 0};
static constexpr lbBase_t cNorm[19] = {1.0, 1.0, 1.0, SQRT2, SQRT2, SQRT2, SQRT2, SQRT2, SQRT2, 1.0, 1.0, 1.0, SQRT2, SQRT2, SQRT2, SQRT2, SQRT2, SQRT2, 0.0};
static constexpr lbBase_t cNormInv[19] = {1.0, 1.0, 1.0, SQRT2INV, SQRT2INV, SQRT2INV, SQRT2INV, SQRT2INV, SQRT2INV, 1.0, 1.0, 1.0, SQRT2INV, SQRT2INV, SQRT2INV, SQRT2INV, SQRT2INV, SQRT2INV, 0.0};
static constexpr int reverseDirection_[19] = {9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 18};
static constexpr lbBase_t B0 = -12.0/54.0;
static constexpr lbBase_t B1 = 1.0/54.0;
static constexpr lbBase_t B2 = 2.0/54.0;
static constexpr lbBase_t B[19] = {B1, B1, B1, B2, B2, B2, B2, B2, B2, B1, B1, B1, B2, B2, B2, B2, B2, B2, B0};

static constexpr lbBase_t UnitMatrixLowTri[6] = {1, 0, 1, 0, 0, 1};

// Functions

inline static int c(const int qDirection, const int dimension)  {return cDMajor_[nD*qDirection + dimension];}
inline static int reverseDirection(const int qDirection) {return reverseDirection_[qDirection];}

inline static std::vector<int> c(const int qDirection) {
std::vector<int> cq(cDMajor_ + nD*qDirection, cDMajor_ + nD*qDirection + nD);
return cq;
}

inline static std::valarray<lbBase_t> cValarray(const int qDirection) {
std::valarray<lbBase_t> cq(nD);
const int dind = nD*qDirection;
for (int d=0; d<nD; ++d)
cq[d] = cDMajor_[dind + d];
return cq;
}

template <typename T1, typename T2>
inline static lbBase_t dot(const T1 &leftVec, const T2 &rightVec);
template<typename T>
inline static T cDot(const int qDir, const T* rightVec);
template<typename T>
inline static lbBase_t cDotRef(const int qDir, const T& rightVec);
template <typename T>
inline static std::valarray<lbBase_t> cDotAll(const T &vec);
template <typename T>
inline static std::valarray<lbBase_t> grad(const T &rho);

template <typename T>
inline static lbBase_t divGrad(const T &rho);

template <typename T>
inline static lbBase_t qSum(const T &dist);
template <typename T>
inline static std::valarray<lbBase_t> qSumC(const T &dist);

inline static int c2q(const std::vector<int> &v);

template <typename T>
inline static std::valarray<lbBase_t> qSumCCLowTri(const T &dist);

template <typename T>
inline static lbBase_t traceLowTri(const T &lowTri);

template <typename T>
inline static lbBase_t traceOfMatrix(const T &mat);

inline static std::valarray<lbBase_t> deltaLowTri();

inline static std::valarray<lbBase_t> deltaMatrix();

template <typename T1, typename T2>
inline static lbBase_t contractionLowTri(const T1 &lowTri1, const T2 &lowTri2);

template <typename T>
inline static lbBase_t contractionRank2(const T &mat1, const T &mat2);

template <typename T>
inline static std::valarray<lbBase_t> matrixMultiplication(const T &mat1, const T &mat2);

template <typename T1, typename T2>
inline static std::valarray<lbBase_t> contractionLowTriVec(const T1 &lowTri, const T2 &vec);

};


template <typename T1, typename T2>
inline lbBase_t D3Q19::dot(const T1 &leftVec, const T2 &rightVec)
{
    return leftVec[0]*rightVec[0] + leftVec[1]*rightVec[1] + leftVec[2]*rightVec[2];
}

template<typename T>
inline T D3Q19::cDot(const int qDir, const T* rightVec)
{
    return c(qDir, 0)*rightVec[0] + c(qDir, 1)*rightVec[1] + c(qDir, 2)*rightVec[2];
}

template<typename T>
inline lbBase_t D3Q19::cDotRef(const int qDir, const T& rightVec)
{
    return c(qDir, 0)*rightVec[0] + c(qDir, 1)*rightVec[1] + c(qDir, 2)*rightVec[2];
}

template <typename T>
inline std::valarray<lbBase_t> D3Q19::cDotAll(const T &vec)
{
std::valarray<lbBase_t> ret(nQ);
ret[0] = +vec[0];
ret[1] = +vec[1];
ret[2] = +vec[2];
ret[3] = +vec[0] +vec[1];
ret[4] = +vec[0] -vec[1];
ret[5] = +vec[0] +vec[2];
ret[6] = +vec[0] -vec[2];
ret[7] = +vec[1] +vec[2];
ret[8] = +vec[1] -vec[2];
ret[9] = -vec[0];
ret[10] = -vec[1];
ret[11] = -vec[2];
ret[12] = -vec[0] -vec[1];
ret[13] = -vec[0] +vec[1];
ret[14] = -vec[0] -vec[2];
ret[15] = -vec[0] +vec[2];
ret[16] = -vec[1] -vec[2];
ret[17] = -vec[1] +vec[2];
ret[18] = 0.0;
return ret;
}

template <typename T>
inline std::valarray<lbBase_t> D3Q19::grad(const T& rho)
{
std::valarray<lbBase_t> ret(nD);
ret[0] =+ w1c2Inv * ( + rho[0] - rho[9] ) + w2c2Inv * ( + rho[3] + rho[4] + rho[5] + rho[6] - rho[12] - rho[13] - rho[14] - rho[15] ) ;
ret[1] =+ w1c2Inv * ( + rho[1] - rho[10] ) + w2c2Inv * ( + rho[3] - rho[4] + rho[7] + rho[8] - rho[12] + rho[13] - rho[16] - rho[17] ) ;
ret[2] =+ w1c2Inv * ( + rho[2] - rho[11] ) + w2c2Inv * ( + rho[5] - rho[6] + rho[7] - rho[8] - rho[14] + rho[15] - rho[16] + rho[17] ) ;
return ret;
}

template <typename T>
inline lbBase_t D3Q19::divGrad(const T& rho)
{
lbBase_t ret;
ret =+ 2*( w0c2Inv - c2Inv ) * ( + rho[18] ) + 2* w1c2Inv * ( + rho[0] + rho[1] + rho[2] + rho[9] + rho[10] + rho[11] ) + 2* w2c2Inv * ( + rho[3] + rho[4] + rho[5] + rho[6] + rho[7] + rho[8] + rho[12] + rho[13] + rho[14] + rho[15] + rho[16] + rho[17] ) ;
return ret;
}

template <typename T>
inline lbBase_t D3Q19::qSum(const T &dist)
{
lbBase_t ret = 0.0;
for (int q = 0; q < nQ; ++q)
ret += dist[q];
return ret;
}

template <typename T>
inline std::valarray<lbBase_t> D3Q19::qSumC(const T &dist)
{
std::valarray<lbBase_t> ret(nD);
ret[0] = + dist[0] + dist[3] + dist[4] + dist[5] + dist[6] - dist[9] - dist[12] - dist[13] - dist[14] - dist[15];
ret[1] = + dist[1] + dist[3] - dist[4] + dist[7] + dist[8] - dist[10] - dist[12] + dist[13] - dist[16] - dist[17];
ret[2] = + dist[2] + dist[5] - dist[6] + dist[7] - dist[8] - dist[11] - dist[14] + dist[15] - dist[16] + dist[17];
return ret;
}

template <typename T>
inline std::valarray<lbBase_t> D3Q19::qSumCCLowTri(const T &dist)
{
std::valarray<lbBase_t> ret(nD*(nD+1)/2);
ret[0] = + dist[0] + dist[3] + dist[4] + dist[5] + dist[6] + dist[9] + dist[12] + dist[13] + dist[14] + dist[15];
ret[1] = + dist[3] - dist[4] + dist[12] - dist[13];
ret[2] = + dist[1] + dist[3] + dist[4] + dist[7] + dist[8] + dist[10] + dist[12] + dist[13] + dist[16] + dist[17];
ret[3] = + dist[5] - dist[6] + dist[14] - dist[15];
ret[4] = + dist[7] - dist[8] + dist[16] - dist[17];
ret[5] = + dist[2] + dist[5] + dist[6] + dist[7] + dist[8] + dist[11] + dist[14] + dist[15] + dist[16] + dist[17];
return ret;
}

inline int D3Q19::c2q(const std::vector<int> &v)
/*
* returns the lattice direction that corresponds to the vector v.
* returns -1 if the vector is not found amongs the lattice vectors.
*/
{
for (int q = 0; q < nQ; ++q) {
std::vector<int> cq(cDMajor_ + nD*q, cDMajor_ + nD*q + nD);
if (cq == v) {
return q;
}
}
return -1;
}

template <typename T>
inline lbBase_t D3Q19::traceLowTri(const T &lowTri)
{
lbBase_t ret;
return ret =+ lowTri[0]+ lowTri[2]+ lowTri[5];
}

template <typename T>
inline lbBase_t D3Q19::traceOfMatrix(const T &mat)
{
lbBase_t ret;
return ret =+ mat[0]+ mat[4]+ mat[8];
}

inline std::valarray<lbBase_t> D3Q19::deltaLowTri()
{
std::valarray<lbBase_t> ret(nD*(nD+1)/2);
ret[0] = 1;
ret[1] = 0;
ret[2] = 1;
ret[3] = 0;
ret[4] = 0;
ret[5] = 1;
return ret;
}

inline std::valarray<lbBase_t> D3Q19::deltaMatrix()
{
std::valarray<lbBase_t> ret(nD*nD);
ret[0] = 1;
ret[1] = 0;
ret[2] = 0;
ret[3] = 0;
ret[4] = 1;
ret[5] = 0;
ret[6] = 0;
ret[7] = 0;
ret[8] = 1;
return ret;
}

template <typename T1, typename T2>
inline lbBase_t D3Q19::contractionLowTri(const T1 &lowTri1, const T2 &lowTri2)
{
lbBase_t ret;
return ret =+ lowTri1[0]*lowTri2[0]+ 2*lowTri1[1]*lowTri2[1]+ lowTri1[2]*lowTri2[2]+ 2*lowTri1[3]*lowTri2[3]+ 2*lowTri1[4]*lowTri2[4]+ lowTri1[5]*lowTri2[5];
}

template <typename T>
inline lbBase_t D3Q19::contractionRank2(const T &mat1, const T &mat2)
{
lbBase_t ret;
return ret =+ mat1[0]*mat2[0]+ mat1[1]*mat2[1]+ mat1[2]*mat2[2]+ mat1[3]*mat2[3]+ mat1[4]*mat2[4]+ mat1[5]*mat2[5]+ mat1[6]*mat2[6]+ mat1[7]*mat2[7]+ mat1[8]*mat2[8];
}

template <typename T>
inline std::valarray<lbBase_t> D3Q19::matrixMultiplication(const T &mat1, const T &mat2)
{
std::valarray<lbBase_t> ret(nD*nD);
ret[0] = + mat1[0]*mat2[0] + mat1[1]*mat2[3] + mat1[2]*mat2[6];
ret[1] = + mat1[0]*mat2[1] + mat1[1]*mat2[4] + mat1[2]*mat2[7];
ret[2] = + mat1[0]*mat2[2] + mat1[1]*mat2[5] + mat1[2]*mat2[8];
ret[3] = + mat1[3]*mat2[0] + mat1[4]*mat2[3] + mat1[5]*mat2[6];
ret[4] = + mat1[3]*mat2[1] + mat1[4]*mat2[4] + mat1[5]*mat2[7];
ret[5] = + mat1[3]*mat2[2] + mat1[4]*mat2[5] + mat1[5]*mat2[8];
ret[6] = + mat1[6]*mat2[0] + mat1[7]*mat2[3] + mat1[8]*mat2[6];
ret[7] = + mat1[6]*mat2[1] + mat1[7]*mat2[4] + mat1[8]*mat2[7];
ret[8] = + mat1[6]*mat2[2] + mat1[7]*mat2[5] + mat1[8]*mat2[8];
return ret;
}

template <typename T1, typename T2>
inline std::valarray<lbBase_t> D3Q19::contractionLowTriVec(const T1 &lowTri, const T2 &vec)
{
std::valarray<lbBase_t> ret(nD);
ret[0] = + lowTri[0]*vec[0] + lowTri[1]*vec[1] + lowTri[3]*vec[2];
ret[1] = + lowTri[1]*vec[0] + lowTri[2]*vec[1] + lowTri[4]*vec[2];
ret[2] = + lowTri[3]*vec[0] + lowTri[4]*vec[1] + lowTri[5]*vec[2];
return ret;
}


#endif // LBD3Q19_H
