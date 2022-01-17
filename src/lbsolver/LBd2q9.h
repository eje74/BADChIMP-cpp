#ifndef LBD2Q9_H
#define LBD2Q9_H

#include "LBglobal.h"
#include <vector>

// See "LBlatticetypes.h" for description of the structure

// TO MAKE CHANGES TO THIS FILE, MAKE THE CHANGES IN "PythonScripts/writeLatticeFile.py",
// RUN THE SCRIPT AND PLACE RESULTING FILES IN "src/"

struct D2Q9{

static constexpr int nD = 2;
static constexpr int nQ = 9;
static constexpr int nDirPairs_ = 4;
static constexpr int nQNonZero_ = 8;

static constexpr lbBase_t c2Inv = 3.0;
static constexpr lbBase_t c4Inv = 9.0;
static constexpr lbBase_t c2 = 1.0 / c2Inv;
static constexpr lbBase_t c4 = 1.0 / c4Inv;
static constexpr lbBase_t c4Inv0_5 = 0.5 * c4Inv;

static constexpr lbBase_t w0 = 16.0/36.0;
static constexpr lbBase_t w0c2Inv = w0*c2Inv;
static constexpr lbBase_t w1 = 4.0/36.0;
static constexpr lbBase_t w1c2Inv = w1*c2Inv;
static constexpr lbBase_t w2 = 1.0/36.0;
static constexpr lbBase_t w2c2Inv = w2*c2Inv;

static constexpr lbBase_t w[9] = {w1, w2, w1, w2, w1, w2, w1, w2, w0};
static constexpr int cDMajor_[18] = {1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 0, 0};
static constexpr lbBase_t cNorm[9] = {1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 0.0};
static constexpr int reverseDirection_[9] = {4, 5, 6, 7, 0, 1, 2, 3, 8};
static constexpr lbBase_t B0 = -16.0/108.0;
static constexpr lbBase_t B1 = 8.0/108.0;
static constexpr lbBase_t B2 = 5.0/108.0;
static constexpr lbBase_t B[9] = {B1, B2, B1, B2, B1, B2, B1, B2, B0};

static constexpr lbBase_t UnitMatrixLowTri[3] = {1, 0, 1};

// Functions

inline static int c(const int qDirection, const int dimension)  {return cDMajor_[nD*qDirection + dimension];}
inline static int reverseDirection(const int qDirection) {return reverseDirection_[qDirection];}

inline static std::vector<int> c(const int qDirection) {
std::vector<int> cq(cDMajor_ + nD*qDirection, cDMajor_ + nD*qDirection + nD);
return cq;
}

template <typename T1, typename T2>
inline static lbBase_t dot(const T1 &leftVec, const T2 &rightVec);
template<typename T>
inline static T cDot(const int qDir, const T* rightVec);
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
inline lbBase_t D2Q9::dot(const T1 &leftVec, const T2 &rightVec)
{
    return leftVec[0]*rightVec[0] + leftVec[1]*rightVec[1];
}

template<typename T>
inline T D2Q9::cDot(const int qDir, const T* rightVec)
{
    return c(qDir, 0)*rightVec[0] + c(qDir, 1)*rightVec[1];
}

template <typename T>
inline std::valarray<lbBase_t> D2Q9::cDotAll(const T &vec)
{
std::valarray<lbBase_t> ret(nQ);
ret[0] = +vec[0];
ret[1] = +vec[0] +vec[1];
ret[2] = +vec[1];
ret[3] = -vec[0] +vec[1];
ret[4] = -vec[0];
ret[5] = -vec[0] -vec[1];
ret[6] = -vec[1];
ret[7] = +vec[0] -vec[1];
ret[8] = 0.0;
return ret;
}

template <typename T>
inline std::valarray<lbBase_t> D2Q9::grad(const T& rho)
{
std::valarray<lbBase_t> ret(nD);
ret[0] =+ w1c2Inv * ( + rho[0] - rho[4] ) + w2c2Inv * ( + rho[1] - rho[3] - rho[5] + rho[7] ) ;
ret[1] =+ w1c2Inv * ( + rho[2] - rho[6] ) + w2c2Inv * ( + rho[1] + rho[3] - rho[5] - rho[7] ) ;
return ret;
}

template <typename T>
inline lbBase_t D2Q9::divGrad(const T& rho)
{
lbBase_t ret;
ret =+ 2*( w0c2Inv - c2Inv ) * ( + rho[8] ) + 2* w1c2Inv * ( + rho[0] + rho[2] + rho[4] + rho[6] ) + 2* w2c2Inv * ( + rho[1] + rho[3] + rho[5] + rho[7] ) ;
return ret;
}

template <typename T>
inline lbBase_t D2Q9::qSum(const T &dist)
{
lbBase_t ret = 0.0;
for (int q = 0; q < nQ; ++q)
ret += dist[q];
return ret;
}

template <typename T>
inline std::valarray<lbBase_t> D2Q9::qSumC(const T &dist)
{
std::valarray<lbBase_t> ret(nD);
ret[0] = + dist[0] + dist[1] - dist[3] - dist[4] - dist[5] + dist[7];
ret[1] = + dist[1] + dist[2] + dist[3] - dist[5] - dist[6] - dist[7];
return ret;
}

template <typename T>
inline std::valarray<lbBase_t> D2Q9::qSumCCLowTri(const T &dist)
{
std::valarray<lbBase_t> ret(nD*(nD+1)/2);
ret[0] = + dist[0] + dist[1] + dist[3] + dist[4] + dist[5] + dist[7];
ret[1] = + dist[1] - dist[3] + dist[5] - dist[7];
ret[2] = + dist[1] + dist[2] + dist[3] + dist[5] + dist[6] + dist[7];
return ret;
}

inline int D2Q9::c2q(const std::vector<int> &v)
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
inline lbBase_t D2Q9::traceLowTri(const T &lowTri)
{
lbBase_t ret;
return ret =+ lowTri[0]+ lowTri[2];
}

template <typename T>
inline lbBase_t D2Q9::traceOfMatrix(const T &mat)
{
lbBase_t ret;
return ret =+ mat[0]+ mat[3];
}

inline std::valarray<lbBase_t> D2Q9::deltaLowTri()
{
std::valarray<lbBase_t> ret(nD*(nD+1)/2);
ret[0] = 1;
ret[1] = 0;
ret[2] = 1;
return ret;
}

inline std::valarray<lbBase_t> D2Q9::deltaMatrix()
{
std::valarray<lbBase_t> ret(nD*nD);
ret[0] = 1;
ret[1] = 0;
ret[2] = 0;
ret[3] = 1;
return ret;
}

template <typename T1, typename T2>
inline lbBase_t D2Q9::contractionLowTri(const T1 &lowTri1, const T2 &lowTri2)
{
lbBase_t ret;
return ret =+ lowTri1[0]*lowTri2[0]+ 2*lowTri1[1]*lowTri2[1]+ lowTri1[2]*lowTri2[2];
}

template <typename T>
inline lbBase_t D2Q9::contractionRank2(const T &mat1, const T &mat2)
{
lbBase_t ret;
return ret =+ mat1[0]*mat2[0]+ mat1[1]*mat2[1]+ mat1[2]*mat2[2]+ mat1[3]*mat2[3];
}

template <typename T>
inline std::valarray<lbBase_t> D2Q9::matrixMultiplication(const T &mat1, const T &mat2)
{
std::valarray<lbBase_t> ret(nD*nD);
ret[0] = + mat1[0]*mat2[0] + mat1[1]*mat2[2];
ret[1] = + mat1[0]*mat2[1] + mat1[1]*mat2[3];
ret[2] = + mat1[2]*mat2[0] + mat1[3]*mat2[2];
ret[3] = + mat1[2]*mat2[1] + mat1[3]*mat2[3];
return ret;
}

template <typename T1, typename T2>
inline std::valarray<lbBase_t> D2Q9::contractionLowTriVec(const T1 &lowTri, const T2 &vec)
{
std::valarray<lbBase_t> ret(nD);
ret[0] = + lowTri[0]*vec[0] + lowTri[1]*vec[1];
ret[1] = + lowTri[1]*vec[0] + lowTri[2]*vec[1];
return ret;
}


#endif // LBD2Q9_H
