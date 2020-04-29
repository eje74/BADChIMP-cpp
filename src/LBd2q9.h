#ifndef LBD2Q9_H
#define LBD2Q9_H

#include "LBglobal.h"
#include "LBfield.h"
#include <vector>

// See "LBlatticetypes.h" for description of the structure

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

template <typename T>
inline static std::valarray<lbBase_t> qSumCCLowTri(const T &dist);

template <typename T>
inline static lbBase_t traceLowTri(const T &lowTri);

template <typename T>
inline static lbBase_t contractionLowTri(const T &lowTri1, const T &lowTri2);

// Two phase
static void gradPush(const lbBase_t& scalarVal, const int* neighList, VectorField<D2Q9>& grad);

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

template <typename T>
inline lbBase_t D2Q9::traceLowTri(const T &lowTri)
{
lbBase_t ret;
return ret =+ lowTri[0]+ lowTri[2];
}

template <typename T>
inline lbBase_t D2Q9::contractionLowTri(const T &lowTri1, const T &lowTri2)
{
lbBase_t ret;
return ret =+ lowTri1[0]*lowTri2[0]+ 2*lowTri1[1]*lowTri2[1]+ lowTri1[2]*lowTri2[2];
}

inline void D2Q9::gradPush(const lbBase_t &scalarVal, const int *neigList, VectorField<D2Q9> &grad)
{
const lbBase_t valTmp1  = scalarVal * c2Inv * w1;
const lbBase_t valTmp2  = scalarVal * c2Inv * w2;

int nodeNeigNo = neigList[0];
grad(0,0,nodeNeigNo) -= valTmp1;

nodeNeigNo = neigList[1];
grad(0,0,nodeNeigNo) -= valTmp2;
grad(0,1,nodeNeigNo) -= valTmp2;

nodeNeigNo = neigList[2];
grad(0,1,nodeNeigNo) -= valTmp1;

nodeNeigNo = neigList[3];
grad(0,0,nodeNeigNo) += valTmp2;
grad(0,1,nodeNeigNo) -= valTmp2;

nodeNeigNo = neigList[4];
grad(0,0,nodeNeigNo) += valTmp1;

nodeNeigNo = neigList[5];
grad(0,0,nodeNeigNo) += valTmp2;
grad(0,1,nodeNeigNo) += valTmp2;

nodeNeigNo = neigList[6];
grad(0,1,nodeNeigNo) += valTmp1;

nodeNeigNo = neigList[7];
grad(0,0,nodeNeigNo) -= valTmp2;
grad(0,1,nodeNeigNo) += valTmp2;

}


#endif // LBD2Q9_H
