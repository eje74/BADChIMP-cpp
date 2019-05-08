#ifndef LBD3Q19_H
#define LBD3Q19_H

#include "LBglobal.h"
#include "LBfield.h"
#include <vector>


// See "LBlatticetypes.h" for description of the structure

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
static constexpr lbBase_t B0 = -12.0/54.0;
static constexpr lbBase_t B1 = 1.0/54.0;
static constexpr lbBase_t B2 = 2.0/54.0;
static constexpr lbBase_t B[19] = {B1, B1, B1, B2, B2, B2, B2, B2, B2, B1, B1, B1, B2, B2, B2, B2, B2, B2, B0};

// Functions

inline static int c(const int qDirection, const int dimension)  {return cDMajor_[nD*qDirection + dimension];}
inline static int reverseDirection(const int qDirection) {return (qDirection + nDirPairs_) % nQNonZero_;}

template <typename T1, typename T2>
inline static lbBase_t dot(const T1 &leftVec, const T2 &rightVec);
template<typename T>
inline static T cDot(const int qDir, const T* rightVec);
template <typename T>
inline static std::vector<lbBase_t> cDotAll(const T &vec);
template <typename T>
inline static std::vector<lbBase_t> grad(const T &rho);

template <typename T>
inline static lbBase_t qSum(const T &dist);
template <typename T>
inline static std::vector<lbBase_t> qSumC(const T &dist);

// Two phase
static void gradPush(const lbBase_t& scalarVal, const int* neighList, VectorField<D3Q19>& grad);

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

template <typename T>
inline std::vector<lbBase_t> D3Q19::cDotAll(const T &vec)
{
std::vector<lbBase_t> ret(nQ);
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
inline std::vector<lbBase_t> D3Q19::grad(const T& rho)
{
std::vector<lbBase_t> ret(nD);
ret[0] =+ w1c2Inv * ( + rho[0] - rho[9] ) + w2c2Inv * ( + rho[3] + rho[4] + rho[5] + rho[6] - rho[12] - rho[13] - rho[14] - rho[15] ) ;
ret[1] =+ w1c2Inv * ( + rho[1] - rho[10] ) + w2c2Inv * ( + rho[3] - rho[4] + rho[7] + rho[8] - rho[12] + rho[13] - rho[16] - rho[17] ) ;
ret[2] =+ w1c2Inv * ( + rho[2] - rho[11] ) + w2c2Inv * ( + rho[5] - rho[6] + rho[7] - rho[8] - rho[14] + rho[15] - rho[16] + rho[17] ) ;
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
inline std::vector<lbBase_t> D3Q19::qSumC(const T &dist)
{
std::vector<lbBase_t> ret(nD);
ret[0] = + dist[0] + dist[3] + dist[4] + dist[5] + dist[6] - dist[9] - dist[12] - dist[13] - dist[14] - dist[15];
ret[1] = + dist[1] + dist[3] - dist[4] + dist[7] + dist[8] - dist[10] - dist[12] + dist[13] - dist[16] - dist[17];
ret[2] = + dist[2] + dist[5] - dist[6] + dist[7] - dist[8] - dist[11] - dist[14] + dist[15] - dist[16] + dist[17];
return ret;
}

inline void D3Q19::gradPush(const lbBase_t &scalarVal, const int *neigList, VectorField<D3Q19> &grad)
{
const lbBase_t valTmp1  = scalarVal * c2Inv * w1;
const lbBase_t valTmp2  = scalarVal * c2Inv * w2;

int nodeNeigNo = neigList[0];
grad(0,0,nodeNeigNo) -= valTmp1;

nodeNeigNo = neigList[1];
grad(0,1,nodeNeigNo) -= valTmp1;

nodeNeigNo = neigList[2];
grad(0,2,nodeNeigNo) -= valTmp1;

nodeNeigNo = neigList[3];
grad(0,0,nodeNeigNo) -= valTmp2;
grad(0,1,nodeNeigNo) -= valTmp2;

nodeNeigNo = neigList[4];
grad(0,0,nodeNeigNo) -= valTmp2;
grad(0,1,nodeNeigNo) += valTmp2;

nodeNeigNo = neigList[5];
grad(0,0,nodeNeigNo) -= valTmp2;
grad(0,2,nodeNeigNo) -= valTmp2;

nodeNeigNo = neigList[6];
grad(0,0,nodeNeigNo) -= valTmp2;
grad(0,2,nodeNeigNo) += valTmp2;

nodeNeigNo = neigList[7];
grad(0,1,nodeNeigNo) -= valTmp2;
grad(0,2,nodeNeigNo) -= valTmp2;

nodeNeigNo = neigList[8];
grad(0,1,nodeNeigNo) -= valTmp2;
grad(0,2,nodeNeigNo) += valTmp2;

nodeNeigNo = neigList[9];
grad(0,0,nodeNeigNo) += valTmp1;

nodeNeigNo = neigList[10];
grad(0,1,nodeNeigNo) += valTmp1;

nodeNeigNo = neigList[11];
grad(0,2,nodeNeigNo) += valTmp1;

nodeNeigNo = neigList[12];
grad(0,0,nodeNeigNo) += valTmp2;
grad(0,1,nodeNeigNo) += valTmp2;

nodeNeigNo = neigList[13];
grad(0,0,nodeNeigNo) += valTmp2;
grad(0,1,nodeNeigNo) -= valTmp2;

nodeNeigNo = neigList[14];
grad(0,0,nodeNeigNo) += valTmp2;
grad(0,2,nodeNeigNo) += valTmp2;

nodeNeigNo = neigList[15];
grad(0,0,nodeNeigNo) += valTmp2;
grad(0,2,nodeNeigNo) -= valTmp2;

nodeNeigNo = neigList[16];
grad(0,1,nodeNeigNo) += valTmp2;
grad(0,2,nodeNeigNo) += valTmp2;

nodeNeigNo = neigList[17];
grad(0,1,nodeNeigNo) += valTmp2;
grad(0,2,nodeNeigNo) -= valTmp2;

}


#endif // LBD3Q19_H

