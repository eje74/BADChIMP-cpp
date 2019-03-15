#ifndef LBD2Q9_H
#define LBD2Q9_H

#include <vector>
#include "LBglobal.h"
#include "LBfield.h"


// See "LBlatticetypes.h" for description of the structure

struct D2Q9 {

static constexpr int nD = 2;
static constexpr int nQ = 9;
static constexpr int nDirPairs_ = 4;
static constexpr int nQNonZero_ = 8;

static constexpr lbBase_t c2Inv = 3.0;
static constexpr lbBase_t c2 = 1.0 / c2Inv;
static constexpr lbBase_t c4Inv = 9.0;
static constexpr lbBase_t c4 = 1.0 / c4Inv;
static constexpr lbBase_t c4Inv0_5 = 0.5 * c4Inv;

static constexpr lbBase_t w0 = 4.0/9.0;
static constexpr lbBase_t w1 = 1.0/9.0;
static constexpr lbBase_t w2 = 1.0/36.0;
static constexpr lbBase_t w0c2Inv = w0*c2Inv;
static constexpr lbBase_t w1c2Inv = w1*c2Inv;
static constexpr lbBase_t w2c2Inv = w2*c2Inv;

// Remember do define the static arrays in the cpp file as well.
static constexpr lbBase_t w[9] = {w1, w2, w1, w2, w1, w2, w1, w2, w0};
static constexpr int cDMajor_[18] = {1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 0, 0};
static constexpr lbBase_t cNorm[9] = {1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 0.0};

// Two phase values
static constexpr lbBase_t B0 = -4.0/27.0;
static constexpr lbBase_t B1 = 2.0/27.0;
static constexpr lbBase_t B2 = 5.0/108.0;

// Remember do define the static arrays in the cpp file as well.
static constexpr lbBase_t B[9] = {B1, B2, B1, B2, B1, B2, B1, B2, B0};

// Functions

inline static int c(const int qDirection, const int dimension)  {return cDMajor_[nD*qDirection + dimension];}

inline static std::vector<int> c(const int qDir)  {return std::vector<int>{cDMajor_[nD*qDir],cDMajor_[nD*qDir+1]};} // NB ??

inline static int reverseDirection(const int qDirection) {return (qDirection + nDirPairs_) % nQNonZero_;}


static lbBase_t dot(const lbBase_t* leftVec, const lbBase_t* rightVec);
/* static lbBase_t cDot(const int qDir, const lbBase_t* rightVec); */

template<typename T>
static T cDot(const int qDir, const T* rightVec);


static void cDotAll(const lbBase_t* vec, lbBase_t* ret);
static void grad(const lbBase_t* rho, lbBase_t* ret);

static void qSum(const lbBase_t* dist, lbBase_t& ret);
static void qSumC(const lbBase_t* dist, lbBase_t* ret);

// Two phase
static void gradPush(const lbBase_t& scalarVal, const int* neighList, VectorField<D2Q9>& grad);

};


inline lbBase_t D2Q9::dot(const lbBase_t* leftVec, const lbBase_t* rightVec)
{
    return leftVec[0]*rightVec[0] + leftVec[1]*rightVec[1];
}

/*inline lbBase_t D2Q9::cDot(const int qDir, const lbBase_t* rightVec)
{
    return c(qDir, 0)*rightVec[0] + c(qDir, 1)*rightVec[1];
} */

template<typename T>
inline T D2Q9::cDot(const int qDir, const T* rightVec)
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

inline void D2Q9::grad(const lbBase_t* rho, lbBase_t* ret)
{
    ret[0] =  w1c2Inv * (rho[0] - rho[4]);
    ret[0] += w2c2Inv * (rho[1] - rho[3] - rho[5] + rho[7]);
    ret[1] =  w1c2Inv * (rho[2] - rho[6]);
    ret[1] += w2c2Inv * (rho[1] + rho[3] - rho[5] - rho[7]);
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
