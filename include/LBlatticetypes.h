#ifndef LBLATTICETYPES_H
#define LBLATTICETYPES_H

#include "LBglobal.h"

struct D2Q9 {

static constexpr int nD = 2;
static constexpr int nQ = 9;
static constexpr int nDirPairs_ = 4; // Number of unsigned bounce back directions
static constexpr int nQNonZero_ = 8; // Number of non zeros basis directions

static constexpr lbBase_t c2Inv = 3.0;
static constexpr lbBase_t c2 = 1.0 / c2Inv;
static constexpr lbBase_t c4Inv = 9.0;
static constexpr lbBase_t c4Inv0_5 = 4.5;
static constexpr lbBase_t c4 = 1.0 / c4Inv;

static constexpr lbBase_t w[9] = {1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 4.0/9.0};
static constexpr int cDMajor_[18] = {1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 0, 0};


inline static int c(const int qDirection, const int dimension)  {return cDMajor_[nD*qDirection + dimension];}
inline static int reverseDirection(const int qDirection) {return (qDirection + nDirPairs_) % nQNonZero_;}

static lbBase_t dot(const lbBase_t* leftVec, const lbBase_t* rightVec);
static lbBase_t cDot(const int qDir, const lbBase_t* rightVec);
static void cDotAll(const lbBase_t* vec, lbBase_t* ret);

static void qSum(const lbBase_t* dist, lbBase_t& ret);
static void qSumC(const lbBase_t* dist, lbBase_t* ret);
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


#endif // LBLATTICETYPES_H
