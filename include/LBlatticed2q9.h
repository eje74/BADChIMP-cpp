#ifndef LBLATTICED2Q9_H
#define LBLATTICED2Q9_H

#include "LBglobal.h"

// LATTICED2Q9_H
class D2Q9
{
public:
    inline int nDremove() const  {return nD;}
    inline int nQremove() const  {return nQ;}

    inline lbbase_t wRemove(const int qDirection) const {return w[qDirection];}
    inline int c(const int qDirection, const int dimension) const  {return cDMajor_[nD*qDirection + dimension];}

    inline int reverseDirection(const int qDirection) const  {return (qDirection + reverseStep_) % nNonZeroDirections_;}

    lbbase_t dot(const lbbase_t* leftVec, const lbbase_t* rightVec) const;
    lbbase_t cDot(const int qDir, const lbbase_t* rightVec) const;
    void cDotAll(const lbbase_t* vec, lbbase_t* ret) const;

    lbbase_t qSum(const lbbase_t* dist) const;
    void qSumC(const lbbase_t* dist, lbbase_t* ret) const;

    static constexpr int nD = 2;
    static constexpr int nQ = 9;

    static constexpr lbbase_t c2Inv = 3.0;
    static constexpr lbbase_t c2 = 1.0 / c2Inv;
    static constexpr lbbase_t c4Inv = 9.0;
    static constexpr lbbase_t c4Inv0_5 = 4.5;
    static constexpr lbbase_t c4 = 1.0 / c4Inv;

    static constexpr lbbase_t w[9] = {1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 4.0/9.0};

private:
    static constexpr int reverseStep_ = 4;
    static constexpr int nNonZeroDirections_ = 8;

    static constexpr int cDMajor_[18] = {1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 0, 0};
    static constexpr int cQMajor_[18] = {1, 1, 0, -1, -1, -1,  0,  1, 0,
                                         0, 1, 1,  1,  0, -1, -1, -1, 0};
};

inline lbbase_t D2Q9::dot(const lbbase_t* leftVec, const lbbase_t* rightVec) const
{
    return leftVec[0]*rightVec[0] + leftVec[1]*rightVec[1];
}

inline lbbase_t D2Q9::cDot(const int qDir, const lbbase_t* rightVec) const
{
    return c(qDir, 0)*rightVec[0] + c(qDir, 1)*rightVec[1];
}

inline void D2Q9::cDotAll(const lbbase_t* vec, lbbase_t* ret) const
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


inline lbbase_t D2Q9::qSum(const lbbase_t* dist) const
{
    lbbase_t ret = 0.0;

    for (int q = 0; q < nQ; ++q)
        ret += dist[q];
    return ret;
}

inline void D2Q9::qSumC(const lbbase_t* dist, lbbase_t* ret) const
{
    ret[0] = dist[0] + dist[1]            - dist[3] - dist[4]  - dist[5]           + dist[7];
    ret[1] =           dist[1] + dist[2]  + dist[3]            - dist[5] - dist[6] - dist[7];
}

#endif
