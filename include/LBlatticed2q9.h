#ifndef LBLATTICED2Q9_H
#define LBLATTICED2Q9_H

#include "LBglobal.h"

// LATTICED2Q9_H
class D2Q9
{
public:
    inline int nD() const {return nD_;}
    inline int nQ() const {return nQ_;}
    inline lbbase_t w(const int qDirection) const  {return w_[qDirection];}
    inline int c(const int qDirection, const int dimension) const  {return cDMajor_[nD_*qDirection + dimension];}
    inline int reverseDirection(const int qDirection) const  {return (qDirection + reverseStep_) % nNonZeroDirections_;}

    lbbase_t dot(const lbbase_t* leftVec, const lbbase_t* rightVec) const;
    lbbase_t cDot(const int qDir, const lbbase_t* rightVec) const;
    void cDotAll(const lbbase_t* rightVec, lbbase_t* ret) const;

    lbbase_t qSum(const lbbase_t* dist) const;
    lbbase_t qSumC(const lbbase_t* dist, lbbase_t* ret) const;


private:
    static constexpr int nD_ = 2;
    static constexpr int nQ_ = 9;
    static constexpr int reverseStep_ = 4;
    static constexpr int nNonZeroDirections_ = 8;

    static constexpr lbbase_t w_[nQ_] = {1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 4.0/9.0};
    static constexpr int cDMajor_[nQ_*nD_] = {1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 0, 0};
    static constexpr int cQMajor_[nQ_*nD_] = {1, 1, 0, -1, -1, -1, 0, 1, 0, 0, 1, 1, 1, 0, -1, -1, -1, 0};

    static constexpr lbbase_t c2Inv_ = 3.0;
    static constexpr lbbase_t c2_ = 1.0 / c2Inv_;
    static constexpr lbbase_t c4Inv_ = 9.0;
    static constexpr lbbase_t c4Inv0_5_ = 4.5;
    static constexpr lbbase_t c4_ = 1.0 / c4Inv_;
};

inline lbbase_t D2Q9::dot(const lbbase_t* leftVec, const lbbase_t* rightVec) const
{
    return leftVec[0]*rightVec[0] + leftVec[1]*rightVec[1];
}

inline lbbase_t D2Q9::cDot(const int qDir, const lbbase_t* rightVec) const
{
    return c(qDir, 0)*rightVec[0] + c(qDir, 1)*rightVec[1];
}

inline lbbase_t D2Q9::qSum(const lbbase_t* dist) const
{
    return dist[0] + dist[1]+ dist[2]+ dist[3]+ dist[4]+ dist[5]+ dist[6]+ dist[7]+ dist[8];
}


// END LATTICED2Q9_H

#endif
