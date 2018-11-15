#ifndef LBLATTICED2Q9_H
#define LBLATTICED2Q9_H

#include "LBglobal.h"

static const int cx_[9] = {1, 1, 0, -1, -1, -1,  0,  1, 0};
static const int cy_[9] = {0, 1, 1,  1,  0, -1, -1, -1, 0};



class D2Q9
{
public:
    inline int nD()  {return nD_;}
    inline int nQ()  {return nQ_;}

private:
    static constexpr int nD_ = 2;
    static constexpr int nQ_ = 9;
    static constexpr lbbase_t w_[nQ_] = {1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 4.0/9.0};

    static constexpr lbbase_t c2Inv_ = 3.0;
    static constexpr lbbase_t c2_ = 1.0 / c2Inv_;
    static constexpr lbbase_t c4Inv_ = 9.0;
    static constexpr lbbase_t c4Inv0_5_ = 4.5;
    static constexpr lbbase_t c4_ = 1.0 / c4Inv_;
};

#endif // LATTICED2Q9_H
