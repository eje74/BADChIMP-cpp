#ifndef LBRANS_H
#define LBRANS_H

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "../lbsolver/LBglobal.h"
#include "../io/Input.h"

//=====================================================================================
//
//                                      R A N S
//
//=====================================================================================
template<typename DXQY>
class Rans
{
public:
    Rans() {};
    Rans(Input &input);

    inline lbBase_t tau() const {return tau_;}
    inline lbBase_t tauK() const {return tauK_;}
    inline lbBase_t tauE() const {return tauE_;}

    inline lbBase_t sourceK() const {return sourceK_;}
    inline lbBase_t sourceE() const {return sourceE_;}

    inline lbBase_t rhoK() const {return rhoK_;}
    inline lbBase_t rhoE() const {return rhoE_;}

private:
    lbBase_t tau_;
    lbBase_t tauK_;
    lbBase_t tauE_;

    lbBase_t sourceK_;
    lbBase_t sourceE_;

    lbBase_t rhoK_;
    lbBase_t rhoE_;
};

//                                        Rans
//------------------------------------------------------------------------------------- Rans
template<typename DXQY>
Rans<DXQY>::Rans(Input &input)
{

}

#endif