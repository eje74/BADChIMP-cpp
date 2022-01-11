#ifndef LBDIFFUSION_H
#define LBDIFFUSION_H

#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBvtk.h"
#include "../lbsolver/LBnodes.h"



template<typename DXQY>
class DiffusionSolver
{
public:
    DiffusionSolver(const lbBase_t tau);
    std::valarray<lbBase_t> setF(const lbBase_t & rho) const;
    template<typename T>
    std::valarray<lbBase_t> omegaBGK(const T & f, const lbBase_t & rho) const;
    lbBase_t diffusionCoefficient() const;
private:
    const lbBase_t tau_;
    const lbBase_t tauInv_;
    const std::valarray<lbBase_t> w_;
};


template<typename DXQY>
DiffusionSolver<DXQY>::DiffusionSolver(const lbBase_t tau): tau_(tau), tauInv_(1.0/tau), w_(DXQY::w, DXQY::nQ) {}

template<typename DXQY>
std::valarray<lbBase_t> DiffusionSolver<DXQY>::setF(const lbBase_t & rho) const
{
    return w_*rho;
}

template<typename DXQY>
template<typename T>
std::valarray<lbBase_t> DiffusionSolver<DXQY>::omegaBGK(const T & f, const lbBase_t & rho) const
{
    return tauInv_*(w_*rho - f);
}

template<typename DXQY>
lbBase_t DiffusionSolver<DXQY>::diffusionCoefficient() const
{
    return DXQY::c2 * (tau_ - 0.5);
}

#endif