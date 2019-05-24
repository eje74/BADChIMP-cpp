#ifndef LBCOLLISION2PHASE_H
#define LBCOLLISION2PHASE_H

#include "LBglobal.h"
#include "LBlatticetypes.h"

template <typename DXQY>
inline std::vector<lbBase_t> calcDeltaOmegaST(const lbBase_t &tau, const lbBase_t &sigma, const lbBase_t &CGNorm, const std::vector<lbBase_t> &cCGNorm)
{
    lbBase_t AF0_5 = 2.25 * CGNorm * sigma / tau;
    std::vector<lbBase_t> ret(DXQY::nQ);
    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
        ret[q] = AF0_5 * (DXQY::w[q] * cCGNorm[q]*cCGNorm[q] - DXQY::B[q] );
    }
    ret[DXQY::nQNonZero_] = -AF0_5 * DXQY::B[DXQY::nQNonZero_];

    return ret;
}


template <typename DXQY>
inline std::vector<lbBase_t> calcDeltaOmegaRC(const lbBase_t &beta, const lbBase_t &rho0, const lbBase_t &rho1, const lbBase_t &rho, const std::vector<lbBase_t> &cCGNorm)
{
    std::vector<lbBase_t> ret(DXQY::nQ);

    lbBase_t rhoFacBeta = beta * rho0 * rho1 / rho;

    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
        ret[q] = rhoFacBeta * DXQY::w[q] * cCGNorm[q] /  DXQY::cNorm[q];
    }
    ret[DXQY::nQNonZero_] = 0.0; // This should be zero by default

    return ret;
}
#endif // LBCOLLISION2PHASE_H