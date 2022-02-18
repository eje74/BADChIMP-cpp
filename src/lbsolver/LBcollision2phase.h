#ifndef LBCOLLISION2PHASE_H
#define LBCOLLISION2PHASE_H

#include "LBglobal.h"
#include "LBlatticetypes.h"

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaST(const lbBase_t &tau, const lbBase_t &sigma, const lbBase_t &CGNorm, const std::valarray<lbBase_t> &cCGNorm)
{
  //A= DXQY::c4Inv /4 * sigma/tau
  //1.125 = 0.5 * DXQY::c4Inv / (2*2)
    lbBase_t AF0_5 = 1.125 * CGNorm * sigma / tau;
    std::valarray<lbBase_t> ret(DXQY::nQ);
    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
        ret[q] = AF0_5 * (DXQY::w[q] * cCGNorm[q]*cCGNorm[q] - DXQY::B[q] );
    }
    ret[DXQY::nQNonZero_] = -AF0_5 * DXQY::B[DXQY::nQNonZero_];

    return ret;
}
template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaSTDiffCorr(const lbBase_t &tau, const std::valarray<lbBase_t> &cTgradphi)
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_factor = (1 - 0.5 / tau);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = DXQY::w[q]*tau_factor * (DXQY::c2Inv*cTgradphi[q]);
    }
    return ret;
    
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaST2(const lbBase_t &tau, const lbBase_t &sigma,  const std::valarray<lbBase_t> &cu, const lbBase_t &CGNorm, const std::valarray<lbBase_t> &cCGNorm)
{
    lbBase_t AF0_5 = 2.25 * CGNorm * sigma / tau;
    std::valarray<lbBase_t> ret(DXQY::nQ);
    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
      ret[q] = AF0_5 * (DXQY::w[q] * cCGNorm[q]*cCGNorm[q] - DXQY::B[q] /**(1+DXQY::c2Inv *cu[q])*/);
    }
    ret[DXQY::nQNonZero_] = -AF0_5 * DXQY::B[DXQY::nQNonZero_];

    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaRC(const lbBase_t &beta, const lbBase_t &rho0, const lbBase_t &rho1, const lbBase_t &rho, const std::valarray<lbBase_t> &cCGNorm)
{
    std::valarray<lbBase_t> ret(DXQY::nQ);

    lbBase_t rhoFacBeta = beta * rho0 * rho1 / rho;
    //lbBase_t rhoFacBeta = beta * rho0 * rho1* rho;
    //lbBase_t rhoFacBeta = beta * rho0 * rho1;

    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
      ret[q] = rhoFacBeta * DXQY::w[q] * cCGNorm[q] /  DXQY::cNorm[q];
      //ret[q] = rhoFacBeta * DXQY::w[q] * cCGNorm[q]; // Removed normalization of individual lattice direction vector.
    }
    ret[DXQY::nQNonZero_] = 0.0; // This should be zero by default

    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaRC2(const lbBase_t &beta, const lbBase_t &rho0, const lbBase_t &rho1, const lbBase_t &rho, const std::valarray<lbBase_t> &cCGNorm)
{
    std::valarray<lbBase_t> ret(DXQY::nQ);

    lbBase_t rhoFacBeta = beta * rho0 * rho1 / rho;
    //lbBase_t rhoFacBeta = beta * rho0 * rho1* rho;
    //lbBase_t rhoFacBeta = beta * rho0 * rho1;

    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
      //ret[q] = rhoFacBeta * DXQY::w[q] * cCGNorm[q] /  DXQY::cNorm[q];
      ret[q] = rhoFacBeta * DXQY::w[q] * cCGNorm[q]; // Removed normalization of individual lattice direction vector.
    }
    ret[DXQY::nQNonZero_] = 0.0; // This should be zero by default

    return ret;
}
template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaRC3(const lbBase_t &beta, const lbBase_t &rho0, const lbBase_t &rho1, const lbBase_t &rho, const std::valarray<lbBase_t> &cu, const lbBase_t &uCGNorm, const std::valarray<lbBase_t> &cCGNorm)
{
    std::valarray<lbBase_t> ret(DXQY::nQ);

    lbBase_t rhoFacBeta = beta * rho0 * rho1 / rho;
    //lbBase_t rhoFacBeta = beta * rho0 * rho1* rho;
    //lbBase_t rhoFacBeta = beta * rho0 * rho1;

    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
      //ret[q] = rhoFacBeta * DXQY::w[q] * cCGNorm[q] /  DXQY::cNorm[q];
      ret[q] = rhoFacBeta * DXQY::w[q] * (cCGNorm[q]+ DXQY::c2Inv * ( cCGNorm[q] * cu[q] - DXQY::c2 * uCGNorm)); // Removed normalization of individual lattice direction vector.
    }
    ret[DXQY::nQNonZero_] = 0.0; // This should be zero by default

    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaRCInd(const lbBase_t &beta, const lbBase_t &indicator, const lbBase_t &rho, const std::valarray<lbBase_t> &cCGNorm)
{
    std::valarray<lbBase_t> ret(DXQY::nQ);

    lbBase_t rhoFacBeta = beta * indicator*(1-indicator)*rho;

    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
        ret[q] = rhoFacBeta * DXQY::w[q] * cCGNorm[q] /  DXQY::cNorm[q];
	//std::cout<<ret[q]<<std::endl;
    }
    ret[DXQY::nQNonZero_] = 0.0; // This should be zero by default

    return ret;
}

#endif // LBCOLLISION2PHASE_H
