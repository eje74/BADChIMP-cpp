#ifndef LBCOLLISION2PHASE_H
#define LBCOLLISION2PHASE_H

#include "LBglobal.h"
#include "LBlatticetypes.h"

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaST(const lbBase_t &tau, const lbBase_t &sigma, const lbBase_t &CGNorm, const std::valarray<lbBase_t> &cCGNorm)
{
  //A= DXQY::c4Inv /4 * sigma/tau
  //1.125 = 0.5 * DXQY::c4Inv / (2*2)
    const lbBase_t AF0_5 = 1.125 * CGNorm * sigma / tau;
    std::valarray<lbBase_t> ret(DXQY::nQ);
    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
        ret[q] = AF0_5 * (DXQY::w[q] * cCGNorm[q]*cCGNorm[q] - DXQY::B[q] );
    }
    ret[DXQY::nQNonZero_] = -AF0_5 * DXQY::B[DXQY::nQNonZero_];

    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaIFPressurePerturb(const lbBase_t &tau, const lbBase_t &sigma, const lbBase_t &CGNorm)
{
  //A= DXQY::c4Inv /4 * sigma/tau
  //1.125 = 0.5 * DXQY::c4Inv / (2*2)
    const lbBase_t AF0_5 = 1.125 * CGNorm * sigma / tau;
    std::valarray<lbBase_t> ret(DXQY::nQ);
    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
        ret[q] = AF0_5 * (DXQY::cNorm[q]*DXQY::cNorm[q]-DXQY::nD*DXQY::c2);
    }
    ret[DXQY::nQNonZero_] = -AF0_5 * DXQY::nD*DXQY::c2;

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

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaRC4(const lbBase_t &beta, const lbBase_t &tau, const lbBase_t &rhoTot, const lbBase_t &phi0, const lbBase_t &phi1, const lbBase_t &kappa, const lbBase_t &FNorm0, const std::valarray<lbBase_t> &cCGNorm, const std::valarray<lbBase_t> &cu, const lbBase_t &uCGNorm, const std::valarray<lbBase_t> &cJphi)
{
    std::valarray<lbBase_t> ret(DXQY::nQ);

    //lbBase_t rhoFacBeta = beta * rho0 * (1-phi0);
    lbBase_t rhoFacBeta = beta * rhoTot * phi0 * phi1;
    lbBase_t rhoInit = 1.0;
    lbBase_t tau_factor = (1 - 0.5 / tau);
    lbBase_t A = 1.0;
    lbBase_t WInv = beta*0.25/tau_factor;
    
    //lbBase_t rhoFacBeta = beta * rho0 * rho1* rho;
    //lbBase_t rhoFacBeta = beta * rho0 * rho1;

    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
      //ret[q] = DXQY::w[q] * rhoFacBeta *  cCGNorm[q] /  DXQY::cNorm[q];
      //ret[q] = DXQY::w[q]*rhoFacBeta * ( ( 1 - 2.0*DXQY::c2*beta*( 1-2*phi0 ) )*cCGNorm[q] -  DXQY::c2Inv * ( cCGNorm[q] * cu[q] - DXQY::c2 * uCGNorm ) + 0*rhoFacBeta * 0.5*DXQY::c2Inv* ( cCGNorm[q]*cCGNorm[q] - 0*DXQY::c2 ) ); // Removed normalization of individual lattice direction vector.

     
      ret[q] = 1.0;
      ret[q] *= rhoFacBeta;

      //ret[q] = rhoTot*(2*phi0-(phi0+0.5)*(phi0+0.5));
      //ret[q] = beta*rhoTot*((1+A)*phi0-(phi0+0.5*A)*(phi0+0.5*A));
      //ret[q] += beta*0.25*rhoInit*A*A;
      //ret[q] += rhoTot*phi0*DXQY::c2Inv;
      //ret[q] = rhoInit*(2*phi0-(phi0+0.5)*(phi0+0.5)) + 0.25*rhoInit;
      ret[q] *= cCGNorm[q];
      

      //ret[q] += 0.125*(beta/(1-0.5/tau))*(beta/(1-0.5/tau))*(1-2*phi0)*(1-phi0)*phi0/(1.1*1.1*(tau-0.5))  *(DXQY::cNorm[q]*DXQY::cNorm[q] - DXQY::nD*DXQY::c2);
      
      //ret[q] += -0.25*rhoInit*DXQY::c2*tau*kappa;
      
      

      
      
      //ret[q] += - 0.5*DXQY::c2 * (2*beta)*(2*beta) * ((1-2*phi0)-2*phi0*(2 - 3*phi0)) * cCGNorm[q];
      //ret[q] += - 2.0*DXQY::c2*beta*( 1-2*phi0 )*cCGNorm[q];
      //ret[q] += beta*( 1-2*phi0 )*cCGNorm[q];
      //ret[q] +=  DXQY::c2Inv * ( cCGNorm[q] * cu[q] - DXQY::c2 * uCGNorm );

      //ret[q] += beta*DXQY::c2*0.6/(1-0.5/0.6)*( 1-2*(1-phi1));
      
      //ret[q] += 0.5 * (2*beta) * (1-2*phi0) *cCGNorm[q]*cCGNorm[q];
      //ret[q] +=  + 2*0.1666666666667 * (2*beta)*(2*beta) * ((1-2*phi0)-2*phi0*(2 - 3*phi0))*cCGNorm[q]*cCGNorm[q]*cCGNorm[q];

      //ret[q] +=  + 0.1666666666667 * (2*beta)*(2*beta) * ((1-2*phi0)-2*phi0*(2 - 3*phi0))*(cCGNorm[q]*cCGNorm[q]*0 - DXQY::c2)*cCGNorm[q];
      
      //ret[q] +=  0.5*DXQY::c2 * (2*beta)*(2*beta) * ((1-2*phi0)-2*phi0*(2 - 3*phi0)) * cCGNorm[q];
      //ret[q] +=  - 0.5/1.5*DXQY::c2 * (2*beta)*(2*beta) * ((1-2*phi0)-2*phi0*(2 - 3*phi0)) * cCGNorm[q];

      
      //ret[q] +=  - 1./4.*DXQY::c2 * (2*beta)*(2*beta) * (1 - 6*phi0 + 6*phi0*phi0) * cCGNorm[q]; //Beste med ett ekstra ledd

      //ret[q] +=  - 1./2.*DXQY::c2 * (2*beta)*(2*beta) * (1 - 6*phi0 + 6*phi0*phi0) * cCGNorm[q]; 

      //ret[q] +=  - 1./8.*DXQY::c4 *DXQY::c2 * (2*beta)*(2*beta)*(2*beta) * (1 - 14*phi0 + 36*phi0*phi0 - 24*phi0*phi0*phi0) * cCGNorm[q];

      //ret[q] +=  + 0.5*DXQY::c2 * (2*beta)*(2*beta) * (1 - 6*phi0 + 6*phi0*phi0) * cCGNorm[q];

      //ret[q] += 0*rhoFacBeta * 0.5*DXQY::c2Inv* ( cCGNorm[q]*cCGNorm[q] - 0*DXQY::c2 );

      //ret[q] += (2*beta)*(1-2*phi0)*DXQY::c2;

      //ret[q] +=  - 0.5 * DXQY::c4 * (2*beta) * (1-2*phi0) * cCGNorm[q];
      
      //ret[q] *= rhoFacBeta;// * (1 + 0.5*DXQY::c2 * (2*beta)*(2*beta) * ((1-2*phi0)-2*phi0*(2 - 3*phi0)));

      //ret[q] += - 0.5 *(1-2*phi0)*0.5*FNorm0*DXQY::c2;
      
      //ret[q] += DXQY::w[q]*(rho0*cu[q] - cJphi[q])*DXQY::c2Inv;

      ret[q] *= DXQY::w[q];

      //ret[q] += -rhoFacBeta*(2*beta) * (1-2*phi0) *(DXQY::w[q] * cCGNorm[q]*cCGNorm[q] - DXQY::B[q] );
    }
    ret[DXQY::nQNonZero_] = 0.0; // This should be zero by default

    return ret;
}


template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaRC5(const lbBase_t &tauSym, const lbBase_t &tauAnti, const lbBase_t &WInv, const lbBase_t &rhoTot, const lbBase_t &phi0, const lbBase_t &phi1, const std::valarray<lbBase_t> &cCGNorm, const std::valarray<lbBase_t> &cu, const lbBase_t &uCGNorm, const lbBase_t &kx2)
{
    std::valarray<lbBase_t> ret(DXQY::nQ);

    
    //lbBase_t rhoFacWInv = 4* WInv * rhoTot * phi0 * phi1;
    lbBase_t rhoFacWInv = 4* WInv * rhoTot * phi0 * phi1;
    lbBase_t rhoInit = 1.0;
    lbBase_t tauSym_factor = (1 - 0.5 / tauSym);
    lbBase_t tauAnti_factor = (1 - 0.5 / tauAnti);
    lbBase_t A = 1.0;
    
    

    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
      

     
      //ret[q] = tauAnti_factor*cCGNorm[q] + tauSym_factor* (DXQY::c2Inv * cCGNorm[q] * cu[q] - uCGNorm);

      //ret[q] = tauAnti_factor*(DXQY::c2/kx2)*cCGNorm[q]/DXQY::cNorm[q] + tauSym_factor* (DXQY::c2Inv * cCGNorm[q] * cu[q] - uCGNorm); //gives less spurious currents

      ret[q] = tauAnti_factor*(DXQY::c2/kx2)*cCGNorm[q]/DXQY::cNorm[q];
      
      
      //ret[q] -= + 0.5*0.5 / tau * rhoFacWInv * DXQY::c2Inv * (DXQY::cNorm[q]*DXQY::cNorm[q] - DXQY::c2 * (DXQY::nD+2))* cCGNorm[q];

      /*
      ret[q] += - 0.5*0.5 / tau * rhoFacWInv * DXQY::c2Inv * (cCGNorm[q] * cCGNorm[q] * cCGNorm[q] - DXQY::c2 * 3 * cCGNorm[q]);
      ret[q] += + 0.5*0.5 / tau * rhoFacWInv  *(DXQY::cNorm[q]*DXQY::cNorm[q]* cCGNorm[q] - DXQY::c2 * 5 * cCGNorm[q]);
      */
      
      //ret[q] += 0.125*(4*WInv)*(4*WInv)*(1-2*phi0)*(1-phi0)*phi0/(1.1*1.1*(tau-0.5))  *(DXQY::cNorm[q]*DXQY::cNorm[q] - DXQY::nD*DXQY::c2);
      
     
      
    
      ret[q] *= DXQY::w[q]*rhoFacWInv;

     
    }
    ret[DXQY::nQNonZero_] = 0.0; // This should be zero by default

    return ret;
}

#endif // LBCOLLISION2PHASE_H
