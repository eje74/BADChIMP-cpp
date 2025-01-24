#ifndef LBCOLLISION_H
#define LBCOLLISION_H

#include "LBglobal.h"
#include "LBlatticetypes.h"


template <typename DXQY>
inline std::valarray<lbBase_t> calcfeq_TEST(const lbBase_t& Gamma0, const lbBase_t& GammaNonZero, const lbBase_t& rho, const lbBase_t& u_sq, const std::valarray<lbBase_t> &cu)
/* calcfeq_TEST : calculates the equilibrium distribution
 *
 * Gamma0       : static mass fraction
 * GammaNonZero : compensation for static mass fraction
 * rho          : density
 * u_sq         : square of the velocity
 * cu           : array of scalar product of all lattice vectors and the velocity.
 * ret          : array of the equilibrium distribution in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    
    for (int q = 0; q < DXQY::nQNonZero_; ++q)
    {
      ret[q] = rho*DXQY::w[q]*(GammaNonZero + (DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq)) ) ;
    }
    ret[DXQY::nQNonZero_] = rho * DXQY::w[DXQY::nQNonZero_]*(Gamma0-0.5*DXQY::c2Inv*u_sq);
    return ret;
}

template <typename DXQY, typename T>
inline std::valarray<lbBase_t> calcOmegaBGK_TEST(const T &f, const T &feq, const lbBase_t &tau)
/* calcOmegaBGK : sets the BGK-collision term in the lattice boltzmann equation
 *
 * f        : pointer to node's lb distribution
 * feq        : pointer to node's lb equilibrium distribution
 * omegaBGK : array of the BGK-collision term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_inv = 1.0 / tau;
    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = -tau_inv * ( f[q] - feq[q]);
    }
    return ret;
}

template <typename DXQY, typename T>
inline std::valarray<lbBase_t> calcOmegaBGK(const T &f, const lbBase_t &tau, const lbBase_t& rho, const lbBase_t& u_sq, const std::valarray<lbBase_t> &cu)
/* calcOmegaBGK : sets the BGK-collision term in the lattice boltzmann equation
 *
 * f        : pointer to node's lb distribution
 * tau      : relaxation time
 * rho      : density
 * u_sq     : square of the velocity
 * cu       : array of scalar product of all lattice vectors and the velocity.
 * omegaBGK : array of the BGK-collision term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_inv = 1.0 / tau;
    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = -tau_inv * ( f[q] - rho * DXQY::w[q]*(1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq) ) );
    }
    return ret;
}



template <typename DXQY, typename T>
inline std::valarray<lbBase_t> calcOmegaBGKTRT(const T &f, const lbBase_t &tauSym, const lbBase_t &tauAnti, const lbBase_t& rho, const lbBase_t& u_sq, const std::valarray<lbBase_t> &cu)
/* calcOmegaBGK : sets the BGK-collision term in the lattice boltzmann equation
 *
 * f        : pointer to node's lb distribution
 * tauSym     : Symmetric relaxation time
 * tauAnti    : Anti-symmetric relaxation time
 * rho      : density
 * u_sq     : square of the velocity
 * cu       : array of scalar product of all lattice vectors and the velocity.
 * omegaBGK : array of the BGK-collision term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    const lbBase_t tauSym_inv = 1.0 / tauSym;
    const lbBase_t tauAnti_inv = 1.0 / tauAnti;
    for (int q = 0; q < DXQY::nQ; ++q)
    {
      const lbBase_t fqSym = 0.5*(f[q]+f[DXQY::reverseDirection(q)]);
      const lbBase_t fqAnti = 0.5*(f[q]-f[DXQY::reverseDirection(q)]);
      
      ret[q] = -tauSym_inv * (fqSym - rho * DXQY::w[q]*(1.0 + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq)))
	-tauAnti_inv * (fqAnti - rho * DXQY::w[q] * DXQY::c2Inv * cu[q]);
    }
    return ret;
}



template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaQ(const lbBase_t &tau, const std::valarray<lbBase_t> &cu, const lbBase_t &u_sq, const lbBase_t &source)
/* calcDeltaOmega : sets the mass source term in the lattice boltzmann equation
 *
 * tau        : relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * u_sq       : square of the velocity
 * source     : source term
 * deltaOmega : array of the force correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_factor = (1 - 0.5 / tau);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = tau_factor * source * DXQY::w[q] * (1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq) );
    }
    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaQTRT(const lbBase_t &tauSym, const lbBase_t &tauAnti, const std::valarray<lbBase_t> &cu, const lbBase_t &u_sq, const lbBase_t &source)
/* calcDeltaOmega : sets the mass source term in the lattice boltzmann equation
 *
 * tauSym     : Symmetric relaxation time
 * tauAnti    : Anti-symmetric relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * u_sq       : square of the velocity
 * source     : source term
 * deltaOmega : array of the force correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tauSym_factor = (1 - 0.5 / tauSym);
    lbBase_t tauAnti_factor = (1 - 0.5 / tauAnti);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
      ret[q] =  source * DXQY::w[q] * ( tauAnti_factor*DXQY::c2Inv*cu[q] + tauSym_factor*(1.0 + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq)) );
    }
    return ret;
}


template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaR(const lbBase_t &tau, const std::valarray<lbBase_t> &cu, const lbBase_t &source)
/* calcDeltaOmega : sets the diffusive mass term in the lattice boltzmann equation
 *
 * tau        : relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * u_sq       : square of the velocity
 * source     : source term
 * deltaOmega : array of the force correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_factor = (1 - 0.5 / tau);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
      ret[q] = tau_factor * source * DXQY::w[q] * (1.0 + DXQY::c2Inv*cu[q]); //ORIGINAL
    }
    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaR2(const lbBase_t &tau, const std::valarray<lbBase_t> &cu, const lbBase_t &source)
/* calcDeltaOmega : sets the diffusive mass term in the lattice boltzmann equation
 *
 * tau        : relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * u_sq       : square of the velocity
 * source     : source term
 * deltaOmega : array of the force correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_factor = (1 - 0.5 / tau);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
      lbBase_t B2q = 0.5 * (DXQY::c2Inv*DXQY::B[q] -DXQY::w[q]);
      //ret[q] = tau_factor * source * DXQY::w[q] * (1.0 + DXQY::c2Inv*cu[q]); //ORIGINAL
      ret[q] = source * tau_factor *  DXQY::w[q];
      //ret[q] = source * tau_factor * ( DXQY::w[q] * (1.0 + DXQY::c2Inv*cu[q]) +  B2q);
    }
    return ret;
}


template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaRTRT(const lbBase_t &tauSym, const lbBase_t &tauAnti, const std::valarray<lbBase_t> &cu, const lbBase_t &source)
/* calcDeltaOmega : sets the diffusive mass term in the lattice boltzmann equation
 *
 * tauSym     : Symmetric relaxation time
 * tauAnti    : Anti-symmetric relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * u_sq       : square of the velocity
 * source     : source term
 * deltaOmega : array of the force correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tauSym_factor = (1 - 0.5 / tauSym);
    lbBase_t tauAnti_factor = (1 - 0.5 / tauAnti);


    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = source * DXQY::w[q] * (1.0*tauSym_factor + tauAnti_factor*DXQY::c2Inv*cu[q]);
    }
    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaF(const lbBase_t &tau, const std::valarray<lbBase_t> &cu, const lbBase_t &uF, const std::valarray<lbBase_t> &cF)
/* calcDeltaOmega : sets the force correction term in the lattice boltzmann equation
 *
 * tau        : relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * uF         : scalar product of velocity and body force.
 * cF         : array of scalar product of all lattice vectors and body force.
 * deltaOmega : array of the force correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_factor = (1 - 0.5 / tau);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = DXQY::w[q]*tau_factor * (DXQY::c2Inv*cF[q] + DXQY::c4Inv * ( cF[q] * cu[q] - DXQY::c2 * uF));
    }
    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaFTRT(const lbBase_t &tauSym, const lbBase_t &tauAnti, const lbBase_t &phi, const std::valarray<lbBase_t> &cu, const lbBase_t &uF, const std::valarray<lbBase_t> &cF)
/* calcDeltaOmega : sets the force correction term in the lattice boltzmann equation
 *
 * tauSym     : Symmetric relaxation time
 * tauAnti    : Anti-symmetric relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * uF         : scalar product of velocity and body force.
 * cF         : array of scalar product of all lattice vectors and body force.
 * deltaOmega : array of the force correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tauSym_factor = (1 - 0.5 / tauSym);
    lbBase_t tauAnti_factor = (1 - 0.5 / tauAnti);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = DXQY::w[q] *phi* (tauAnti_factor*DXQY::c2Inv*cF[q] + tauSym_factor*DXQY::c4Inv * ( cF[q] * cu[q] - DXQY::c2 * uF));
    }
    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaFDiff(const lbBase_t &tau, const lbBase_t &phi, const std::valarray<lbBase_t> &cu, const lbBase_t &uF, const std::valarray<lbBase_t> &cF)
/* calcDeltaOmega : sets the force correction term in the lattice boltzmann equation
 *
 * tau        : relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * uF         : scalar product of velocity and body force.
 * cF         : array of scalar product of all lattice vectors and body force.
 * deltaOmega : array of the force correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_factor = (1 - 0.5 / tau);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = DXQY::w[q]*tau_factor * (DXQY::c2Inv*cF[q]*phi);
    }
    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaFDiffTRT(const lbBase_t &tauSym, const lbBase_t &tauAnti, const lbBase_t &phi, const std::valarray<lbBase_t> &cu, const lbBase_t &uF, const std::valarray<lbBase_t> &cF)
/* calcDeltaOmega : sets the force correction term in the lattice boltzmann equation
 *
 * tauSym     : Symmetric relaxation time
 * tauAnti    : Anti-symmetric relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * uF         : scalar product of velocity and body force.
 * cF         : array of scalar product of all lattice vectors and body force.
 * deltaOmega : array of the force correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tauAnti_factor = (1 - 0.5 / tauAnti);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = DXQY::w[q]*tauAnti_factor * (DXQY::c2Inv*cF[q]*phi);
    }
    return ret;
}


template <typename DXQY, typename T>
inline std::valarray<lbBase_t> calcOmegaPureDiff(const T &f, const lbBase_t &tau, const lbBase_t &rho)
/* calcOmegaBGK : sets the BGK-collision term in the lattice boltzmann equation
 *
 * f        : pointer to node's lb distribution
 * tau      : relaxation time
 * rho      : density
 * omegaBGK : array of the BGK-collision term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_inv = 1.0 / tau;
    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = -tau_inv * ( f[q] - rho * DXQY::w[q] ) ;
    }
    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaGradDiff(const lbBase_t &tau, const std::valarray<lbBase_t> &cu, const lbBase_t &uGrad, const std::valarray<lbBase_t> &cGrad)
/* calcDeltaOmega : sets the force correction term in the lattice boltzmann equation
 *
 * tau        : relaxation time
 * cu         : array of scalar product of all lattice vectors and the velocity.
 * uGrad         : scalar product of velocity and grad phi.
 * cGrad         : array of scalar product of all lattice vectors and grad phi.
 * deltaOmega : array of the grad phi correction term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t tau_factor = (1 - 0.5 / tau);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = DXQY::w[q]*tau_factor * DXQY::c2Inv * ( cGrad[q] * cu[q] - DXQY::c2 * uGrad);
    }
    return ret;
}


#endif // LBCOLLISION_H
