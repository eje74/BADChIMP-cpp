#ifndef LBCOLLISION_H
#define LBCOLLISION_H

#include "LBglobal.h"
#include "LBlatticetypes.h"


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


template <typename DXQY>
inline std::valarray<lbBase_t> calcDeltaOmegaQ(const lbBase_t &tau, const std::valarray<lbBase_t> &cu, const lbBase_t &u_sq, const lbBase_t &source)
/* calcDeltaOmega : sets the force correction term in the lattice boltzmann equation
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

#endif // LBCOLLISION_H
