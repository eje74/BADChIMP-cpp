#ifndef LBMACROSCOPIC_H
#define LBMACROSCOPIC_H

#include "LBglobal.h"
#include "LBlatticetypes.h"

template <typename DXQY>
inline void calcRho(const lbBase_t* f, lbBase_t &rho)
{
     DXQY::qSum(f, rho);
}

template <typename DXQY>
inline void calcVel(const lbBase_t* f, const lbBase_t &rho, lbBase_t* vel, const lbBase_t* F)
{
    DXQY::qSumC(f, vel);
    for (int d = 0; d < DXQY::nD; ++d) {
      vel[d] = (vel[d] + 0.5 * F[d]) /rho;
    }

}

#endif // LBMACROSCOPIC_H
