#ifndef LBSNIPPETS_H
#define LBSNIPPETS_H

#include "LBfield.h"
#include "LBglobal.h"

// Constant pressure with oil
inline void setConstDensity(lbBase_t &q0, lbBase_t &q1, const lbBase_t rho0, const lbBase_t rho1, const lbBase_t rhoConst)
// Assuming index 0: water and 1:oil
// Pressure is kept constant by adding oil and removing oil/water mixture
{

    lbBase_t rhoTot = rho0 + rho1;
    lbBase_t qTmp = 2*(rhoConst - rhoTot);

    lbBase_t x1 = 1.0  - (rho0/rhoTot)*(qTmp < 0.0);

    q1 = x1 * qTmp;
    q0 = qTmp - q1;
}

inline void setConstSource(lbBase_t &q0, lbBase_t &q1, const lbBase_t rho0, const lbBase_t rho1, const lbBase_t rate)
// Assuming index 0: water and 1:oil
// adding oil/water mixture at constant rate
{
    q0 = rate * rho0/(rho0 + rho1);
    q1 = rate - q0;
}


void printAsciiToScreen(int nX, int nY, ScalarField& val, int** labels);

void printAsciiToScreen(int nX, int nY, ScalarField& val, int** lab, double pauseInSeconds);

void printAsciiToScreen(int nX, int nY, int nZ, int itr, ScalarField& val, int*** lab, double pauseInSeconds);


#endif // LBSNIPPETS_H
