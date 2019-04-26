#ifndef LBSNIPPETS_H
#define LBSNIPPETS_H

#include "LBfield.h"


// Constant pressure with oil
void setConstPressure(const lbBase_t rho0, const lbBase_t rho1, const lbBase_t rhoConst);

void printAsciiToScreen(int nX, int nY, ScalarField& val, int** labels);


void printAsciiToScreen(int nX, int nY, ScalarField& val, int** lab, double pauseInSeconds);

void printAsciiToScreen(int nX, int nY, int nZ, int itr, ScalarField& val, int*** lab, double pauseInSeconds);

#endif // LBSNIPPETS_H
