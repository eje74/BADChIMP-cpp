#ifndef LBSNIPPETS_H
#define LBSNIPPETS_H

#include "LBfield.h"

void printAsciiToScreen(int nX, int nY, ScalarField& val, int** labels);


void printAsciiToScreen(int nX, int nY, ScalarField& val, int** lab, double pauseInSeconds);

void printAsciiToScreen(int nX, int nY, int nZ, int itr, ScalarField& val, int*** lab, double pauseInSeconds);

#endif // LBSNIPPETS_H
