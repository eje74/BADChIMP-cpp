#ifndef RANS_MAIN_H
#define RANS_MAIN_H

#include "../LBSOLVER.h"
#include "../IO.h"
#include "LBrans.h"
#include <chrono>
#include <numeric>

template <typename DXQY>
void collisionPropagation(int nodeNo, lbBase_t tau, lbBase_t rampup, LbField<DXQY> &f, LbField<DXQY> &fTmp, ScalarField &rho, VectorField<DXQY> &vel, VectorField<DXQY> &bodyForce, Grid<DXQY> &grid)
{

    const std::valarray<lbBase_t> force = rampup * bodyForce(0, nodeNo);
    const std::valarray<lbBase_t> fNode = f(0, nodeNo);

    // Macroscopic values
    //--------------------------------------------------------------------------------- Macroscopic values
    const lbBase_t rhoNode = calcRho<DXQY>(fNode);
    const auto velNode = calcVel<DXQY>(fNode, rhoNode, force);
    rho(0, nodeNo) = rhoNode;
    vel.set(0, nodeNo) = velNode;
    // Collision
    //--------------------------------------------------------------------------------- Collision
    const lbBase_t u2 = DXQY::dot(velNode, velNode);
    const auto cu = DXQY::cDotAll(velNode);
    const auto omegaBGK = calcOmegaBGK<DXQY>(fNode, tau, rhoNode, u2, cu);

    const lbBase_t uF = DXQY::dot(velNode, force);
    const auto cF = DXQY::cDotAll(force);
    const auto deltaOmegaF = calcDeltaOmegaF<DXQY>(tau, cu, uF, cF);

    fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);
}

#endif