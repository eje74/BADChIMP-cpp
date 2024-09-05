#ifndef LBBOUNDARYRANS
#define LBBOUNDARYRANS

#include "../LBSOLVER.h"
#include "../IO.h"
#include "LBBoundaryInterpolation.h"
#include <numeric>

template <typename DXQY>
inline std::valarray<lbBase_t> calcfeq(const lbBase_t& rho, const std::valarray<lbBase_t> &velNode);

template <typename DXQY>
std::valarray<lbBase_t> calcfneq(const lbBase_t rhoNode, const std::valarray<lbBase_t> &velNode, const std::valarray<lbBase_t> &fNode);

template <typename DXQY>
std::valarray<lbBase_t> calcfneqWall(const lbBase_t tau, const lbBase_t rho, const lbBase_t &dut_dn, const std::valarray<lbBase_t> &F, const InterpolationElement &ie);

template <typename DXQY>
std::valarray<lbBase_t> calcDeltafneqwall(const std::valarray<lbBase_t> &fneq, const InterpolationElement &ie);

lbBase_t interpolateScalar(const ScalarField &rho, const std::vector<lbBase_t> &w, const std::vector<int> & pnts);

template<typename DXQY>
std::valarray<lbBase_t> interpolateVector(const VectorField<DXQY> &val, const std::vector<lbBase_t> &w, const std::vector<int> &pnts);

template<typename DXQY>
std::valarray<lbBase_t> interpolateLBFieldNeq(const ScalarField &rho, const VectorField<DXQY> &vel, const LbField<DXQY> &f, const std::vector<lbBase_t> &w, const std::vector<int> &pnts);

//========================================================================================================== wall_boundary_rans works for newtonian
template<typename DXQY>
void wall_boundary_rans_old(
    lbBase_t tau, 
    lbBase_t rampup, 
    LbField<DXQY> &f, 
    LbField<DXQY> &fTmp,
    ScalarField &rho, 
    VectorField<DXQY> &vel,
    VectorField<DXQY> &bodyForce, 
    std::vector<InterpolationElement> &boundarynodes, 
    Grid<DXQY> &grid)
{
    for (auto &bn: boundarynodes) {
        const int nodeNo = bn.nodeNo;

        const lbBase_t rhoNode = 1.0; //interpolateScalar(rho, bn.wb, bn.pnts);
        const std::valarray<lbBase_t> velNodeIntrp = interpolateVector(vel, bn.wb, bn.pnts);
        const std::valarray<lbBase_t> velNode = bn.gamma*velNodeIntrp/(bn.gamma + bn.gamma2);
        const std::valarray<lbBase_t> fneq = interpolateLBFieldNeq(rho, vel, f, bn.wb, bn.pnts);
        const std::valarray<lbBase_t> tvec = {-bn.normal[1], bn.normal[0]};
        const lbBase_t dvel_dn = DXQY::dot(velNodeIntrp, tvec)/(bn.gamma + bn.gamma2);
        const std::valarray<lbBase_t> forceNode = rampup*bodyForce(0, nodeNo);
        const std::valarray<lbBase_t> fneq_wall = calcfneqWall<DXQY>(tau, rhoNode, dvel_dn, forceNode, bn);

        const std::valarray<lbBase_t> delta_fneq = bn.gamma2*calcDeltafneqwall<DXQY>(fneq, bn)/(bn.gamma + bn.gamma2);

        f.set(0, nodeNo) = calcfeq<DXQY>(rhoNode, velNode) + (bn.gamma*fneq + bn.gamma2*fneq_wall)/(bn.gamma + bn.gamma2);

        collisionPropagation(nodeNo, tau, rampup, f, fTmp, rho, vel, bodyForce, grid);

        rho(0, nodeNo) = interpolateScalar(rho, bn.wb, bn.pnts);


    }
}

//========================================================================================================== wall_boundary_rans

template<typename DXQY>
void wall_boundary_rans_new(
    lbBase_t tau, 
    lbBase_t rampup, 
    LbField<DXQY> &f, 
    LbField<DXQY> &fTmp,
    ScalarField &rho, 
    VectorField<DXQY> &vel,
    VectorField<DXQY> &bodyForce, 
    std::vector<InterpolationElement> &boundarynodes, 
    Grid<DXQY> &grid)
{
    for (auto &bn: boundarynodes) {
        const int nodeNo = bn.nodeNo;

        const lbBase_t rhoNode = 1.0; //interpolateScalar(rho, bn.wb, bn.pnts);
        const std::valarray<lbBase_t> velNodeIntrp = interpolateVector(vel, bn.wb, bn.pnts);
        std::valarray<lbBase_t> velWall = interpolateVector(vel, bn.wa, bn.pnts);
        const std::valarray<lbBase_t> tvec = {-bn.normal[1], bn.normal[0]};

        velWall = DXQY::dot(velWall, tvec)*tvec;

        if (DXQY::dot(velWall, tvec)*DXQY::dot(velNodeIntrp, tvec) < 0) {
                velWall = 0*DXQY::dot(velWall, tvec)*tvec;

        }  

        if (  std::abs(DXQY::dot(velWall, tvec)) >  std::abs(DXQY::dot(velNodeIntrp, tvec)) ) {
            velWall = DXQY::dot(velNodeIntrp, tvec)*tvec;
        }

        const std::valarray<lbBase_t> velNode = (bn.gamma*velNodeIntrp + bn.gamma2*velWall)/(bn.gamma + bn.gamma2);

        //const std::valarray<lbBase_t> velNode = bn.gamma*velNodeIntrp/(bn.gamma + bn.gamma2);
        const std::valarray<lbBase_t> fneq = interpolateLBFieldNeq(rho, vel, f, bn.wb, bn.pnts);
        const lbBase_t dvel_dn = DXQY::dot(velNodeIntrp, tvec)/(bn.gamma + bn.gamma2);
        const std::valarray<lbBase_t> forceNode = rampup*bodyForce(0, nodeNo);
        const std::valarray<lbBase_t> fneq_wall = calcfneqWall<DXQY>(tau, rhoNode, dvel_dn, forceNode, bn);

        const std::valarray<lbBase_t> delta_fneq = bn.gamma2*calcDeltafneqwall<DXQY>(fneq, bn)/(bn.gamma + bn.gamma2);

        f.set(0, nodeNo) = calcfeq<DXQY>(rhoNode, velNode) - 3*0.5*D2Q9::cDotAll(forceNode) + 0*(bn.gamma*fneq + bn.gamma2*fneq_wall)/(bn.gamma + bn.gamma2);

        collisionPropagation(nodeNo, tau, rampup, f, fTmp, rho, vel, bodyForce, grid);

        rho(0, nodeNo) = interpolateScalar(rho, bn.wb, bn.pnts);


    }
}


template<typename DXQY>
void wall_boundary_rans(
    lbBase_t tau, 
    lbBase_t rampup, 
    LbField<DXQY> &f, 
    LbField<DXQY> &fTmp,
    ScalarField &rho, 
    VectorField<DXQY> &vel,
    VectorField<DXQY> &bodyForce, 
    std::vector<InterpolationElement> &boundarynodes, 
    Grid<DXQY> &grid)
{
    for (auto &bn: boundarynodes) {
        const int nodeNo = bn.nodeNo;

        const lbBase_t rhoNode = 1.0; //interpolateScalar(rho, bn.wb, bn.pnts);
        const std::valarray<lbBase_t> velNodeIntrp = interpolateVector(vel, bn.wb, bn.pnts);
        std::valarray<lbBase_t> velWall = interpolateVector(vel, bn.wa, bn.pnts);
        const std::valarray<lbBase_t> tvec = {-bn.normal[1], bn.normal[0]};

        velWall = DXQY::dot(velWall, tvec)*tvec;

        if (DXQY::dot(velWall, tvec)*DXQY::dot(velNodeIntrp, tvec) < 0) {
                velWall = 0*DXQY::dot(velWall, tvec)*tvec;

        }  

        if (  std::abs(DXQY::dot(velWall, tvec)) >  std::abs(DXQY::dot(velNodeIntrp, tvec)) ) {
            velWall = DXQY::dot(velNodeIntrp, tvec)*tvec;
        }

        const std::valarray<lbBase_t> velNode = (bn.gamma*velNodeIntrp + bn.gamma2*velWall)/(bn.gamma + bn.gamma2);

        //const std::valarray<lbBase_t> velNode = bn.gamma*velNodeIntrp/(bn.gamma + bn.gamma2);
        const std::valarray<lbBase_t> fneq = interpolateLBFieldNeq(rho, vel, f, bn.wb, bn.pnts);
        const lbBase_t dvel_dn = DXQY::dot(velNodeIntrp, tvec)/(bn.gamma + bn.gamma2);
        const std::valarray<lbBase_t> forceNode = rampup*bodyForce(0, nodeNo);
        const std::valarray<lbBase_t> fneq_wall = calcfneqWall<DXQY>(tau, rhoNode, dvel_dn, forceNode, bn);

        const std::valarray<lbBase_t> delta_fneq = bn.gamma2*calcDeltafneqwall<DXQY>(fneq, bn)/(bn.gamma + bn.gamma2);

        f.set(0, nodeNo) = calcfeq<DXQY>(rhoNode, velNode) - 3*0.5*D2Q9::cDotAll(forceNode) + 0*(bn.gamma*fneq + bn.gamma2*fneq_wall)/(bn.gamma + bn.gamma2);

        collisionPropagation(nodeNo, tau, rampup, f, fTmp, rho, vel, bodyForce, grid);

        rho(0, nodeNo) = interpolateScalar(rho, bn.wb, bn.pnts);


    }
}


template <typename DXQY>
inline std::valarray<lbBase_t> calcfeq(const lbBase_t& rho, const std::valarray<lbBase_t> &velNode)
/* 
 * rho      : density
 * u_sq     : square of the velocity
 * cu       : array of scalar product of all lattice vectors and the velocity.
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);

    const lbBase_t u_sq = DXQY::dot(velNode, velNode);
    const std::valarray<lbBase_t> cu = DXQY::cDotAll(velNode);

    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = rho * DXQY::w[q]*(1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq) );
    }
    return ret;
}

template <typename DXQY>
std::valarray<lbBase_t> calcfneq(const lbBase_t rhoNode, const std::valarray<lbBase_t> &velNode, const std::valarray<lbBase_t> &fNode)
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    ret = fNode - calcfeq<DXQY>(rhoNode, velNode);
    return ret;
}

template <typename DXQY>
std::valarray<lbBase_t> calcfneqWall(const lbBase_t tau, const lbBase_t rho, const lbBase_t &dut_dn, const std::valarray<lbBase_t> &F, const InterpolationElement &ie)
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    const std::vector nvec = ie.normal;
    const std::vector tvec = {-nvec[1], nvec[0]};

    const auto cn = DXQY::cDotAll(nvec);
    const auto ct = DXQY::cDotAll(tvec);
    const auto cF = DXQY::cDotAll(F);

    for (int q = 0; q < DXQY::nQ; ++q) {
        const lbBase_t Qnt =  cn[q]*ct[q];
        ret[q] = -DXQY::w[q]*DXQY::c2Inv*(tau*rho*dut_dn*Qnt + cF[q]);
    }
    return ret;
}


template <typename DXQY>
std::valarray<lbBase_t> calcDeltafneqwall(const std::valarray<lbBase_t> &fneq, const InterpolationElement &ie)
{
    std::valarray<lbBase_t> ret(0.0, DXQY::nQ);
    const std::vector<lbBase_t> nvec = ie.normal;
    const std::vector<lbBase_t> tvec = {-nvec[1], nvec[0]};

    // Calculate Pi_nn and Pi_tt
    const auto cn = DXQY::cDotAll(nvec);
    const auto ct = DXQY::cDotAll(tvec);

    lbBase_t Pi_nn = 0.0;
    lbBase_t Pi_tt = 0.0;

    for (int q = 0; q < DXQY::nQ; ++q) {
        Pi_nn += (cn[q]*cn[q] - DXQY::c2)*fneq[q];
        Pi_tt += (ct[q]*ct[q] - DXQY::c2)*fneq[q];
    }

    // Calculate the correction to fneq at the wall
    for (int q = 0; q < DXQY::nQ; ++q) {
        ret[q] = -DXQY::c4Inv0_5*DXQY::w[q]*((cn[q]*cn[q] - DXQY::c2)*Pi_nn + (ct[q]*ct[q] - DXQY::c2)*Pi_tt);
    }

    return ret;
}



lbBase_t interpolateScalar(const ScalarField &rho, const std::vector<lbBase_t> &w, const std::vector<int> & pnts)
{
    lbBase_t ret = 0.0;
        for (std::size_t i = 0; i < pnts.size(); ++i) {
            const int ni = pnts[i];
            ret += w[i]*rho(0, ni);
        }
    return ret;
}

template<typename DXQY>
std::valarray<lbBase_t> interpolateVector(const VectorField<DXQY> &val, const std::vector<lbBase_t> &w, const std::vector<int> &pnts)
{
    std::valarray<lbBase_t> ret(0.0, DXQY::nD);

        for (std::size_t i = 0; i < pnts.size(); ++i) {
            const int ni = pnts[i];
            ret += w[i]*val(0, ni);
        }
    return ret;
}

template<typename DXQY>
std::valarray<lbBase_t> interpolateLBFieldNeq(const ScalarField &rho, const VectorField<DXQY> &vel, const LbField<DXQY> &f, const std::vector<lbBase_t> &w, const std::vector<int> &pnts)
{
    std::valarray<lbBase_t> ret(0.0, DXQY::nQ);

    for (std::size_t i = 0; i < pnts.size(); ++i)
    {
        const int n = pnts[i];
        const lbBase_t rhoNode = rho(0, n);
        const std::valarray<lbBase_t> velNode = vel(0, n);
        const std::valarray<lbBase_t> fNode = f(0, n);

        ret += w[i]*(fNode - calcfeq<DXQY>(rhoNode, velNode));
    }

    return ret;
}

#endif