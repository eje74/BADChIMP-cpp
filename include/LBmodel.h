#ifndef LBMODEL_H
#define LBMODEL_H

#include <math.h>
#include "LBglobal.h"
#include "LBfield.h"
#include "LBgrid.h"
#include "LBboundary.h"
#include "LBlatticetypes.h"


template <typename LT>
class TwoPhaseCG
{
public:
    TwoPhaseCG(int nNodes):
        rho(ScalarField(2, nNodes)),
        vel(VectorField<LT>(1, nNodes)),
        f(LbField<LT>(2, nNodes)),
        fTmp(LbField<LT>(2, nNodes)),
        bodyForce(VectorField<LT>(1, 1)),
        cgField(ScalarField(1, nNodes))
    {
    }

    // FUNCTIONS
    void setBodyForce(lbBase_t* bF) {for(int d = 0; d < LT::nD; ++d) bodyForce(0, d, 0) = bF[d];}
    void setRhoToConst(const Boundary<LT> &bnd, const lbBase_t rhoConst, const int fieldNo);
    void setTau(const lbBase_t tau0, const lbBase_t tau1) {tau0_ = tau0; tau1_ = tau1;}
    void setNuInv() {nu0Inv_ = 1.0 / (LT::c2 * (tau0_ - 0.5)); nu1Inv_ = 1.0 / (LT::c2 * (tau1_ - 0.5));}
    void setSigma(lbBase_t const sigma) {sigma_ = sigma;}
    void setBeta(lbBase_t const beta) {beta_ = beta;}

    inline LbField<LT>& getF() {return f;}
    inline LbField<LT>& getFtmp() {return fTmp;}

    // -- set local variables
    inline void setLocalFtot(const int nodeNo) { for (int q=0; q < LT::nQ; ++q) {fTotNode[q] = f(0, q, nodeNo) + f(1, q, nodeNo);} }
    // -- fetch local variables
    inline void fetchRho(const int nodeNo) {rho0Node = rho(0, nodeNo); rho1Node = rho(1, nodeNo); rhoNode = rho0Node + rho1Node;}
    // -- init
    void initRho(Grid<LT> &grid);  // Should we just sent a list of node tags. Should we make the node list as a std::vector
    void initVelocity(Grid<LT> &grid);
    void initLbField(Grid<LT> &grid);
    // -- macroscopic
    inline void calcRho(const int &nodeNo);
    inline void calcCgKernel(const int &nodeNo);
    inline void calcForce();
    inline void calcVel(const int &nodeNo);
    // -- mesoscopic
    inline void calcTauEff() {tau_ = LT::c2Inv * rhoNode / (rho0Node * nu0Inv_ + rho1Node * nu1Inv_) + 0.5;}

    // -- collision
    inline void calcOmegaBGK();
    inline void calcDeltaOmega();
    inline void calcDeltaOmegaST(const int &nodeNo, const Grid<LT> &grid);
    // -- propagation
    inline void CollisionPropagation(const int &nodeNo, const Grid<LT> &grid);


    // SETUP MACROSCOPIC FIELDS
    ScalarField rho; //(2, nNodes); // LBfield
    VectorField<LT> vel; //(1, nNodes); // LBfield

private:
    LbField<LT> f; //(2, nNodes);  // LBfield
    LbField<LT> fTmp; //(2, nNodes);  // LBfield
    VectorField<LT> bodyForce;
    ScalarField cgField; //(1, nNodes); // LBfield

    lbBase_t tau_, tau0_, tau1_;
    lbBase_t nu0Inv_, nu1Inv_;
    lbBase_t sigma_;
    lbBase_t beta_;

    lbBase_t cu[LT::nQ], uu;
    lbBase_t cF[LT::nQ], uF;

    lbBase_t rhoNode, rho0Node, rho1Node;
    lbBase_t velNode[LT::nD];
    lbBase_t fTotNode[LT::nQ];
    lbBase_t forceNode[LT::nD];

    // Local collision perturbation
    lbBase_t omegaBGK[LT::nQ]; // Holds the collision values
    lbBase_t deltaOmega[LT::nQ];
    // -- Color Gradient
    lbBase_t deltaOmegaST[LT::nQ];
    lbBase_t bwCos[LT::nQ];

};

template <typename LT>
void TwoPhaseCG<LT>::setRhoToConst(const Boundary<LT> &bnd, const lbBase_t rhoConst, const int fieldNo)
{
    for (int n = 0; n < bnd.getNumNodes(); ++n)
        rho(fieldNo, bnd.nodeNo(n)) = rhoConst;
}


template <typename LT>
void TwoPhaseCG<LT>::initRho(Grid<LT> &grid)
{
    std::srand(8549389);

    for (int n = 1; n < grid.nNodes(); ++n) {
        rho(0, n) = 1.0*(1.0 * std::rand()) / (RAND_MAX * 1.0);
        rho(1, n) = 1 - rho(0, n);
    }
}


template <typename LT>
void TwoPhaseCG<LT>::initVelocity(Grid<LT> &grid)
{
    // Calculate Velocity
    for (int n = 1; n < grid.nNodes(); ++n)
        for (int d = 0; d < LT::nD; ++d)
            vel(0, d, n) = 0.0;

}


template <typename LT>
void TwoPhaseCG<LT>::initLbField(Grid<LT> &grid)
{
    for (int p = 0; p < 2; ++p) { // Phase
        for (int n = 1; n < grid.nNodes(); ++n ) { // Nodes
            // Correct velocity for force
            lbBase_t rhoNode = rho(0, n) + rho(1, n);
            lbBase_t velNode[LT::nD];
            for (int d = 0; d < LT::nD; ++d)
                velNode[d] = vel(0, d, n) - 0.5 * bodyForce(0, d, 0) / (rhoNode + 1.0*(rhoNode < lbBaseEps));

            LT::cDotAll(velNode, cu);
            uu = LT::dot(velNode, velNode);
            for (int q = 0; q < LT::nQ; q++) { // Set the lb field to its equlibrium distribution
                f(p, q, n) = LT::w[q] * rho(p, n) * (1.0 + LT::c2Inv*cu[q] + LT::c4Inv0_5*(cu[q]*cu[q] - LT::c2*uu));
            }
        } // Nodes end
    }  // Phase end
}


template <typename LT>
void inline TwoPhaseCG<LT>::calcRho(const int &nodeNo)
{
    LT::qSum(&f(0, 0, nodeNo), rho(0, nodeNo));
    LT::qSum(&f(1, 0, nodeNo), rho(1, nodeNo));
}


template <typename LT>
void inline TwoPhaseCG<LT>::calcCgKernel(const int &nodeNo)
{
    lbBase_t rho0Node = rho(0, nodeNo);
    lbBase_t rho1Node = rho(1, nodeNo);

    cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
}

template <typename LT>
inline void TwoPhaseCG<LT>::calcForce()
{
    // Calculate local force
    lbBase_t c0 = rho0Node / rhoNode;
    for (int d=0; d < LT::nD; ++d)
        forceNode[d] = c0 * bodyForce(0, d, 0);

    // Calculate force dependent lattice variables
    LT::cDotAll(forceNode, cF);
}


template <typename LT>
inline void TwoPhaseCG<LT>::calcVel(const int &nodeNo)
{
    LT::qSumC(fTotNode, velNode);
    for (int d = 0; d < LT::nD; ++d) {
      velNode[d] = (velNode[d] + 0.5 * forceNode[d]) /rhoNode;
      vel(0, d, nodeNo) = velNode[d];
    }

    // Calcualte velocity dependent lattice variables
    uu = LT::dot(velNode, velNode);  // Square of the velocity
    LT::cDotAll(velNode, cu);  // velocity dotted with lattice vectors

    // Calculate Force velocity product
    uF = LT::dot(velNode, forceNode);
}

template <typename LT>
inline void TwoPhaseCG<LT>::calcOmegaBGK()
/* calcOmegaBGK : sets the BGK-collision term in the lattice boltzmann equation
 * omegaBGK : array of the BGK-collision term in each lattice direction
 */
{
    lbBase_t tau_inv = 1.0 / tau_;
    for (int q = 0; q < LT::nQ; ++q)
    {
        omegaBGK[q] = -tau_inv * ( fTotNode[q] - rhoNode * LT::w[q]*(1.0 + LT::c2Inv*cu[q] + LT::c4Inv0_5*(cu[q]*cu[q] - LT::c2*uu) ) );
    }
}


template <typename LT>
inline void TwoPhaseCG<LT>::calcDeltaOmega()
{
    lbBase_t tau_factor = (1 - 0.5 / tau_);

    for (int q = 0; q < LT::nQ; ++q)
    {
        deltaOmega[q] = LT::w[q]*tau_factor * (LT::c2Inv*cF[q] + LT::c4Inv * ( cF[q] * cu[q] - LT::c2 * uF));
    }
}

template <typename LT>
inline void TwoPhaseCG<LT>::calcDeltaOmegaST(const int &nodeNo, const Grid<LT> &grid)
{
    // -- calculate the normelaized color gradient
    lbBase_t CGNorm, CG2;
    lbBase_t colorGradNode[LT::nD];

    lbBase_t scalarTmp[LT::nQ];
    for (int q = 0; q < LT::nQ; ++q) {
        int neigNode = grid.neighbor(q, nodeNo);
        scalarTmp[q] = cgField(0, neigNode);
    }
    LT::grad(scalarTmp, colorGradNode);
    CG2 = LT::dot(colorGradNode, colorGradNode);
    CGNorm = sqrt(CG2);

    // Normalization of colorGrad
    for (int d = 0; d < LT::nD; ++d)
        colorGradNode[d] /= CGNorm + (CGNorm < lbBaseEps);

    // CALCULATE SURFACE TENSION PERTURBATION
    lbBase_t cCGNorm[LT::nQ];
    LT::cDotAll(colorGradNode, cCGNorm);
    lbBase_t AF0_5 = 2.25 * CGNorm * sigma_ / tau_;
    lbBase_t rhoFacBeta = beta_ * rho0Node * rho1Node / rhoNode;

    for (int q = 0; q < LT::nQNonZero_; ++q) {
        deltaOmegaST[q] = AF0_5 * (LT::w[q] * cCGNorm[q]*cCGNorm[q] - LT::B[q] );
        bwCos[q] = rhoFacBeta * LT::w[q] * cCGNorm[q] /  LT::cNorm[q];
    }
    deltaOmegaST[LT::nQNonZero_] = -AF0_5 * LT::B[LT::nQNonZero_];
    bwCos[LT::nQNonZero_] = 0.0; // This should be zero by default
}


template <typename LT>
inline void TwoPhaseCG<LT>::CollisionPropagation(const int &nodeNo, const Grid<LT> &grid)
{
    lbBase_t c0, c1;
    c0 = (rho0Node/rhoNode);  // Concentration of phase 0
    c1 = (rho1Node/rhoNode);  // Concentration of phase 1
    for (int q = 0; q < LT::nQ; q++) {  // Collision should provide the right hand side must be
        fTmp(0, q,  grid.neighbor(q, nodeNo)) = c0 * (fTotNode[q] + omegaBGK[q] + deltaOmega[q] +  deltaOmegaST[q]) +  bwCos[q];
        fTmp(1, q,  grid.neighbor(q, nodeNo)) = c1 * (fTotNode[q] + omegaBGK[q] + deltaOmega[q] +  deltaOmegaST[q]) -  bwCos[q];
    }
}


#endif // LBMODEL_H
