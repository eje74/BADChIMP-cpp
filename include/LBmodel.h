#ifndef LBMODEL_H
#define LBMODEL_H

#include "LBglobal.h"
#include "LBfield.h"
#include "LBgrid.h"

/*
            const int nodeNo = bulk.nodeNo(bulkNo); // Find current node number
            // UPDATE MACROSCOPIC DENSITIES
            lbBase_t rho0Node, rho1Node;
            // Calculate rho for each phase
            calcRho<LT>(&f(0,0,nodeNo), rho0Node);  // LBmacroscopic
            rho(0, nodeNo) = rho0Node; // save to global field
            calcRho<LT>(&f(1,0,nodeNo), rho1Node);  // LBmacroscopic
            rho(1, nodeNo) = rho1Node; // save to global field

            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
*/

/*
    LbField<LT> f(2, nNodes);  // LBfield
    LbField<LT> fTmp(2, nNodes);  // LBfield

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(2, nNodes); // LBfield
    VectorField<LT> vel(1, nNodes); // LBfield
    ScalarField cgField(1, nNodes); // LBfield

 */

template <typename LT>
class TwoPhaseCG
{
public:
    TwoPhaseCG(int nNodes):
        rho(ScalarField(2, nNodes)),
        vel(VectorField<LT>(1, nNodes)),
        f(LbField<LT>(2, nNodes)),
        fTmp(LbField<LT>(2, nNodes)),
        cgField(ScalarField(1, nNodes))
    {
    }

    // FUNCTIONS
    void initRho(Grid<LT>& grid);

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho; //(2, nNodes); // LBfield
    VectorField<LT> vel; //(1, nNodes); // LBfield

private:
    LbField<LT> f; //(2, nNodes);  // LBfield
    LbField<LT> fTmp; //(2, nNodes);  // LBfield
    ScalarField cgField; //(1, nNodes); // LBfield

    lbBase_t rho0Node, rho1Node;
};


template <typename LT>
void TwoPhaseCG<LT>::initRho(Grid<LT> &grid)
{
    for (int n = 1; n < grid.nNodes_; ++n) {
        rho(0, n) = 1.0*(1.0 * std::rand()) / (RAND_MAX * 1.0);
        rho(1, n) = 1 - rho(0, n);
    }
}

#endif // LBMODEL_H
