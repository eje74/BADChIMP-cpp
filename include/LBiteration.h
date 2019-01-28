#ifndef LBITERATION_H
#define LBITERATION_H

#include "LBglobal.h"
#include "LBbulk.h"
#include "LBgrid.h"
#include "LBfield.h"
#include "LBmacroscopic.h"
#include "LBcollision.h"
#include "LBhalfwaybb.h"

template <typename DXQY>
inline void oneIteration(const lbBase_t &tau, const VectorField<DXQY> &force,
                         const Bulk &bulk, const Grid<DXQY> &grid, const HalfWayBounceBack<DXQY> boundary,
                         ScalarField &rho,  VectorField<DXQY> &vel , LbField<DXQY> &f, LbField<DXQY> &fTmp)
{
    for (int bulkNo = 0; bulkNo < bulk.nElements(); bulkNo++ ) {
        // Find current node number
        const int nodeNo = bulk.nodeNo(bulkNo);
        lbBase_t *forceNode = &force(0, 0, 0); // ASSUME A CONSTANT FORCE

        // UPDATE MACROSCOPIC VARAIBLES
        lbBase_t rhoNode;
        lbBase_t velNode[DXQY::nD];


        calcRho<DXQY>(&f(0,0,nodeNo), rhoNode);
        calcVel<DXQY>(&f(0,0,nodeNo), rhoNode, velNode, forceNode);

        rho(0, nodeNo)  = rhoNode;
        for (int d = 0; d < 2; ++d)
            vel(0, d, nodeNo) = velNode[d];


        // Calculate the BGK collision part
        lbBase_t omegaBGK[DXQY::nQ]; // Holds the collision values
        lbBase_t uu, cu[DXQY::nQ];  // pre-calculated values
        uu = DXQY::dot(velNode, velNode);  // Square of the velocity
        DXQY::cDotAll(velNode, cu);  // velocity dotted with lattice vectors

        calcOmegaBGK<DXQY>(&f(0, 0, nodeNo), tau, rhoNode, uu, cu, omegaBGK);


        lbBase_t deltaOmega[DXQY::nQ];
        lbBase_t  uF, cF[DXQY::nQ];
        uF = DXQY::dot(velNode, forceNode);
        DXQY::cDotAll(forceNode, cF);

        calcDeltaOmega<DXQY>(tau, cu, uF, cF, deltaOmega);

        // * Collision and propagation:
        for (int q = 0; q < DXQY::nQ; q++) {  // Collision should provide the right hand side must be
            fTmp(0, q,  grid.neighbor(q, nodeNo)) = f(0, q, nodeNo) + omegaBGK[q] + deltaOmega[q];
        }


    } // End nodes

    // Swap data_ from fTmp to f;
    f.swapData(fTmp);

    // BOUNDARY CONDITIONS
    boundary.apply(0, f, grid);
}

#endif // LBITERATION_H
