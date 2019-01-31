// //////////////////////////////////////////////
//
// TWO PHASE SOLVER
//
//
// Using the color gradient method (Rothman-Keller
//  type) with the Reis Phillips surface pertubation
//  and / Latva-Kokko recolor step. (As suggested
//  by Leclaire et al and Yu-Hang Fu et al)
//
// This is an explicit two-phase implementation. A
// full multiphase implementation will be added later.
//
// TO RUN PROGRAM: type bin/runner on command
// line in main directory
//
// //////////////////////////////////////////////

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "LBlatticetypes.h"
#include "LBgrid.h"
#include "LBfield.h"
#include "LBboundary.h"
#include "LBhalfwaybb.h"
#include "LBbulk.h"
#include "LBgeometry.h"
#include "LBmacroscopic.h"
#include "LBcollision.h"
#include "LBiteration.h"
#include "LBinitiatefield.h"


// SET THE LATTICE TYPE
#define LT D2Q9

int main()
{
    std::cout << "Begin test Two phase";
    std::cout << std::endl;


    // INPUT DATA
    int nIterations;
    int nX, nY;
    lbBase_t tau0, tau1;
    lbBase_t nu0Inv, nu1Inv;

    lbBase_t force[2] = {1.0e-8, 0.0};
    VectorField<LT> bodyForce(1,1);
    bodyForce(0, 0, 0) = force[0];
    bodyForce(0, 1, 0) = force[1];

    nIterations = 10000;
    nX = 250; nY = 101;
    //nX = 1; nY = 101;

    tau0 = 0.7;
    tau1 = 0.7;
    nu0Inv = 1.0 / (LT::c2 * (tau0 - 0.5));
    nu1Inv = 1.0 / (LT::c2 * (tau1 - 0.5));

    // SETUP GEOMETRY
    int ** geo;
    newGeometry(nX, nY, geo); // LBgeometry
    inputGeometry(nX, nY, geo); // LBgeometry
    analyseGeometry<LT>(nX, nY, geo); // LBgeometry

    int ** labels;
    newNodeLabel(nX, nY, labels); // LBgeometry

    int nBulkNodes = 0;
    nBulkNodes = setBulkLabel(nX, nY, geo, labels); // LBgeometry
    int nNodes = 0;
    nNodes = setNonBulkLabel(nBulkNodes, nX, nY, geo, labels); // LBgeometry

    // USE NODENO 0 AS DUMMY NODE.
    // Add one extra storage for the default dummy node
    nNodes += 1;

    // SETUP GRID
    std::cout << "NUMBER OF NODES = " << nNodes << std::endl;
    Grid<LT> grid(nNodes); // object declaration: neigList_ stores node numbers of neighbors;  pos_ stores Cartesian coordinates of given node
    setupGrid(nX, nY, labels, grid); // LBgeometry, maybe move to LBgrid?

    // SETUP BULK
    std::cout << "NUMBER OF Bulk NODES = " << nBulkNodes << std::endl;
    Bulk bulk(nBulkNodes); // object declaration: nBulkNodes_ stores number of Bulk nodes; bulkNode_ stores node numbers of all bulk nodes
    setupBulk(nX, nY, geo, labels, bulk); // LBgeometry

    // SETUP BOUNDARY
    int nBoundary = nBoundaryNodes(1, nX, nY, geo);
    std::cout << "NUMBER OF BOUNDARY NODES = " << nBoundary << std::endl;
    HalfWayBounceBack<LT> boundary( nBoundary ); // LBhalfwaybb
    setupBoundary(1, nX, nY, geo, labels, grid, boundary); // LBgeometry

    // SETUP LB FIELDS
    LbField<LT> f(2, nNodes);  // LBfield
    LbField<LT> fTmp(2, nNodes);  // LBfield

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(2, nNodes); // LBfield
    VectorField<LT> vel(1, nNodes); // LBfield

    // FILL MACROSCOPIC FIELDS
    // -- Phase 0
    setFieldToConst(0.5, 0, rho); // LBmacroscopic
    // -- Phase 1
    setFieldToConst(0.5, 1, rho); // LBmacroscopic
    // -- Phase total
    lbBase_t velTmp[LT::nD] = {0.0, 0.0};
    setFieldToConst(velTmp, 0, vel);  // LBmacroscopic


    // INITIATE LB FIELDS
    // -- phase 0
    initiateLbField(0, 0, 0, bulk, rho, vel, f);  // LBinitiatefield
    // -- phase 1
    initiateLbField(1, 1, 0, bulk, rho, vel, f);  // LBinitiatefield

    // -----------------MAIN LOOP------------------
    /* Comments to main loop:
     * Calculation of cu is kept outside of calcOmega and calcDeltaOmega
     *  to avoid recalculation  of the dot-product.
     * It is therefore natural to give all scalar products as inputs to the
     *  calcOmega and calcDeltaOmega.
     * Note that we have a oneIteration function in 'LBiteration.h'. This seemed
     *   to affect the speed an gave a slow down of 10-15 %.
     */

    for (int i = 0; i < nIterations; i++) {
        for (int bulkNo = 0; bulkNo < bulk.nElements(); bulkNo++ ) {
            const int nodeNo = bulk.nodeNo(bulkNo); // Find current node number

            // UPDATE MACROSCOPIC VARIABLES
            lbBase_t rhoNode, rho0Node, rho1Node;
            lbBase_t velNode[LT::nD];
            lbBase_t *forceNode = &bodyForce(0,0,0);
            lbBase_t fTot[LT::nQ];

            // Set the local total lb distribution
            for (int q = 0; q < LT::nQ; ++q)
                fTot[q] = f(0, q, nodeNo) + f(1, q, nodeNo);

            // Calculate rho for each phase
            calcRho<LT>(&f(0,0,nodeNo), rho0Node);  // LBmacroscopic
            rho(0, nodeNo) = rho0Node; // save to global field
            calcRho<LT>(&f(1,0,nodeNo), rho1Node);  // LBmacroscopic
            rho(1, nodeNo) = rho1Node; // save to global field
            // Total density
            rhoNode = rho0Node + rho1Node;
            // Total velocity
            calcVel<LT>(fTot, rhoNode, velNode, forceNode);  // LBmacroscopic
            vel(0, 0, nodeNo) = velNode[0];
            vel(0, 1, nodeNo) = velNode[1];

            // CALCULATE BGK COLLISION TERM
            // Mean collision time /rho_tot/\nu_tot = \sum_s \rho_s/\nu_s
            lbBase_t tau;
            tau = LT::c2Inv * rhoNode / (rho0Node*nu0Inv + rho1Node*nu1Inv) + 0.5;

            lbBase_t omegaBGK[LT::nQ]; // Holds the collision values
            lbBase_t uu, cu[LT::nQ];  // pre-calculated values
            uu = LT::dot(velNode, velNode);  // Square of the velocity
            LT::cDotAll(velNode, cu);  // velocity dotted with lattice vectors
            calcOmegaBGK<LT>(fTot, tau, rhoNode, uu, cu, omegaBGK);  // LBcollision

            // CALCULATE FORCE CORRECTION TERM
            lbBase_t deltaOmega[LT::nQ];
            lbBase_t  uF, cF[LT::nQ];
            uF = LT::dot(velNode, forceNode);
            LT::cDotAll(forceNode, cF);

            calcDeltaOmega<LT>(tau, cu, uF, cF, deltaOmega);  // LBcollision

            // COLLISION AND PROPAGATION
            lbBase_t c0, c1;
            c0 = (rho0Node/rhoNode);  // Concentration of phase 0
            c1 = (rho1Node/rhoNode);  // Concentration of phase 1
            for (int q = 0; q < LT::nQ; q++) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = c0 * (fTot[q] + omegaBGK[q] + deltaOmega[q]);
                fTmp(1, q,  grid.neighbor(q, nodeNo)) = c1 * (fTot[q] + omegaBGK[q] + deltaOmega[q]);
            }
        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // BOUNDARY CONDITIONS
        boundary.apply(0, f, grid);  // LBboundary
        boundary.apply(1, f, grid);
    } // End iterations
    // -----------------END MAIN LOOP------------------
    
    std::cout << std::setprecision(5) << vel(0, 1, labels[1][0])  << " " << vel(0, 0, labels[50][0])  << std::endl;


    // CLEANUP
    deleteNodeLabel(nX, nY, labels);
    deleteGeometry(nX, nY, geo);

    return 0;
}

