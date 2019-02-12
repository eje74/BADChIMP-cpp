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
#include <math.h>
#include <stdlib.h>
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
#include "LBsnippets.h"


// SET THE LATTICE TYPE
//#define LT D2Q9
#define LT D3Q19


int main()
{
    std::cout << "Begin test Two phase";
    std::cout << std::endl;


    // INPUT DATA
    int nIterations;
    int nX, nY;
    int nZ;
    lbBase_t tau0, tau1;
    lbBase_t nu0Inv, nu1Inv;

    lbBase_t beta;
    lbBase_t sigma;

    //lbBase_t force[2] = {1.0e-5, 1.e-4};// {1.0e-8, 0.0};
    lbBase_t force[3] = {0.0, 0.0, 0.0}; //{1.0e-5, 1.e-4, 0.0};
    VectorField<LT> bodyForce(1,1);
    bodyForce(0, 0, 0) = force[0];
    bodyForce(0, 1, 0) = force[1];
    bodyForce(0, 2, 0) = force[2];
    
    nIterations = 100000;//100000;
    //nX = 250; nY = 101;
    //nX = 190; nY = 120;
    nX = 140;
    nY = 101;
    nZ = 2;

    
    tau0 = 1.0;
    tau1 = 1.0;


    sigma = 0.0001;

    beta = 1.0;

    nu0Inv = 1.0 / (LT::c2 * (tau0 - 0.5));
    nu1Inv = 1.0 / (LT::c2 * (tau1 - 0.5));

    // SETUP GEOMETRY
    /*
    int ** geo;
    newGeometry(nX, nY, geo); // LBgeometry
    inputGeometry(nX, nY, geo); // LBgeometry
    analyseGeometry<LT>(nX, nY, geo); // LBgeometry
    */
    int *** geo;
    newGeometry(nX, nY, nZ, geo); // LBgeometry
    inputGeometry(nX, nY, nZ, geo); // LBgeometry
    analyseGeometry<LT>(nX, nY, nZ, geo); // LBgeometry
    
    /*
    int ** labels;
    newNodeLabel(nX, nY, labels); // LBgeometry
    */
    int *** labels;
    newNodeLabel(nX, nY, nZ, labels); // LBgeometry
    

    int nBulkNodes = 0;
    // int bulkLabel[] = {0, 1, 2};
    // Note to Self: add an list of input that defines what
    //  values in geo that are treated as bulk

    //nBulkNodes = setBulkLabel(nX, nY, geo, labels); // LBgeometry
    nBulkNodes = setBulkLabel(nX, nY, nZ, geo, labels); // LBgeometry
    
    int nNodes = 0;
    // Note to Self: add an list of input that defines what
    //  values in geo that are treated as nonbulk

    //nNodes = setNonBulkLabel(nBulkNodes, nX, nY, geo, labels); // LBgeometry
    nNodes = setNonBulkLabel(nBulkNodes, nX, nY, nZ, geo, labels); // LBgeometry
    
    // USE NODENO 0 AS DUMMY NODE.
    // Add one extra storage for the default dummy node
    nNodes += 1;

    // SETUP GRID
    std::cout << "NUMBER OF NODES = " << nNodes << std::endl;
    Grid<LT> grid(nNodes); // object declaration: neigList_ stores node numbers of neighbors;  pos_ stores Cartesian coordinates of given node

    //setupGrid(nX, nY, labels, grid); // LBgeometry, maybe move to LBgrid?
    setupGrid(nX, nY, nZ, labels, grid); // LBgeometry, maybe move to LBgrid?
	
    
    // SETUP BULK
    std::cout << "NUMBER OF Bulk NODES = " << nBulkNodes << std::endl;
    Bulk bulk(nBulkNodes); // object declaration: nBulkNodes_ stores number of Bulk nodes; bulkNode_ stores node numbers of all bulk nodes

    //setupBulk(nX, nY, geo, labels, bulk); // LBgeometry
    setupBulk(nX, nY, nZ, geo, labels, bulk); // LBgeometry

    
    // SETUP BOUNDARY
    int nBoundary = 0;
    
    //nBoundary = nBoundaryNodes(1, nX, nY, geo);
    nBoundary = nBoundaryNodes(1, nX, nY, nZ, geo);

    std::cout << "NUMBER OF BOUNDARY NODES = " << nBoundary << std::endl;
    HalfWayBounceBack<LT> boundary( nBoundary ); // LBhalfwaybb
    
    //setupBoundary(1, nX, nY, geo, labels, grid, boundary); // LBgeometry
    setupBoundary(1, nX, nY, nZ, geo, labels, grid, boundary); // LBgeometry

    
    // SETUP SOLID BOUNDARY
    int nSolidBoundary = 0;

    //nSolidBoundary = nBoundaryNodes(3, nX, nY, geo);
    nSolidBoundary = nBoundaryNodes(3, nX, nY, nZ, geo);

    
    Boundary<LT> solidBoundary( nSolidBoundary );
    //setupBoundary(3, nX, nY, geo, labels, grid, solidBoundary); // LBgeometry
    setupBoundary(3, nX, nY, nZ, geo, labels, grid, solidBoundary); // LBgeometry


    // SETUP LB FIELDS
    LbField<LT> f(2, nNodes);  // LBfield
    LbField<LT> fTmp(2, nNodes);  // LBfield

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(2, nNodes); // LBfield
    VectorField<LT> vel(1, nNodes); // LBfield
    VectorField<LT> colorGrad(1, nNodes); // LBfield


    // FILL MACROSCOPIC FIELDS
    // -- Phase 0
    setFieldToConst(0.6, 0, rho); // LBmacroscopic
    // -- Phase 1
    setFieldToConst(0.4, 1, rho); // LBmacroscopic


    std::srand(8549389);
    /*
    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            rho(0, labels[y][x]) = 1.0*(1.0 * std::rand()) / (RAND_MAX * 1.0);
            rho(1, labels[y][x]) = 1 - rho(0, labels[y][x]);
        }
    }
    */
    
    for (int z = 0; z < nZ; ++z) {
        for (int y = 0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x) {
                rho(0, labels[z][y][x]) = 1.0*(1.0 * std::rand()) / (RAND_MAX * 1.0);
                rho(1, labels[z][y][x]) = 1 - rho(0, labels[z][y][x]);
            }
        }
    }
    
    
    lbBase_t rho0Tot=0.0;
    lbBase_t rho1Tot=0.0;
    for (int noNo = 1; noNo < nNodes; ++noNo) {
      rho0Tot+=rho(0,noNo);
      rho1Tot+=rho(1,noNo);
      //std::cout << noNo << " " << rho(0,noNo) << std::endl;
    }
    std::cout << "before init" << std::endl;
    std::cout << "rho0Tot = " <<rho0Tot<< std::endl;
    std::cout << "rho1Tot = " <<rho1Tot<< std::endl;
    std::cout << "rhoTot = " <<rho0Tot+rho1Tot<< std::endl;
    
    // rho(0, 1) = 0.7;
    // rho(1, 1) = 0.3;
    // -- Phase total velocity
    
    //lbBase_t zeroVec[LT::nD] = {0.0, 0.0};
    lbBase_t zeroVec[LT::nD] = {0.0, 0.0, 0.0};

    setFieldToConst(zeroVec, 0, vel);  // LBmacroscopic
    // -- Global color gradient
    setFieldToConst(zeroVec, 0, colorGrad);

    // -- set solid boundary
    setFieldToConst(solidBoundary, 0.7, 0, rho);
    setFieldToConst(solidBoundary, 0.3, 1, rho);

    // INITIATE LB FIELDS
    // -- phase 0
    initiateLbField(0, 0, 0, bulk, rho, vel, f);  // LBinitiatefield
    // -- phase 1
    initiateLbField(1, 1, 0, bulk, rho, vel, f);  // LBinitiatefield

    std::cout << "after init" << std::endl;
    printAsciiToScreen(nX, nY, nZ, 0, rho, labels, 0);
    std::cout << "nQ = " << LT::nQ <<std::endl;
    //printAsciiToScreen(nX, nY, rho, labels, 0.0);

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
            // UPDATE MACROSCOPIC DENSITIES
            lbBase_t rho0Node, rho1Node;
            // Calculate rho for each phase
            calcRho<LT>(&f(0,0,nodeNo), rho0Node);  // LBmacroscopic
            rho(0, nodeNo) = rho0Node; // save to global field
            calcRho<LT>(&f(1,0,nodeNo), rho1Node);  // LBmacroscopic
            rho(1, nodeNo) = rho1Node; // save to global field

            // Calculate color gradient
/*            lbBase_t cgTerm = (rho0Node - rho1Node)/(rho0Node + rho1Node);
            LT::gradPush(cgTerm, grid.neighbor(nodeNo), colorGrad); */
        }  // End for all bulk nodes

/*        for (int bndNo = 0; bndNo < solidBoundary.getNumNodes(); ++bndNo) {
            const int nodeNo = solidBoundary.nodeNo(bndNo);
            // Calculate color gradient
            lbBase_t rho0Node, rho1Node;
            rho0Node = rho(0, nodeNo);
            rho1Node = rho(1, nodeNo);
            lbBase_t cgTerm = (rho0Node - rho1Node)/(rho0Node + rho1Node);
            LT::gradPush(cgTerm, grid.neighbor(nodeNo), colorGrad);
        } */


        for (int bulkNo = 0; bulkNo < bulk.nElements(); bulkNo++ ) {
            const int nodeNo = bulk.nodeNo(bulkNo); // Find current node number

            // UPDATE MACROSCOPIC VARIABLES
            lbBase_t rhoNode, rho0Node, rho1Node;
            //lbBase_t forceNode[LT::nD] = {0.0, 0.0};
            lbBase_t forceNode[LT::nD] = {0.0, 0.0, 0.0};
            lbBase_t velNode[LT::nD];
            lbBase_t fTot[LT::nQ];


            // Set the local total lb distribution
            for (int q = 0; q < LT::nQ; ++q)
                fTot[q] = f(0, q, nodeNo) + f(1, q, nodeNo);

            // Set densities
            rho0Node = rho(0, nodeNo);
            rho1Node = rho(1, nodeNo);
            // Total density
            rhoNode = rho0Node + rho1Node;
            // Total velocity and force
            forceNode[0] = 0.0;//(rho0Node - rho1Node) / rhoNode * bodyForce(0, 0, 0);
            forceNode[1] = rho0Node/rhoNode * bodyForce(0, 1, 0);
            forceNode[2] = rho0Node/rhoNode * bodyForce(0, 2, 0);

            calcVel<LT>(fTot, rhoNode, velNode, forceNode);  // LBmacroscopic
            vel(0, 0, nodeNo) = velNode[0];
            vel(0, 1, nodeNo) = velNode[1];
            vel(0, 2, nodeNo) = velNode[2];

            // CALCULATE BGK COLLISION TERM
            // Mean collision time /rho_tot/\nu_tot = \sum_s \rho_s/\nu_s
            lbBase_t tau;
            tau = LT::c2Inv * rhoNode / (rho0Node*nu0Inv + rho1Node*nu1Inv) + 0.5;

            lbBase_t omegaBGK[LT::nQ]; // Holds the collision values
            lbBase_t uu, cu[LT::nQ];  // pre-calculated values
            uu = LT::dot(velNode, velNode);  // Square of the velocity
            LT::cDotAll(velNode, cu);  // velocity dotted with lattice vectors
            calcOmegaBGK<LT>(fTot, tau, rhoNode, uu, cu, omegaBGK);  // LBcollision


            lbBase_t tmpSum = 0.0;
            for (int q=0; q < LT::nQ; ++q)
                tmpSum += omegaBGK[q];
            //std::cout << tmpSum << std::endl;

            // CALCULATE FORCE CORRECTION TERM
            lbBase_t deltaOmega[LT::nQ];
            lbBase_t  uF, cF[LT::nQ];
            uF = LT::dot(velNode, forceNode);
            LT::cDotAll(forceNode, cF);

            calcDeltaOmega<LT>(tau, cu, uF, cF, deltaOmega);  // LBcollision

            // -- calculate the normelaized color gradient
            lbBase_t CGNorm, CG2;
            // lbBase_t* colorGradNode = &colorGrad(0, 0, nodeNo);
            lbBase_t colorGradNode[LT::nD];

            lbBase_t scalarTmp[LT::nQ];
            for (int q = 0; q < LT::nQ; ++q) {
                int neigNode = grid.neighbor(q, nodeNo);
                lbBase_t rho0Neig = rho(0, neigNode);
                lbBase_t rho1Neig = rho(1, neigNode);
                scalarTmp[q] = (rho0Neig - rho1Neig) / (rho0Neig + rho1Neig);
            }
            LT::grad(scalarTmp, colorGradNode);
            CG2 = LT::dot(colorGradNode, colorGradNode);
            CGNorm = sqrt(CG2);

            // Normalization of colorGrad
            for (int d = 0; d < LT::nD; ++d)
                colorGradNode[d] /= CGNorm + (CGNorm < lbBaseEps);

            // CALCULATE SURFACE TENSION PERTURBATION
            lbBase_t deltaOmegaST[LT::nQ];
            lbBase_t bwCos[LT::nQ];
            lbBase_t cCGNorm[LT::nQ];
            LT::cDotAll(colorGradNode, cCGNorm);
            lbBase_t AF0_5 = 2.25 * CGNorm * sigma / tau;
            lbBase_t rhoFacBeta = beta * rho0Node * rho1Node / rhoNode;

            for (int q = 0; q < LT::nQNonZero_; ++q) {
                deltaOmegaST[q] = AF0_5 * (LT::w[q] * cCGNorm[q]*cCGNorm[q] - LT::B[q] );
                bwCos[q] = rhoFacBeta * LT::w[q] * cCGNorm[q] /  LT::cNorm[q];
            }
            deltaOmegaST[LT::nQNonZero_] = -AF0_5 * LT::B[LT::nQNonZero_];
            bwCos[LT::nQNonZero_] = 0.0; // This should be zero by default


            // COLLISION AND PROPAGATION
            lbBase_t c0, c1;
            c0 = (rho0Node/rhoNode);  // Concentration of phase 0
            c1 = (rho1Node/rhoNode);  // Concentration of phase 1
            for (int q = 0; q < LT::nQ; q++) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = c0 * (fTot[q] + omegaBGK[q] + deltaOmega[q] +  deltaOmegaST[q]) +  bwCos[q];
                fTmp(1, q,  grid.neighbor(q, nodeNo)) = c1 * (fTot[q] + omegaBGK[q] + deltaOmega[q] +  deltaOmegaST[q]) -  bwCos[q];
            }

	    
            // CLEAR GLOBAL FIELDS
            // Color gradient
/*            colorGrad(0,0,nodeNo)=0.0;
            colorGrad(0,1,nodeNo)=0.0;
            colorGrad(0,2,nodeNo)=0.0; */


        } // End nodes

        // PRINT
        if ((i % 50)  == 0)
//            printAsciiToScreen(nX, nY, nZ, i+1, rho, labels, 0);
        printAsciiToScreen(nX, nY, rho, labels[0], 0.1);

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // BOUNDARY CONDITIONS
        boundary.apply(0, f, grid);  // LBboundary
        boundary.apply(1, f, grid);
    } // End iterations
    // -----------------END MAIN LOOP------------------

    std::cout << nNodes << std::endl;
    
    
    //    std::cout << std::setprecision(5) << vel(0, 1, labels[1][0])  << " " << vel(0, 0, labels[50][0])  << std::endl;
    //std::cout << rho(0, labels[1][1][1]) << " " << rho(1, labels[0][10][1])  <<  " " << rho(0, labels[1][10][1]) <<  " " << rho(1, labels[1][0][1]) << std::endl;
    /*    int tmp;
    for (int y = 0; y < 30; ++y) {
        for (int x = 0; x < 30; ++x) {
            tmp = 10*rho(1, labels[y][x]);
            std::cout << std::setw(4) << std::setprecision(2) << tmp*0.1;
        }
        std::cout << std::endl << std::endl;
    } */

    // CLEANUP
    /*
    deleteNodeLabel(nX, nY, labels);
    deleteGeometry(nX, nY, geo);
    */
    deleteNodeLabel(nX, nY, nZ, labels);
    deleteGeometry(nX, nY, nZ, geo);
    
    return 0;
}

