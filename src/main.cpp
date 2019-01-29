//
// TO RUN PROGRAM: type bin/runner in command line in main directory
//
// TO DO :
//        1) [OK] Define lattice information
//        2) Boundary nodes
//        3) What about ghost nodes. Do we need them?
//        3b) [OK] uansett ghost nodes info bør ligge i grid
//        4) [OK] Change function variables to longer descreptive names, (so we can easly see the difference input-variables and loop-iterators)
//        5) Se på denne: https://upload.wikimedia.org/wikipedia/commons/7/7d/OptimizingCpp.pdf
//        6) Legger bouncback etter hverandre
//         |b bHat
//
// GRID :
// /////////////////////////////////////////

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
    std::cout << "Begin test";
    std::cout << std::endl;


    // INPUT DATA
    int nIterations;
    int nX, nY;
    lbBase_t tau;

    lbBase_t force[2] = {1.0e-8, 0.0};
    VectorField<LT> bodyForce(1,1);
    bodyForce(0, 0, 0) = force[0];
    bodyForce(0, 1, 0) = force[1];

    nIterations = 10000;
    nX = 250; nY = 101;

    tau = 0.7;

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
    std::cout << "NUMBER OF NODES = " << nNodes << std::endl;

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
    LbField<LT> f(1, nNodes); // Bør f-field vite at det er nQ. Dette kan nok gjøre at ting blir gjort raskere
    LbField<LT> fTmp(1, nNodes); // Do not need to be a LbField (as f), since its just a memory holder, that is immediately swaped after propagation.

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(1, nNodes);
    VectorField<LT> vel(1, nNodes);

    // FILL MACROSCOPIC FIELDS
    setFieldToConst(1.0, 0, rho);
    lbBase_t velTmp[LT::nD] = {0.0, 0.0};
    setFieldToConst(velTmp, 0, vel);

    // INITIATE LB FIELDS
    initiateLbField(0, bulk, rho, vel, f);

    // -----------------MAIN LOOP------------------
    /* Comments to main loop:
     * Calculation of cu is kept outside of calcOmega and calcDeltaOmega
     *  to avoid recalculation  of the dot-product.
     * It is therefore natural to give all scalar products as inputs to the
     *  calcOmega and calcDeltaOmega.
     */

    for (int i = 0; i < nIterations; i++) {
      // oneIteration(tau, bodyForce, bulk, grid, boundary, rho, vel, f, fTmp);
        for (int bulkNo = 0; bulkNo < bulk.nElements(); bulkNo++ ) {
            const int nodeNo = bulk.nodeNo(bulkNo); // Find current node number

            // UPDATE MACROSCOPIC VARAIBLES
            lbBase_t rhoNode;
            lbBase_t velNode[LT::nD];
            lbBase_t *forceNode = &bodyForce(0,0,0);

            calcRho<LT>(&f(0,0,nodeNo), rhoNode);
            calcVel<LT>(&f(0,0,nodeNo), rhoNode, velNode, forceNode);

            // Sets the global macroscopic values
            storeGlobalValues(0, nodeNo, rhoNode, velNode, rho, vel);

            // CALCULATE BGK COLLISION TERM
            lbBase_t omegaBGK[LT::nQ]; // Holds the collision values
            lbBase_t uu, cu[LT::nQ];  // pre-calculated values
            uu = LT::dot(velNode, velNode);  // Square of the velocity
            LT::cDotAll(velNode, cu);  // velocity dotted with lattice vectors
            calcOmegaBGK<LT>(&f(0, 0, nodeNo), tau, rhoNode, uu, cu, omegaBGK);

            // CALCULATE FORCE CORRECTION TERM
            lbBase_t deltaOmega[LT::nQ];
            lbBase_t  uF, cF[LT::nQ];
            uF = LT::dot(velNode, forceNode);
            LT::cDotAll(forceNode, cF);

            calcDeltaOmega<LT>(tau, cu, uF, cF, deltaOmega);

            // COLLISION AND PROPAGATION
            for (int q = 0; q < LT::nQ; q++) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = f(0, q, nodeNo) + omegaBGK[q] + deltaOmega[q];
            }


        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);

        // BOUNDARY CONDITIONS
        boundary.apply(0, f, grid); 

    } // End iterations
    // -----------------END MAIN LOOP------------------
    
    std::cout << std::setprecision(3) << vel(0, 0, labels[10][8])  << std::endl;


   // CLEANUP
    deleteNodeLabel(nX, nY, labels);
    deleteGeometry(nX, nY, geo);

    return 0;
}

