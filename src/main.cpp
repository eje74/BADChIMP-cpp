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

#define dxqy D2Q9


int main()
{
    std::cout << "Begin test";
    std::cout << std::endl;


    // INPUT DATA
    int nIterations;
    int nX, nY;
    lbBase_t tau;

    lbBase_t force[2] = {1.0e-8, 0.0};
    VectorField<D2Q9> F(1,1);
    F(0, 0, 0) = force[0];
    F(0, 1, 0) = force[1];

    nIterations = 10000;
    nX = 250; nY = 101;

    tau = 0.7;

    // SETUP GEOMETRY
    int ** geo;
    newGeometry(nX, nY, geo); // LBgeometry
    inputGeometry(nX, nY, geo); // LBgeometry
    analyseGeometry<D2Q9>(nX, nY, geo); // LBgeometry

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
    Grid<D2Q9> grid(nNodes); // object declaration: neigList_ stores node numbers of neighbors;  pos_ stores Cartesian coordinates of given node
    setupGrid(nX, nY, labels, grid); // LBgeometry, maybe move to LBgrid?

    // SETUP BULK
    std::cout << "NUMBER OF Bulk NODES = " << nBulkNodes << std::endl;
    Bulk bulk(nBulkNodes); // object declaration: nBulkNodes_ stores number of Bulk nodes; bulkNode_ stores node numbers of all bulk nodes
    setupBulk(nX, nY, geo, labels, bulk); // LBgeometry

    // SETUP BOUNDARY
    int nBoundary = nBoundaryNodes(1, nX, nY, geo);
    std::cout << "NUMBER OF BOUNDARY NODES = " << nBoundary << std::endl;
    HalfWayBounceBack<D2Q9> boundary( nBoundary ); // LBhalfwaybb
    setupBoundary(1, nX, nY, geo, labels, grid, boundary); // LBgeometry

    // SETUP FIELDS
    LbField<D2Q9> f(1, nNodes); // Bør f-field vite at det er nQ. Dette kan nok gjøre at ting blir gjort raskere
    LbField<D2Q9> fTmp(1, nNodes); // Do not need to be a LbField (as f), since its just a memory holder, that is immediately swaped after propagation.
    VectorField<D2Q9> vel(1, nNodes); // Bør vel-field vite at det er 2D. Dette kan nok gjøre at ting blir gjort raskere
    ScalarField rho(1, nNodes); // Dette kan fikses for 'single rho fields' slik at man slipper å skrive rho(0, pos)

    // INIT FILEDS
    lbBase_t velTmp[D2Q9::nD] = {0.0, 0.0};
    for (int n = 0; n < bulk.nElements(); n++ ) {
        const int nodeNo = bulk.nodeNo(n);
        lbBase_t cu[D2Q9::nQ], uu;
        D2Q9::cDotAll(velTmp, cu);
        uu = D2Q9::dot(velTmp, velTmp);
        for (int q = 0; q < D2Q9::nQ; q++) {
            f(0, q, nodeNo) = D2Q9::w[q] * (1.0 + D2Q9::c2Inv*cu[q] + D2Q9::c4Inv0_5*(cu[q]*cu[q] - D2Q9::c2*uu)); // f_eq
        }
    }

    // -----------------MAIN LOOP------------------
    for (int i = 0; i < nIterations; i++) {
        oneIteration(tau, F, bulk, grid, boundary, rho, vel, f, fTmp);

        /* for (int bulkNo = 0; bulkNo < bulk.nElements(); bulkNo++ ) {
            // Find current node number
            const int nodeNo = bulk.nodeNo(bulkNo);

            // UPDATE MACROSCOPIC VARAIBLES
            lbBase_t rhoNode;
            lbBase_t velNode[D2Q9::nD];

            calcRho<D2Q9>(&f(0,0,nodeNo), rhoNode);
            calcVel<D2Q9>(&f(0,0,nodeNo), rhoNode, velNode, force);

            for (int d = 0; d < 2; ++d)
                vel(0, d, nodeNo) = velNode[d];


            // Calculate the BGK collision part
            lbBase_t omegaBGK[D2Q9::nQ]; // Holds the collision values
            lbBase_t uu, cu[D2Q9::nQ];  // pre-calculated values
            uu = D2Q9::dot(velNode, velNode);  // Square of the velocity
            D2Q9::cDotAll(velNode, cu);  // velocity dotted with lattice vectors

            calcOmegaBGK<D2Q9>(&f(0, 0, nodeNo), tau, rhoNode, uu, cu, omegaBGK);


            lbBase_t deltaOmega[D2Q9::nQ];
            lbBase_t  uF, cF[D2Q9::nQ];
            uF = D2Q9::dot(velNode, force);
            D2Q9::cDotAll(force, cF);

            calcDeltaOmega<D2Q9>(tau, cu, uF, cF, deltaOmega);

            // * Collision and propagation:
            for (int q = 0; q < D2Q9::nQ; q++) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = f(0, q, nodeNo) + omegaBGK[q] + deltaOmega[q];
            }


        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);

        // BOUNDARY CONDITIONS
        boundary.apply(0, f, grid); */


    } // End iterations (LOOP TYPE 1)
    // -----------------END MAIN LOOP------------------
    
    std::cout << std::setprecision(3) << vel(0, 0, labels[10][8])  << std::endl;

//    for (int y = 1; y < nY; ++y) {
//        std::cout << vel(0, 0, labels[y][8]) << std::endl;
//    }


   // CLEANUP
    deleteNodeLabel(nX, nY, labels);
    deleteGeometry(nX, nY, geo);

    return 0;
}

