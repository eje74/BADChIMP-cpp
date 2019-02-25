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

#include "Input.h"
#include "Output.h"
#include "Mpi.h"
#include "Geo.h"

#include "LBmodel.h"


// SET THE LATTICE TYPE
//#define LT D2Q9
#define LT D3Q19


int main(int argc, char *argv[])
{
    std::cout << "Begin test Two phase" << std::endl;

    // read input files
    Input input("/home/ejette/Programs/GITHUB/badchimpp/input.dat");
    // initialize MPI
    Mpi mpi(&argc, &argv, input["mpi"]["procs"]);

    // read geo and create node-array
    Geo geo2("/home/ejette/Programs/GITHUB/badchimpp/geo30x30x30wWall.dat", mpi);

    geo2.print_limits();

    lbBase_t force[3] = {0.0, 0.0, 0.0}; //{1.0e-5, 1.e-4, 0.0};
    
    int nIterations = static_cast<int>( input["iterations"]["max"]);
    int nX = geo2.get_N(0);
    int nY = geo2.get_N(1);
    int nZ = geo2.get_N(2);

    lbBase_t tau0 = input["fluid"]["tau"][0];
    lbBase_t tau1 = input["fluid"]["tau"][1];
    lbBase_t sigma = input["fluid"]["sigma"];
    lbBase_t beta = input["fluid"]["beta"];

    // SETUP GEOMETRY
    int *** geo;
    newGeometry(nX, nY, nZ, geo); // LBgeometry
    geo2.set_node_values_v2<LT>();
    geo2.export_geo_to_3D(geo);
    
    int *** labels;
    newNodeLabel(nX, nY, nZ, labels); // LBgeometry
    geo2.set_labels();
    geo2.export_labels_to_3D(labels);
    int nBulkNodes = geo2.get_num_bulk_nodes();

    int nNodes = geo2.get_num_nodes();
    // USE NODENO 0 AS DUMMY NODE.
    // Add one extra storage for the default dummy node
    nNodes += 1;


    // SETUP GRID
    std::cout << "NUMBER OF NODES = " << nNodes << std::endl;
    Grid<LT> grid(nNodes); // object declaration: neigList_ stores node numbers of neighbors;  pos_ stores Cartesian coordinates of given node

    setupGrid(nX, nY, nZ, labels, grid); // LBgeometry, maybe move to LBgrid?
	
    
    // SETUP BULK
    std::cout << "NUMBER OF Bulk NODES = " << nBulkNodes << std::endl;
    Bulk bulk(nBulkNodes); // object declaration: nBulkNodes_ stores number of Bulk nodes; bulkNode_ stores node numbers of all bulk nodes

    setupBulk(nX, nY, nZ, geo, labels, bulk); // LBgeometry
    
    // SETUP BOUNDARY
    int nBoundary = 0;
    
    nBoundary = nBoundaryNodes(1, nX, nY, nZ, geo);

    std::cout << "NUMBER OF BOUNDARY NODES = " << nBoundary << std::endl;
    HalfWayBounceBack<LT> boundary( nBoundary ); // LBhalfwaybb
    setupBoundary(1, nX, nY, nZ, geo, labels, grid, boundary); // LBgeometry
    
    // SETUP SOLID BOUNDARY
    int nSolidBoundary = 0;
    nSolidBoundary = nBoundaryNodes(3, nX, nY, nZ, geo);

    
    Boundary<LT> solidBoundary( nSolidBoundary );
    setupBoundary(3, nX, nY, nZ, geo, labels, grid, solidBoundary); // LBgeometry


    // BEGIN MAKE MODEL CLASS
    TwoPhaseCG<LT> fluidModel(nNodes);
    // -- Fill macroscopic fields
    fluidModel.setTau(tau0, tau1);
    fluidModel.setNuInv();
    fluidModel.setSigma(sigma);
    fluidModel.setBeta(beta);

    fluidModel.setBodyForce(force);
    fluidModel.initRho(grid);  // Instead of grid should we send a tag list
    fluidModel.initVelocity(grid);
    fluidModel.initLbField(grid);
    fluidModel.setRhoToConst(solidBoundary, 0.7, 0);
    fluidModel.setRhoToConst(solidBoundary, 0.3, 1);
    // END MAKE MODEL CLASS

    // initialize output
    geo2.labels_ = labels;
    Output output("out", mpi, geo2);
    output.add_file("fluid");
    output["fluid"].add_variables({"rho0","rho1"}, {&fluidModel.rho(0,0),&fluidModel.rho(1,0)}, {sizeof(fluidModel.rho(0,0)),sizeof(fluidModel.rho(1,0))}, {1,1}, {fluidModel.rho.nFields_,fluidModel.rho.nFields_});
    output.write("fluid",0);



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
      
        for (int bulkNo = 0; bulkNo < bulk.nElements(); bulkNo++ ) {  // Change nElements to nNodes?
            const int nodeNo = bulk.nodeNo(bulkNo); // Find current node number
            // UPDATE MACROSCOPIC DENSITIES
            fluidModel.calcRho(nodeNo);
            // Calculate color gradient kernel
            fluidModel.calcCgKernel(nodeNo);
        }  // End for all bulk nodes

        for (int bndNo = 0; bndNo < solidBoundary.getNumNodes(); ++bndNo) { // Change getNumNodes to nNodes ?
            const int nodeNo = solidBoundary.nodeNo(bndNo);
            // Calculate color gradient kernel
            fluidModel.calcCgKernel(nodeNo);
        }  // End for solid rim


        for (int bulkNo = 0; bulkNo < bulk.nElements(); bulkNo++ ) { // Change nElements to nNodes?
            const int nodeNo = bulk.nodeNo(bulkNo); // Find current node number

            // MODEL
            fluidModel.setLocalFtot(nodeNo);
            fluidModel.fetchRho(nodeNo);
            fluidModel.calcForce();
            fluidModel.calcVel(nodeNo); // Can be split into calc local and set global velocity. Test

            fluidModel.calcTauEff();
            fluidModel.calcOmegaBGK();
            fluidModel.calcDeltaOmega();
            fluidModel.calcDeltaOmegaST(nodeNo, grid);
            fluidModel.CollisionPropagation(nodeNo, grid);

        } // End nodes


/*        if ( (i % static_cast<int>(input["iterations"]["write"])) == 0) {
          output.write_all(i);
          std::cout << output.get_filename("fluid") << std::endl;
        } */

        // Swap data_ from fTmp to f;
        fluidModel.getF().swapData(fluidModel.getFtmp());  // LBfield

        // BOUNDARY CONDITIONS
        boundary.apply(0, fluidModel.getF(), grid);  // LBboundary
        boundary.apply(1, fluidModel.getF(), grid);

    } // End iterations
    // -----------------END MAIN LOOP------------------

    std::cout << fluidModel.rho(0, labels[1][1][1]) << " " << fluidModel.rho(1, labels[1][10][1])  <<  " " << fluidModel.rho(0, labels[1][10][1]) <<  " " << fluidModel.rho(1, labels[1][0][1]) << std::endl;

    // CLEANUP
    deleteNodeLabel(nX, nY, nZ, labels);
    deleteGeometry(nX, nY, nZ, geo);
    
    return 0;
}

