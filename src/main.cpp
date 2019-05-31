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
// TO RUN PROGRAM: type "mpirun -np <#procs> badchimpp" in command
// line in main directory
//
// //////////////////////////////////////////////

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "LBbndmpi.h"
#include "LBboundary.h"
#include "LBcollision.h"
#include "LBcollision2phase.h"
#include "LBlatticetypes.h"
#include "LBfield.h"
#include "LBgeometry.h"
#include "LBgrid.h"
#include "LBhalfwaybb.h"
#include "LBinitiatefield.h"
#include "LBmacroscopic.h"
#include "LBnodes.h"
#include "LBsnippets.h"
#include "LButilities.h"


#include "Input.h"
#include "Output.h"
//#include "Mpi_class.h"
//#include "Geo.h"
//#include "Field.h"

#include<algorithm> // std::max

// SET THE LATTICE TYPE
// #define LT D2Q9
#define LT D3Q19


int main()
{
    std::cout << "Begin test Two phase new" << std::endl;

    // SETUP THE INPUT AND OUTPUT PATHS
    std::string chimpDir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/";


    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outDir = chimpDir + "output/rho_val_";

    // read input files
    Input input(inputDir + "input.dat");

    // initialize MPI
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // READ BOUNDARY FILES WITH MPI INFORMATION
    MpiFile<LT> rankFile(mpiDir + "rank.mpi");
    MpiFile<LT> localFile(mpiDir + "rank_" + std::to_string(myRank) + "_labels.mpi");
    MpiFile<LT> globalFile(mpiDir + "node_labels.mpi");

    // SETUP GRID
    auto grid  = Grid<LT>::makeObject(localFile, rankFile);
    // SETUP NODE
    auto nodes = Nodes<LT>::makeObject(localFile, rankFile, myRank, grid);

    // SETUP MPI BOUNDARY
    BndMpi<LT> mpiBoundary(myRank);

    mpiBoundary.setup(localFile, globalFile, rankFile, nodes,  grid);

    // SETUP BOUNCE BACK BOUNDARY (fluid boundary)
    HalfWayBounceBack<LT> bbBnd = makeFluidBoundary<HalfWayBounceBack>(myRank, grid);

    // SETUP SOLID BOUNDARY
    std::vector<int> solidBnd = findSolidBndNodes(myRank, nodes, grid);

    // SETUP BULK NODES
    std::vector<int> bulkNodes = findBulkNodes(myRank, nodes);

    // SETUP CONST PRESSURE AND FLUID SINK
    std::vector<int> constDensNodes, sourceNodes;
    for (auto bulkNode: bulkNodes) {
        int y = grid.pos(bulkNode, 1);
        if (nodes.getRank(bulkNode) == myRank) {
            if ( y < 3) // sink
                sourceNodes.push_back(bulkNode);
            if ( y > (rankFile.dim(1) - 2*1 - 4) ) // Const pressure
                constDensNodes.push_back(bulkNode);
        }
    }


    // READ INPUT FILE
    // Scalar source
    ScalarField Q(2, grid.size());
    for (int n = 0; n < Q.size(); ++n) {
        Q(0, n) = 0.0;
        Q(1, n) = 0.0;
    }
    // Vector source
    VectorField<LT> bodyForce(1, 1);
    std::vector<lbBase_t> tmpVec = input["fluid"]["bodyforce"];

    bodyForce(0, 0) = tmpVec;

    int nIterations = static_cast<int>( input["iterations"]["max"]);

    lbBase_t tau0 = input["fluid"]["tau"][0];
    lbBase_t tau1 = input["fluid"]["tau"][1];
    lbBase_t sigma = input["fluid"]["sigma"];
    lbBase_t beta = input["fluid"]["beta"];

    // SET DERIVED VARAIBLES
    lbBase_t nu0Inv = 1.0 / (LT::c2 * (tau0 - 0.5));
    lbBase_t nu1Inv = 1.0 / (LT::c2 * (tau1 - 0.5));

    // SETUP LB FIELDS
    LbField<LT> f(2, grid.size());  // LBfield
    LbField<LT> fTmp(2, grid.size());  // LBfield

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(2, grid.size()); // LBfield
    VectorField<LT> vel(1, grid.size()); // LBfield
    ScalarField cgField(1, grid.size()); // LBfield


    // FILL MACROSCOPIC FIELDS
    //   Fluid densities and velocity
    std::srand(8549388);
    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 1.0; // (1.0 * std::rand()) / (RAND_MAX * 1.0);
        rho(1, nodeNo) = 1 - rho(0, nodeNo);

        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }
    //   Solid boundary (Wettability)
    for (auto nodeNo: solidBnd) {
        rho(0, nodeNo) = 1.0;
        rho(1, nodeNo) = 1.0 - rho(0, nodeNo);
    }


    // INITIATE LB FIELDS
    // -- phase 0
    //f(0,bulk).initializeLB(rho(0,bulk), vel(0,bulk))
    initiateLbField(0, 0, 0, bulkNodes, rho, vel, f);  // LBinitiatefield
    // -- phase 1
    initiateLbField(1, 1, 0, bulkNodes, rho, vel, f);  // LBinitiatefield

    // JLV
    // set up output
    std::vector<std::vector<int>> node_pos; node_pos.reserve(bulkNodes.size());
    for (const auto& node:bulkNodes) {
      node_pos.push_back(grid.pos(node));
    }
//    for (auto i=0; i<bulkNodes.size(); ++i) {
//      node_pos[i] = grid.pos(bulkNodes[i]);
//    }
    Output output({globalFile.dim(0), globalFile.dim(1), globalFile.dim(2)}, "out", myRank, nProcs-1, node_pos);
    output.add_file("fluid");
    output["fluid"].add_variable("rho0", rho.get_iterator(0, bulkNodes), sizeof(lbBase_t), 1);
    output["fluid"].add_variable("rho1", rho.get_iterator(1, bulkNodes), sizeof(lbBase_t), 1);
    output.write("fluid", 0);
    // JLV

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

        for (auto nodeNo : bulkNodes) {
            // UPDATE MACROSCOPIC DENSITIES
            // Calculate rho for each phase
            lbBase_t rho0Node = rho(0, nodeNo) = calcRho<LT>(&f(0,0,nodeNo));  // LBmacroscopic
            lbBase_t rho1Node = rho(1, nodeNo) = calcRho<LT>(&f(1,0,nodeNo));  // LBmacroscopic

            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }  // End for all bulk nodes

        // Need to set the density when when calculating the densities
        // Predefine two regions:
        //   1) Pressure is constant
        //   2) Constant sink
        // Strategy - When mass is removed, then remove both oil and water
        //          - If mass is added at the top then then only add oil.

        for (auto nodeNo: constDensNodes) {
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);

            setConstDensity(Q(0, nodeNo), Q(1, nodeNo), rho0Node, rho1Node, 1.0);
            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }

        for (auto nodeNo: sourceNodes) {
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);
            setConstSource(Q(0, nodeNo), Q(1, nodeNo), rho0Node, rho1Node, -1e-4);
            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }


        for (auto nodeNo: solidBnd) {
            const lbBase_t rho0Node = rho(0, nodeNo);
            const lbBase_t rho1Node = rho(1, nodeNo);
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }


        //  MPI: COMMUNCATE SCALAR 'cgField'
        mpiBoundary.communciateScalarField(cgField);


        for (auto nodeNo: bulkNodes) {

            // Set the local total lb distribution
            std::vector<lbBase_t> fTot(LT::nQ);
            for (int q = 0; q < LT::nQ; ++q)
                fTot[q] = f(0, q, nodeNo) + f(1, q, nodeNo);

            // UPDATE MACROSCOPIC VARIABLES
            // -- densities
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);
            // -- total density
            lbBase_t rhoNode = rho0Node + rho1Node;
            // -- force
            std::vector<lbBase_t> forceNode = setForceGravity(rho0Node, rho1Node, bodyForce, 0);

            // -- velocity
            std::vector<lbBase_t> velNode = calcVel<LT>(fTot, rhoNode, forceNode);  // LBmacroscopic
            vel(0, nodeNo) = velNode;

            // Correct mass density for mass source
            lbBase_t q0Node = Q(0, nodeNo);
            lbBase_t q1Node = Q(1, nodeNo);
            rho(0, nodeNo) = rho0Node += 0.5*q0Node;
            rho(1, nodeNo) = rho1Node += 0.5*q1Node;


            // CALCULATE BGK COLLISION TERM
            // Mean collision time /rho_tot/\nu_tot = \sum_s \rho_s/\nu_s
            lbBase_t tau = LT::c2Inv * rhoNode / (rho0Node*nu0Inv + rho1Node*nu1Inv) + 0.5;

            lbBase_t uu = LT::dot(velNode, velNode);  // Square of the velocity
            std::vector<lbBase_t> cu = LT::cDotAll(velNode);  // velocity dotted with lattice vectors
            std::vector<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fTot, tau, rhoNode, uu, cu);  // LBcollision


            // CALCULATE FORCE CORRECTION TERM
            lbBase_t  uF = LT::dot(velNode, forceNode);
            std::vector<lbBase_t>  cF = LT::cDotAll(forceNode);
            std::vector<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);  // LBcollision

            // CALCULATE MASS SOURCE CORRECTION TERM
            std::vector<lbBase_t> deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tau, cu, uu, q0Node);
            std::vector<lbBase_t> deltaOmegaQ1 = calcDeltaOmegaQ<LT>(tau, cu, uu, q1Node);

            // CALCULATE SURFACE TENSION PERTURBATION
            // -- calculate the normalized color gradient
            std::vector<lbBase_t> colorGradNode = grad(cgField, 0, nodeNo, grid);
            lbBase_t CGNorm = vecNorm<LT>(colorGradNode);
            lbBase_t CGNormInv = 1.0/(CGNorm + (CGNorm < lbBaseEps));
            for (auto &cg: colorGradNode)
                cg *= CGNormInv;

            std::vector<lbBase_t> cCGNorm = LT::cDotAll(colorGradNode);

            std::vector<lbBase_t> deltaOmegaST = calcDeltaOmegaST<LT>(tau, sigma, CGNorm, cCGNorm);
            std::vector<lbBase_t> deltaOmegaRC = calcDeltaOmegaRC<LT>(beta, rho0Node, rho1Node, rhoNode, cCGNorm);

            // COLLISION AND PROPAGATION
            lbBase_t c0, c1;
            c0 = (rho0Node/rhoNode);  // Concentration of phase 0
            c1 = (rho1Node/rhoNode);  // Concentration of phase 1
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = c0 * (fTot[q] + omegaBGK[q] + deltaOmegaF[q] +  deltaOmegaST[q]) +  deltaOmegaRC[q]  + deltaOmegaQ0[q];
                fTmp(1, q,  grid.neighbor(q, nodeNo)) = c1 * (fTot[q] + omegaBGK[q] + deltaOmegaF[q] +  deltaOmegaST[q]) -  deltaOmegaRC[q]  + deltaOmegaQ1[q];
            }
        } // End nodes

        // PRINT

        if ( (i % static_cast<int>(input["iterations"]["write"])) == 0) {
          if (myRank==0)
            std::cout << "PLOT AT ITERATION : " << i << std::endl;
          std::string tmpName(outDir);
          tmpName += std::to_string(myRank) + "_" + std::to_string(i);
          tmpName += ".dat";
          std::ofstream ofs;
          ofs.open(tmpName);
          if (!ofs) {
            std::cout << "Error: could not open file: " << tmpName << std::endl;
            return 1;
          }
          for (auto nodeNo: bulkNodes) {
            ofs << std::setprecision(23) << grid.pos(nodeNo, 0) << " " << grid.pos(nodeNo, 1) << " " << grid.pos(nodeNo, 2) << " " << rho(0, nodeNo) << " " << rho(1, nodeNo)
                            << " " << vel(0, 0, nodeNo) << " " << vel(0, 1, nodeNo) << " " << vel(0, 2, nodeNo) << std::endl;
          }
          ofs.close();

          // JLV
          output.write("fluid", i);
          // JLV
        }

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // MPI Boundary
        // Hente verdier hente fra ghost
        // Sette i bulk
        mpiBoundary.communicateLbField(grid, f, 0);
        mpiBoundary.communicateLbField(grid, f, 1);

        // BOUNDARY CONDITIONS
        bbBnd.apply(0, f, grid);  // LBboundary
        bbBnd.apply(1, f, grid);


    } // End iterations
    // -----------------END MAIN LOOP------------------

    //mpi.end();

    MPI_Finalize();

    return 0;
}

