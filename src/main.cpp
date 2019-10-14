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

#include<algorithm> // std::max

// SET THE LATTICE TYPE
// #define LT D2Q9
#define LT D3Q19


int main()
{
    // SETUP MPI
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // TEST PRINTOUT
    std::cout << "Begin test Two phase new" << std::endl;

    // SETUP THE INPUT AND OUTPUT PATHS
    // std::string chimpDir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/";
    //std::string chimpDir = "/home/ejette/Programs/GITHUB/badchimpp/";
    std::string chimpDir = "./";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outDir = chimpDir + "output/rho_val_";



    std::cout << "before read files" << std::endl;
    // READ BOUNDARY FILES WITH MPI INFORMATION
    MpiFile<LT> rankFile(mpiDir + "rank_" + std::to_string(myRank) + "_rank.mpi");
    MpiFile<LT> localFile(mpiDir + "rank_" + std::to_string(myRank) + "_local_labels.mpi");
    MpiFile<LT> globalFile(mpiDir + "rank_" + std::to_string(myRank) + "_global_labels.mpi");

    std::cout << "Sizes = " << rankFile.size() << " " << localFile.size() << " " << globalFile.size() << std::endl;

    // SETUP GRID
    std::cout << "grid" << std::endl;
    auto grid  = Grid<LT>::makeObject(localFile, rankFile, myRank);
    std::cout << "grid.size = " << grid.size() << std::endl;

    // SETUP NODE
    std::cout << "nodes" << std::endl;
    auto nodes = Nodes<LT>::makeObject(localFile, rankFile, myRank, grid);

    // SETUP MPI BOUNDARY
    std::cout << "mpi boundary" << std::endl;
    BndMpi<LT> mpiBoundary(myRank);
    mpiBoundary.setup(localFile, globalFile, rankFile, nodes,  grid);

    // SETUP BOUNCE BACK BOUNDARY (fluid boundary)
    std::cout << "bbBnd" << std::endl;
    HalfWayBounceBack<LT> bbBnd(findBulkNodes(nodes), nodes, grid); // = makeFluidBoundary<HalfWayBounceBack>(nodes, grid);

    // SETUP SOLID BOUNDARY
    std::vector<int> solidBnd = findSolidBndNodes(nodes);

    // SETUP BULK NODES
    std::vector<int> bulkNodes = findBulkNodes(nodes);

    // SETUP CONST PRESSURE AND FLUID SINK
    std::vector<int> constDensNodes, sourceNodes;
    
    for (auto bulkNode: bulkNodes) {
        int y = grid.pos(bulkNode, 1) - 1;
        if (nodes.getRank(bulkNode) == myRank) {
	  //if ( y < 3) // sink
	  if (y < 0){ // Never true
	    sourceNodes.push_back(bulkNode);
	  }
	  //if ( y > (rankFile.dim_global(1) - 2*1 - 4) ) // Const pressure
	  if ( y > (rankFile.dim_global(1)+10) ){ // Never true 
	    constDensNodes.push_back(bulkNode);
	  }
        }
    }
    
    
    // READ INPUT FILE
    Input input(inputDir + "input.dat");

    // Vector source
    VectorField<LT> bodyForce(1, 1);
    bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);

    int nIterations = static_cast<int>( input["iterations"]["max"]);

    lbBase_t tau0 = input["fluid"]["tau"][0];
    lbBase_t tau1 = input["fluid"]["tau"][1];
    lbBase_t sigma = input["fluid"]["sigma"];
    lbBase_t beta = input["fluid"]["beta"];

    // SET DERIVED VARAIBLES
    lbBase_t nu0Inv = 1.0 / (LT::c2 * (tau0 - 0.5));
    lbBase_t nu1Inv = 1.0 / (LT::c2 * (tau1 - 0.5));

    // Scalar source
    ScalarField Q(2, grid.size());
    for (int n = 0; n < Q.size(); ++n) {
        Q(0, n) = 0.0;
        Q(1, n) = 0.0;
    }
    // SETUP LB FIELDS
    LbField<LT> f(1, grid.size());  // LBfield
    LbField<LT> fTmp(1, grid.size());  // LBfield
    LbField<LT> fInd(1, grid.size()); // LB indicator field
    LbField<LT> fIndTmp(1, grid.size()); // LB indicator field

    LbField<LT> g(1, grid.size()); // LBfield diffusion
    LbField<LT> gTmp(1, grid.size()); // LBfield diffusion

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(1, grid.size()); // LBfield
    VectorField<LT> vel(1, grid.size()); // LBfield
    ScalarField cgField(1, grid.size()); // LBfield
    ScalarField indField(1, grid.size()); // LB indicator field

    // FILL MACROSCOPIC FIELDS
    //   Fluid densities and velocity
    std::srand(8549388);
    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 1.0; //(1.0 * std::rand()) / (RAND_MAX * 1.0);
	int y = grid.pos(nodeNo, 1) - 1;
	if ( y < 100)
	  indField(0, nodeNo) = 1.0; // LB indicator field
	else
	  indField(0, nodeNo) = 0.0;
        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }
    //   Solid boundary (Wettability)
    for (auto nodeNo: solidBnd) {
        rho(0, nodeNo) = 1.0;
	indField(0, nodeNo) = 1.0;
    }

    // INITIATE LB FIELDS
    // -- phase 0
    initiateLbField(0, 0, 0, bulkNodes, rho, vel, f);  // LBinitiatefield
    // -- phase 1
    //initiateLbField(1, 1, 0, bulkNodes, rho, vel, f);  // LBinitiatefield
    // -- indicator field
    initiateLbField(0, 0, 0, bulkNodes, indField, vel, fInd);  // LBinitiatefield

    // JLV
    // set up output
    std::vector<std::vector<int>> node_pos; node_pos.reserve(bulkNodes.size());
    for (const auto& node:bulkNodes) {
      node_pos.push_back(grid.pos(node));
    }

    Output output(globalFile.dim_global(), "out", myRank, nProcs-1, node_pos);
    output.add_file("fluid");
    output["fluid"].add_variable("rho", rho.get_data(), rho.get_field_index(0, bulkNodes), 1);
    output["fluid"].add_variable("vel", vel.get_data(), vel.get_field_index(0, bulkNodes), LT::nD);
    output["fluid"].add_variable("indicator", indField.get_data(), indField.get_field_index(0, bulkNodes), 1);
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
            lbBase_t rhoNode = rho(0, nodeNo) = calcRho<LT>(f(0,nodeNo));  // LBmacroscopic
	    lbBase_t indicator0Node = indField(0, nodeNo) = calcRho<LT>(fInd(0,nodeNo));  // LBmacroscopic

            // Calculate color gradient kernel
            //cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
	    cgField(0, nodeNo) = 2*indicator0Node-1;
        }  // End for all bulk nodes

        // Need to set the density when when calculating the densities
        // Predefine two regions:
        //   1) Pressure is constant
        //   2) Constant sink
        // Strategy - When mass is removed, then remove both oil and water
        //          - If mass is added at the top then then only add oil.

	/*Phase sources Q at the Constant-Density (Pressure) nodes are calculated*/
	/*
	for (auto nodeNo: constDensNodes) {
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);

            setConstDensity(Q(0, nodeNo), Q(1, nodeNo), rho0Node, rho1Node, 1.0);
            // Calculate color gradient kernel
            cgField(0, nodeNo) = 0;//(rho0Node - rho1Node)/(rho0Node + rho1Node);
        }
	*/
	/*Phase sources Q at the Constant-rate nodes are calculated*/
	/*
	for (auto nodeNo: sourceNodes) {
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);
            //setConstSource(Q(0, nodeNo), Q(1, nodeNo), rho0Node, rho1Node, -1e-4);
	    setConstSource(Q(0, nodeNo), Q(1, nodeNo), rho0Node, rho1Node, -1e-6);
            // Calculate color gradient kernel
            cgField(0, nodeNo) = 0;//(rho0Node - rho1Node)/(rho0Node + rho1Node);
        }
	*/

        for (auto nodeNo: solidBnd) {
	  //const lbBase_t rhoNode = rho(0, nodeNo);
	  //const lbBase_t rho1Node = rho(1, nodeNo);
	  const lbBase_t indicator0Node = indField(0, nodeNo);
	  //cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
	  cgField(0, nodeNo) = 2*indicator0Node-1;
        }

        //  MPI: COMMUNCATE SCALAR 'cgField'
        mpiBoundary.communciateScalarField(cgField);


        for (auto nodeNo: bulkNodes) {

            // Set the local total lb distribution
            // Use std:valarray insatead of auto as valarray uses expression templates that do not work well with auto
            //std::valarray<lbBase_t> fTot = f(0, nodeNo) + f(1, nodeNo);
	    std::valarray<lbBase_t> fTot = f(0, nodeNo);

            // UPDATE MACROSCOPIC VARIABLES
            // -- densities
            lbBase_t rho0Node = rho(0, nodeNo);
            //lbBase_t rho1Node = rho(1, nodeNo);
	    lbBase_t indicator0Node = indField(0, nodeNo);
            // -- total density
            //lbBase_t rhoNode = rho0Node + rho1Node;
	    lbBase_t rhoNode = rho0Node;
            // -- force
            std::valarray<lbBase_t> forceNode = setForceGravity(rho0Node, rhoNode, bodyForce, 0);
            // -- velocity
            std::valarray<lbBase_t> velNode = calcVel<LT>(fTot, rhoNode, forceNode);  // LBmacroscopic

            vel.set(0, nodeNo) = velNode;

            // Correct mass density for mass source
            lbBase_t q0Node = Q(0, nodeNo);
            lbBase_t q1Node = Q(1, nodeNo);
            rho(0, nodeNo) = rho0Node += 0.5*q0Node;
            //rho(1, nodeNo) = rho1Node += 0.5*q1Node;


            // CALCULATE BGK COLLISION TERM
            // Mean collision time /rho_tot/\nu_tot = \sum_s \rho_s/\nu_s
            //lbBase_t tau = LT::c2Inv * rhoNode / (rho0Node*nu0Inv + rho1Node*nu1Inv) + 0.5;
	    lbBase_t tau = LT::c2Inv * 1 / (indicator0Node*nu0Inv + (1-indicator0Node)*nu1Inv) + 0.5;

            lbBase_t uu = LT::dot(velNode, velNode);  // Square of the velocity
            std::valarray<lbBase_t> cu = LT::cDotAll(velNode);  // velocity dotted with lattice vectors
            std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fTot, tau, rhoNode, uu, cu);  // LBcollision
	    std::valarray<lbBase_t> omegaBGKInd = calcOmegaBGK<LT>(fInd(0, nodeNo), tau, indicator0Node, uu, cu);  // LBcollisionIndicator

            // CALCULATE FORCE CORRECTION TERM
            lbBase_t  uF = LT::dot(velNode, forceNode);
            std::valarray<lbBase_t>  cF = LT::cDotAll(forceNode);
            std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);  // LBcollision

            // CALCULATE MASS SOURCE CORRECTION TERM
            std::valarray<lbBase_t> deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tau, cu, uu, q0Node);
            std::valarray<lbBase_t> deltaOmegaQ1 = calcDeltaOmegaQ<LT>(tau, cu, uu, q1Node);

            // CALCULATE SURFACE TENSION PERTURBATION
            // -- calculate the normalized color gradient
            std::valarray<lbBase_t> colorGradNode = grad(cgField, 0, nodeNo, grid);
            lbBase_t CGNorm = vecNorm<LT>(colorGradNode);
            colorGradNode *= 1.0/(CGNorm + (CGNorm < lbBaseEps));

            std::valarray<lbBase_t> cCGNorm = LT::cDotAll(colorGradNode);

            std::valarray<lbBase_t> deltaOmegaST = calcDeltaOmegaST<LT>(tau, sigma, CGNorm, cCGNorm);
            //std::valarray<lbBase_t> deltaOmegaRC = calcDeltaOmegaRC<LT>(beta, rhoNode, rhoNode, rhoNode, cCGNorm);
	    std::valarray<lbBase_t> deltaOmegaRCInd = calcDeltaOmegaRCInd<LT>(beta, indicator0Node, rhoNode, cCGNorm);

            // COLLISION AND PROPAGATION
            lbBase_t c0, c1;
            //c0 = (rho0Node/rhoNode);  // Concentration of phase 0
            //c1 = (rho1Node/rhoNode);  // Concentration of phase 1
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be
	      //fTmp(0, q,  grid.neighbor(q, nodeNo)) = c0 * (fTot[q] + omegaBGK[q] + deltaOmegaF[q] +  deltaOmegaST[q]) +  deltaOmegaRC[q]  + deltaOmegaQ0[q];
		fTmp(0, q,  grid.neighbor(q, nodeNo)) = fTot[q] + omegaBGK[q] +  deltaOmegaST[q] + deltaOmegaF[q];
                //fTmp(1, q,  grid.neighbor(q, nodeNo)) = c1 * (fTot[q] + omegaBGK[q] + deltaOmegaF[q] +  deltaOmegaST[q]) -  deltaOmegaRC[q]  + deltaOmegaQ1[q];
		fIndTmp(0, q,  grid.neighbor(q, nodeNo)) = fInd(0, q, nodeNo) + omegaBGKInd[q] +  deltaOmegaRCInd[q];
            }
        } // End nodes

        // PRINT

        if ( (i % static_cast<int>(input["iterations"]["write"])) == 0) {
	  
	  if (myRank==0)
            std::cout << "PLOT AT ITERATION : " << i << std::endl;
	  /*
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
                            << " " << vel(0, 0, nodeNo) << " " << vel(0, 1, nodeNo) << " " << vel(0, 2, nodeNo) << " " << nodes.getType(nodeNo) << std::endl;
          }
          ofs.close();
	  */
          // JLV
          output.write("fluid", i);
          // JLV
        }

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield
	fInd.swapData(fIndTmp);  // LBfield

        // MPI Boundary
        // Hente verdier hente fra ghost
        // Sette i bulk
        mpiBoundary.communicateLbField(grid, f, 0);
        //mpiBoundary.communicateLbField(grid, f, 1);
	mpiBoundary.communicateLbField(grid, fInd, 0);

        // BOUNDARY CONDITIONS
        bbBnd.apply(0, f, grid);  // LBboundary
        //bbBnd.apply(1, f, grid);
	bbBnd.apply(0, fInd, grid);

    } // End iterations
    // -----------------END MAIN LOOP------------------

    //mpi.end();

    MPI_Finalize();

    return 0;
}

