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
    //std::string chimpDir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/";
    std::string chimpDir = "./";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

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


    //---------------------END OF SETUP OF BOUNDARY CONDITION AND MPI---------------------

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

    lbBase_t tauDiff0 = input["diff-field"]["tauDiff"][0];
    lbBase_t tauDiff1 = input["diff-field"]["tauDiff"][1];
    lbBase_t tauDiff2 = input["diff-field"]["tauDiff"][2];

    lbBase_t H = input["Partial-Misc"]["H"];
    lbBase_t kinConst = input["Partial-Misc"]["kinConst"];

    std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));

    std::string outDir2 = outputDir+"out"+dirNum;

    //std::string outDir2 = <std::string>(input["out"]["directory"][0]); //virker ikke!

    // SET DERIVED VARAIBLES
    lbBase_t nu0Inv = 1.0 / (LT::c2 * (tau0 - 0.5));
    lbBase_t nu1Inv = 1.0 / (LT::c2 * (tau1 - 0.5));

    lbBase_t DwInv = 1.0 / (LT::c2 * (tauDiff0 - 0.5));
    lbBase_t DdwInv = 1.0 / (LT::c2 * (tauDiff1 - 0.5));
    lbBase_t DsInv = DwInv; //1.0 / (LT::c2 * (tauDiff0 - 0.5));
    lbBase_t DoilInv = 1.0 / (LT::c2 * (tauDiff2 - 0.5));
    lbBase_t Dw = (LT::c2 * (tauDiff0 - 0.5));
    lbBase_t Ddw = (LT::c2 * (tauDiff1 - 0.5));
    lbBase_t Ds = Dw;//(LT::c2 * (tauDiff0 - 0.5));
    lbBase_t Doil = (LT::c2 * (tauDiff2 - 0.5));


    // SETUP LB FIELDS
    LbField<LT> f(1, grid.size());  // LBfield
    LbField<LT> fTmp(1, grid.size());  // LBfield

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(1, grid.size()); // LBfield
    VectorField<LT> vel(1, grid.size()); // LBfield
    ScalarField eff_nu(1, grid.size()); // LBfield

    // FILL MACROSCOPIC FIELDS
    //   Fluid densities and velocity

    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 1.0;
        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
	vel(0, 1, nodeNo) = 0.03;
    }
    //---------------------END OF INPUT TO INITIALIZATION OF FIELDS---------------------


    // INITIATE LB FIELDS
    // -- fluid field
    initiateLbField(0, 0, 0, bulkNodes, rho, vel, f);  // LBinitiatefield

    // JLV
    //---------------------SETUP OUTPUT---------------------
    std::vector<std::vector<int>> node_pos;
    node_pos.reserve(bulkNodes.size());
    for (const auto& node:bulkNodes) {
        node_pos.push_back(grid.pos(node));
    }

    Output output(globalFile.dim_global(), outDir2, myRank, nProcs-1, node_pos);
    output.add_file("fluid");
    output["fluid"].add_variable("rho", rho.get_data(), rho.get_field_index(0, bulkNodes), 1);
    output["fluid"].add_variable("vel", vel.get_data(), vel.get_field_index(0, bulkNodes), LT::nD);
    output["fluid"].add_variable("eff_nu", eff_nu.get_data(), eff_nu.get_field_index(0, bulkNodes), 1);

    output.write("fluid", 0);

    //---------------------END OF SETUP OUTPUT---------------------


    // -----------------MAIN LOOP------------------
    /* Comments to main loop:
     * Calculation of cu is kept outside of calcOmega and calcDeltaOmega
     *  to avoid recalculation  of the dot-product.
     * It is therefore natural to give all scalar products as inputs to the
     *  calcOmega and calcDeltaOmega.
     * Note that we have a oneIteration function in 'LBiteration.h'. This seemed
     *   to affect the speed an gave a slow down of 10-15 %.
     */



    for (int i = 0; i <= nIterations; i++) {

        for (auto nodeNo : bulkNodes) {
            // UPDATE MACROSCOPIC DENSITIES
            lbBase_t rhoTotNode = rho(0, nodeNo) = calcRho<LT>(f(0,nodeNo));  // LBmacroscopic
        }  // End for all bulk nodes



        for (auto nodeNo: bulkNodes) {

            // Set the local total lb distribution
            std::valarray<lbBase_t> fTot = f(0, nodeNo);

            // -- total density
            lbBase_t rhoTotNode = rho(0, nodeNo);
            // -----------------------------------------------

            // -- force
            std::valarray<lbBase_t> forceNode = bodyForce(0, 0);

            // -- velocity
            std::valarray<lbBase_t> velNode = calcVel<LT>(fTot, rhoTotNode, forceNode);  // LBmacroscopic
            vel.set(0, nodeNo) = velNode;

	    std::valarray<lbBase_t> Stilde = calcShearRateTilde<LT>(fTot, rhoTotNode, velNode, forceNode, 0); // LBmacroscopic
	    /*
	    if(i == 1000 && nodeNo == 5){
	      std::cout << "Sij = [";
	      for (int d = 0; d < LT::nD*LT::nD; ++d)
		std::cout << " "<< Stilde[d];
	      std::cout << " ]"<<std::endl;
	    }
	    */
	    std::valarray<lbBase_t> StildeLowTri = calcShearRateTildeLowTri<LT>(fTot, rhoTotNode, velNode, forceNode, 0); // LBmacroscopic
	    lbBase_t StildeAbs= sqrt(2*LT::contractionLowTri(StildeLowTri, StildeLowTri));
	    if(i == 1000 && nodeNo == 5)
	      std::cout << "|Stilde| = "<<StildeAbs<<std::endl;
	    
	    lbBase_t CSmagorinsky = 0.16;
            // CALCULATE BGK COLLISION TERM
            // Mean collision time /rho_tot/\nu_tot = \sum_s \rho_s/\nu_s
            lbBase_t tau = 0.5*(tau0+sqrt(tau0*tau0+2*CSmagorinsky*CSmagorinsky*LT::c4Inv*StildeAbs/rhoTotNode));

	    eff_nu(0, nodeNo) = LT::c2Inv*tau + 0.5;

            lbBase_t uu = LT::dot(velNode, velNode);  // Square of the velocity
            std::valarray<lbBase_t> cu = LT::cDotAll(velNode);  // velocity dotted with lattice vectors
            std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fTot, tau, rhoTotNode, uu, cu);  // LBcollision
  
            // CALCULATE FORCE CORRECTION TERM
            lbBase_t  uF = LT::dot(velNode, forceNode);
            std::valarray<lbBase_t>  cF = LT::cDotAll(forceNode);
            std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);  // LBcollision


            // COLLISION AND PROPAGATION


            //-----------------------
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be

                fTmp(0, q,  grid.neighbor(q, nodeNo)) =    fTot[q]            + omegaBGK[q]           + deltaOmegaF[q];
                //-----------------------
            }
        } // End nodes

        // PRINT

        if ( (i % static_cast<int>(input["iterations"]["write"])) == 0) {

            // JLV
            output.write("fluid", i);
            // JLV
            if (myRank==0)
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
        }

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield


        // MPI Boundary
        // Hente verdier hente fra ghost
        // Sette i bulk
        mpiBoundary.communicateLbField(0, f, grid);


        // BOUNDARY CONDITIONS
        bbBnd.apply(0, f, grid);  // LBboundary



    } // End iterations
    // -----------------END MAIN LOOP------------------

    //mpi.end();

    MPI_Finalize();

    return 0;
}
