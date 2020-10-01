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
#include "LBfreeSlipCartesian.h"
#include "LBfreeFlowCartesian.h"
#include "LBinitiatefield.h"
#include "LBmacroscopic.h"
#include "LBnodes.h"
#include "LBsnippets.h"
#include "LButilities.h"

#include "Input.h"
#include "Output.h"

#include "LBvtk.h"
#include<algorithm> // std::max


// Pressure boundary condition
template<typename DXQY>
void applyPressureBoundary(lbBase_t rho, LbField<DXQY> &f, std::vector<int> &bndNodes, Grid<DXQY> &grid)
{
    for(auto nodeNo: bndNodes) {
        for (int q = 0; q < DXQY::nQ; ++q) {
            int neigNo = grid.neighbor(q, nodeNo);
            f(0, q, neigNo) = DXQY::w[q]*rho;
        }   
    }
}

// SET THE LATTICE TYPE
#define LT D2Q9



int main()
{
    // SETUP MPI
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // SETUP THE INPUT AND OUTPUT PATHS
    //std::string chimpDir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/";
    std::string chimpDir = "./";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

    // *******************
    // READ THE VTKLB FILE
    // *******************
    LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");

    Grid<LT> grid(vtklb.endNodeNo());
    Nodes<LT> nodes(vtklb.endNodeNo(), myRank);
    ScalarField val( 1, grid.size() );

    // Set grid's positions
    vtklb.toPos();
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        std::vector<int> pos = vtklb.getPos<int>();
        grid.addNodePos(pos, n);
    }

    // Set node neighborhoods
    vtklb.toNeighbors();
    for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        grid.addNeighbors(vtklb.getNeighbors<int>(), n);
    }

    // Set rank for processors shared nodes
    for (int n = 0; n < vtklb.getNumNeigProc(); ++n) {
        int nodeRank = vtklb.getNeigRank(n);
        std::vector<int> neigNodes = vtklb.getNeigNodesNo(n);
        for (auto nodeNo: neigNodes) {
            nodes.addNodeRank(nodeRank, nodeNo);
        }
    }

    // Set node type
    vtklb.toAttribute("nodetype");
    for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        int nodeType = vtklb.getScalar<int>();
        nodes.addNodeType(nodeType, n);
    }
    nodes.setupNodeType(grid);

    // Read the test attribute
    vtklb.toAttribute("val");
    for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        val(0, n) =  vtklb.getScalar<int>();
    }

    // Setup the MPI communication object.
    BndMpi<LT> mpiBoundary(myRank);
    mpiBoundary.setup(vtklb, nodes, grid);
    // *******************
    // END TEST VTKLB FILE
    // *******************

    // WALL Boundary
    std::vector<int> bulkNodes = findBulkNodes(nodes);
    std::vector<int> fluidBndNodes = findFluidBndNodes(nodes);
    HalfWayBounceBack<LT> bbBnd(fluidBndNodes, nodes, grid);

    //---------------------END OF SETUP OF BOUNDARY CONDITION AND MPI---------------------

    // READ INPUT FILE
    Input input(inputDir + "input.dat");
    int nIterations = static_cast<int>( input["iterations"]["max"]);

    lbBase_t tau = input["diff"]["tau"];

    std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));

    std::string outDir2 = outputDir+"out"+dirNum;


    // SETUP LB FIELDS
    LbField<LT> f(1, grid.size());  // LBfield
    LbField<LT> fTmp(1, grid.size());  // LBfield

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(1, grid.size()); // LBfield
    // FILL MACROSCOPIC FIELDS


    //---------------------END OF INPUT TO INITIALIZATION OF FIELDS---------------------


    // INITIATE LB FIELDS
    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 0;
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = LT::w[q]*rho(0, nodeNo);
        }
    }

    // JLV
    //---------------------SETUP OUTPUT---------------------
    std::vector<std::vector<int>> node_pos;
    node_pos.reserve(bulkNodes.size());
    for (const auto& node:bulkNodes) {
        node_pos.push_back(grid.pos(node));
    }


    for (int i = 0; i <= nIterations; i++) {

        for (auto nodeNo : bulkNodes) {
            // UPDATE MACROSCOPIC DENSITIES
            lbBase_t rhoTotNode = rho(0, nodeNo) = calcRho<LT>(f(0,nodeNo));  // LBmacroscopic
        }  // End for all bulk nodes

        for (auto nodeNo: bulkNodes) {
            std::valarray<lbBase_t> fNode = f(0, nodeNo);
            lbBase_t rhoNode = calcRho<LT>(fNode);
            // CALCULATE BGK COLLISION TERM
            std::valarray<lbBase_t> omegaBGK = calcOmegaPureDiff<LT>(fNode, tau, rhoNode); // LBcollision

            // COLLISION AND PROPAGATION
            //-----------------------
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) =    fNode[q]            + omegaBGK[q];
            }
        } // End nodes

        // PRINT

        if ( (i % static_cast<int>(input["iterations"]["write"])) == 0) {

            // JLV
            // output.write("fluid", i);
            // JLV
            if (myRank==0)
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
        }

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield


        // MPI Boundary
        // Hente verdier fra ghost
        // Sette i bulk
        mpiBoundary.communicateLbField(0, f, grid);


        // BOUNDARY CONDITIONS
        bbBnd.apply(0, f, grid);  // LBboundary


    } // End iterations
    // -----------------END MAIN LOOP------------------

    MPI_Finalize();

    return 0;
}
