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
#include "LBpressurebnd.h"

#include "Input.h"
#include "Output.h"

#include "LBvtk.h"

#include "LBbounceback.h"
#include<algorithm> // std::max


// SET THE LATTICE TYPE
#define LT D2Q9


int main()
{
    // *********
    // SETUP MPI
    // *********    
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // ********************************
    // SETUP THE INPUT AND OUTPUT PATHS
    // ********************************
    std::string chimpDir = "./";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

    // ***********************
    // SETUP GRID AND GEOMETRY
    // ***********************
    Input input(inputDir + "input.dat");
    LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
    Grid<LT> grid(vtklb);
    Nodes<LT> nodes(vtklb, grid);
    BndMpi<LT> mpiBoundary(vtklb, nodes, grid);

    // ****************
    // SETUP BOUNDARIES
    // ****************
    // Boundary markers:
    std::vector<int> boundaryMarker(grid.size());
    //   read from file
    vtklb.toAttribute("boundary");
    std::vector<int> inletNodes, outletNodes, freeSlipNodes, solidBounceBackNodes;
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo)
    {
        int val = vtklb.getScalar<int>();
        boundaryMarker[nodeNo] = val;
        if (val == 1) {  // Free slip
            freeSlipNodes.push_back(nodeNo);
            solidBounceBackNodes.push_back(nodeNo);
        } else if (val == 2) { // Inlet
            inletNodes.push_back(nodeNo);
            solidBounceBackNodes.push_back(nodeNo);
        } else if (val == 3) { // Outlet
            outletNodes.push_back(nodeNo);
            solidBounceBackNodes.push_back(nodeNo);
        }
    } 
    // Solid bounce back boundary
    SolidBounceBack<LT> solidBounceBackBnd(solidBounceBackNodes, nodes, grid);
 
    // *****************
    // SETUP BULK NODES
    // *****************
    std::vector<int> bulkNodes;
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        if ( nodes.isFluid(nodeNo) &&       // is a fluid
             nodes.isMyRank(nodeNo) &&      // is on this rank
            (boundaryMarker[nodeNo] == 0) ) // is marked as a standard fluid             
            {
                bulkNodes.push_back(nodeNo);
            }
    }

    // *************
    // SET LB VALUES
    // *************
    // Relaxation time
    lbBase_t tau = input["diff"]["tau"];

    // ******************
    // MACROSCOPIC FIELDS
    // ******************    
    // Density
    ScalarField rho(1, grid.size());
    // Velocity
    VectorField<LT> vel(1, grid.size());
    // Initiate values
    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 1.0;
        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }

    // *********
    // LB FIELDS
    // *********    
    LbField<LT> f(1, grid.size());  // LBfield
    LbField<LT> fTmp(1, grid.size());  // LBfield
    // initieate values
    for (auto nodeNo: bulkNodes) {
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = LT::w[q]*rho(0, nodeNo);
        }
    }
        
    // **********
    // OUTPUT VTK
    // **********
    auto node_pos = grid.getNodePos(bulkNodes); // Need a named variable as Outputs constructor takes a reference as input    
    auto global_dimensions = vtklb.getGlobaDimensions();
    // Setup output file
    Output output(global_dimensions, outputDir, myRank, nProcs, node_pos);
    output.add_file("diff");
    // Add density to the output file
    output["diff"].add_variable("rho", rho.get_data(), rho.get_field_index(0, bulkNodes), 1);
    // Print geometry and boundary marker
    outputGeometry("geo", outputDir, myRank, nProcs, nodes, grid, vtklb);
    outputStdVector("boundary", boundaryMarker, outputDir, myRank, nProcs, grid, vtklb);
 
    // *********
    // MAIN LOOP
    // *********    
    // Number of iterations
    int nIterations = static_cast<int>( input["iterations"]["max"]);
    // Write interval
    int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    // For all time steps
    for (int i = 0; i <= nIterations; i++) {
        // For all bulk nodes     
        for (auto nodeNo: bulkNodes) {
            // Copy of local velocity diestirubtion
            std::valarray<lbBase_t> fNode = f(0, nodeNo);
            // MACROSCOPIC VALUES
            lbBase_t rhoNode = calcRho<LT>(fNode);
            std::valarray<lbBase_t> velNode = calcVel<LT>(fNode, rhoNode);
            rho(0, nodeNo) = rhoNode;
            vel(0, nodeNo) = velNode;
            // BGK COLLISION TERM
            // Simple SRT
            lbBase_t u2 = LT::dot(velNode, velNode);
            std::valarray<lbBase_t> cu = LT::cDotAll(velNode);            
            std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);
            // COLLISION AND PROPAGATION
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) =    fNode[q]            + omegaBGK[q];
            }
        } // End nodes 
       
        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield
        
        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(0, f, grid);        
        // Solid bounce back boundary
        solidBounceBackBnd.apply(0, f, grid);        
    
        // *************
        // WRITE TO FILE
        // *************
        if ( (i % nItrWrite) == 0) {
            output.write("diff", i);
            if (myRank==0)
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
        }
    } // End iterations

    MPI_Finalize();

    return 0;
}
