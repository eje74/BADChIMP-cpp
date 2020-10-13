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
#include<algorithm> // std::max


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

    // ***********************
    // SETUP GRID AND GEOMETRY
    // ***********************
    LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
    Grid<LT> grid(vtklb);
    Nodes<LT> nodes(vtklb, grid);
    BndMpi<LT> mpiBoundary(vtklb, nodes, grid);



    // ****************
    // SETUP BOUNDARIES
    // ****************
    // *** inlet, outlet and free slip ***
    //  Boundaries set in the vtklb file
    std::vector<int> boundaryMarker;
    boundaryMarker.reserve(grid.size());
    boundaryMarker.push_back(0);
    
    vtklb.toAttribute("boundary");
    std::vector<int> inletBoundary, outletBoundary, freeSlipBoundary;
    for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
    {
        int val = vtklb.getScalar<int>();
        boundaryMarker.push_back(val);
        if (val == 1) {
            freeSlipBoundary.push_back(n);
        } else if (val == 2) {
            inletBoundary.push_back(n);
        } else if (val == 3) {
            outletBoundary.push_back(n);
        }
    } 
    // *** wall boundaries ***
    
    
    // WALL Boundary
    std::vector<int> bulkNodes =  findBulkNodes(nodes, boundaryMarker);
    std::vector<int> fluidBndNodes = findFluidBndNodes(nodes, boundaryMarker);
    HalfWayBounceBack<LT> bounceBackBnd(fluidBndNodes, nodes, grid);
    HalfWayBounceBack<LT> freeSlipBnd(freeSlipBoundary, nodes, grid);
    InletOutlet<LT> inletBnd(inletBoundary, nodes, grid);
    InletOutlet<LT> outletBnd(outletBoundary, nodes, grid);
    // PressureBnd<LT> pBnd(pressureBndNodes, nodes, grid);


    //---------------------END OF SETUP OF BOUNDARY CONDITION AND MPI---------------------

    // SETUP the boundary markers



    // READ INPUT FILE
    Input input(inputDir + "input.dat");

    int nIterations = static_cast<int>( input["iterations"]["max"]);
    int nItrWrite = static_cast<int>( input["iterations"]["write"]);
        

    lbBase_t tau = input["diff"]["tau"];

    std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));

    std::string outDir2 = outputDir+"out"+dirNum;



    // SETUP LB FIELDS
    LbField<LT> f(1, grid.size());  // LBfield
    LbField<LT> fTmp(1, grid.size());  // LBfield

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(1, grid.size()); // LBfield
    VectorField<LT> vel(1, grid.size());
    // FILL MACROSCOPIC FIELDS


    // INITIATE LB FIELDS
    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 1;        
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = LT::w[q]*rho(0, nodeNo);
        }
    }


    // JLV
    //---------------------SETUP OUTPUT---------------------
    auto node_pos = grid.getNodePos(bulkNodes); // Need a named variable as Outputs constructor takes a reference as input
 //   std::cout << "node pos size " << node_pos.size() << " " << bulkNodes.size() << std::endl;

/*    std::cout << "SIZE = " << node_pos.size() << " " << bulkNodes.size() << std::endl;
    for (int c = 0; c < node_pos.size(); ++c)
    {
        std::cout << grid.pos(bulkNodes[c], 0) << " " << grid.pos(bulkNodes[c], 1) << " | ";
        std::cout << node_pos[c][0] << " " << node_pos[c][1] << " " << node_pos[c][2] << std::endl;
    } */

/*    std::vector<std::vector<int>> node_pos;
    node_pos.reserve(bulkNodes.size());
    for (auto nodeNo: bulkNodes)
    {
        std::vector<int> tmp(3, 1);
        for (int i=0; i < LT::nD; ++i)
            tmp[i] = grid.pos(nodeNo, i);
        node_pos.push_back(tmp);
    } */
    //    return nodePos;
    
    auto global_dimensions = vtklb.getGlobaDimensions();


    
    Output output(global_dimensions, outputDir, myRank, nProcs, node_pos);
    output.add_file("diff");
    output["diff"].add_variable("rho", rho.get_data(), rho.get_field_index(0, bulkNodes), 1);
//    output["diff"].add_variable("vel", vel.get_data(), vel.get_field_index(0, bulkNodes), LT::nD);
    // Write geo
    outputGeometry("geo", outputDir, myRank, nProcs, nodes, grid, vtklb);
    outputStdVector("boundary", boundaryMarker, outputDir, myRank, nProcs, grid, vtklb);
    
    // Write boundary
//    std::vector<std::vector<int>> node_pos_all = grid.getNodePos(vtklb.beginNodeNo(), vtklb.endNodeNo());
//    Output bndout(vtklb.getGlobaDimensions(), outputDir, myRank, nProcs, node_pos_all);
   
    
    for (int i = 0; i <= nIterations; i++) {

        for (auto nodeNo: bulkNodes) {
            std::valarray<lbBase_t> fNode = f(0, nodeNo);
            lbBase_t rhoNode = calcRho<LT>(fNode);
            std::valarray<lbBase_t> velNode = calcVel<LT>(fNode, rhoNode);
            rho(0, nodeNo) = rhoNode;
            vel(0, nodeNo) = velNode;
            // CALCULATE BGK COLLISION TERM
            lbBase_t u2 = LT::dot(velNode, velNode);
            std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
            
            std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);

            // COLLISION AND PROPAGATION
            //-----------------------
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) =    fNode[q]            + omegaBGK[q];
            }
        } // End nodes

        // PRINT

        if ( (i % nItrWrite) == 0) {
            
            // JLV
            output.write("diff", i);
            // JLV
 //           for (auto nodeNo: bulkNodes) {
 //               std::cout << rho(0, nodeNo) << " (" <<  myRank << ") " << std::endl;
 //           }
            if (myRank==0)
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
        }

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield


        // BOUNDARY CONDITIONS
 //       mpiBoundary.communicateLbField(0, f, grid);
        
        
         //       pBnd.apply(0, f, grid, rho); // Pressure boundary
 //       bounceBackBnd.apply(0, f, grid);  // LBboundary
 //       freeSlipBnd.apply(0, f, grid);

        // MPI Boundary
        // Hente verdier fra ghost
        // Sette i bulk
        lbBase_t rhoBnd = 1.0;
        std::vector<lbBase_t> velBnd = {0, 0, 0};
//        inletBnd.apply(0, f, grid, 1.0, velBnd);  
//        outletBnd.apply(0, f, grid, 1.0, velBnd);
        


    } // End iterations
    // -----------------END MAIN LOOP------------------

    MPI_Finalize();

    return 0;
}
