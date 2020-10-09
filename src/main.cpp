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
#define LT D3Q19



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

    // Set grid's positions
    vtklb.toPos();
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        std::vector<int> pos = vtklb.getPos<int>();
        grid.addNodePos(pos, n);
    }

    // Set node neighborhoods
    vtklb.toNeighbors();
    for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        std::vector<int> neigNodes = vtklb.getNeighbors<int>();
        grid.addNeighbors(neigNodes, n);
    }


    // Set rank for processors shared nodes
    nodes.addNodeRank(-1, 0); // Set default to -1
    for (int n = 0; n < vtklb.getNumNeigProc(); ++n) {
        int nodeRank = vtklb.getNeigRank(n);
        std::vector<int> neigNodes = vtklb.getNeigNodesNo(n);
        for (auto nodeNo: neigNodes) {
            nodes.addNodeRank(nodeRank, nodeNo);
        }
    }
 
    // Setup node types
    vtklb.toAttribute("nodetype");
    for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        int nodeType = vtklb.getScalar<int>();
        nodes.addNodeType(nodeType, n);
    }
    nodes.setupNodeType(grid);

    // Setup inlet outlet boundary markers
//    std::vector<int> boundaryMarker(grid.size());
    ScalarField boundaryMarker(1, grid.size());
    std::vector<int> bndM(grid.size());
//    std::vector<int> pressureBndNodes; // Pressure boundary
    boundaryMarker(0, 0) = 0;
    bndM[0]  = 0;
    
    vtklb.toAttribute("boundary");
    std::vector<int> inletBoundary, outletBoundary, freeSlipBoundary;
    for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
    {
        boundaryMarker(0, n) = vtklb.getScalar<int>();
        bndM[n] = boundaryMarker(0, n);
        if (bndM[n] == 1) {
            inletBoundary.push_back(n);
        } else if (bndM[n] == 2) {
            outletBoundary.push_back(n);
 //           std::cout << " " << n << " " <<std::endl;
        } else if (bndM[n] > 2) {
            freeSlipBoundary.push_back(n);
        }

    } 

    // Setup the MPI communication object.
    BndMpi<LT> mpiBoundary(myRank);
    mpiBoundary.setup(vtklb, nodes, grid);

    // *******************
    // END TEST VTKLB FILE
    // *******************


    // WALL Boundary
    std::vector<int> bulkNodes =  findBulkNodes(nodes, bndM);
    std::vector<int> fluidBndNodes = findFluidBndNodes(nodes, bndM);
    HalfWayBounceBack<LT> bounceBackBnd(fluidBndNodes, nodes, grid);
    HalfWayBounceBack<LT> freeSlipBnd(freeSlipBoundary, nodes, grid);
    InletOutlet<LT> inletBnd(inletBoundary, nodes, grid);
    InletOutlet<LT> outletBnd(outletBoundary, nodes, grid);
    // PressureBnd<LT> pBnd(pressureBndNodes, nodes, grid);


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
    VectorField<LT> vel(1, grid.size());
    // FILL MACROSCOPIC FIELDS


    // INITIATE LB FIELDS
    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 1;        
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = LT::w[q]*rho(0, nodeNo);
        }
    }

    // INIT INLET
/*    for (auto nodeNo: inletBoundary)
        rho(0, nodeNo) = 1.1;
    for (auto nodeNo: outletBoundary)
        rho(0, nodeNo) = 0.9;
*/
    // Setup pressure boundaries
/*    for (auto nodeNo: pressureBndNodes) {
        rho(0, nodeNo)  = 0.0;
        if (boundaryMarker[nodeNo] == 1)
            rho(0, nodeNo) = 1.0;
    } */

    // JLV
    //---------------------SETUP OUTPUT---------------------
    std::vector<int> globalDim(3, 1); // Set as default to 3 dimensions, as prescribed by the Output class
    for (int d = 0; d < LT::nD; ++d)
        globalDim[d] = vtklb.getGlobaDimensions(d);

    std::vector<std::vector<int>> node_pos;
    node_pos.reserve(bulkNodes.size());
    for (const auto& node:bulkNodes) {
        std::vector<int> pos(3, 0);
        for (int d=0; d < LT::nD; ++d) pos[d] = grid.pos(node, d);
        node_pos.push_back(pos);
    }
      Output output(globalDim, outputDir, myRank, nProcs, node_pos);  
      output.add_file("diff");
      output["diff"].add_variable("rho", rho.get_data(), rho.get_field_index(0, bulkNodes), 1); 
      output["diff"].add_variable("vel", vel.get_data(), vel.get_field_index(0, bulkNodes), LT::nD);

    // Write geo      
    outputGeometry("geo", outputDir, myRank, nProcs, nodes, grid, vtklb);
    
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

        if ( (i % static_cast<int>(input["iterations"]["write"])) == 0) {
            
            // JLV
            output.write("diff", i);
            // JLV
            if (myRank==0)
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
        }

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield


        // BOUNDARY CONDITIONS
        lbBase_t rhoBnd = 1.0;
        std::vector<lbBase_t> velBnd = {0, 0, 0};
        
        
        inletBnd.apply(0, f, grid, 1.1, velBnd);  
        outletBnd.apply(0, f, grid, 0.9, velBnd);
         //       pBnd.apply(0, f, grid, rho); // Pressure boundary
        bounceBackBnd.apply(0, f, grid);  // LBboundary

        freeSlipBnd.apply(0, f, grid);

        // MPI Boundary
        // Hente verdier fra ghost
        // Sette i bulk
        mpiBoundary.communicateLbField(0, f, grid);
        


    } // End iterations
    // -----------------END MAIN LOOP------------------

    MPI_Finalize();

    return 0;
}
