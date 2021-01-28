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
#include <ctime>

#include "../lbsolver/LBbndmpi.h"
#include "../lbsolver/LBboundary.h"
#include "../lbsolver/LBcollision.h"
#include "../lbsolver/LBcollision2phase.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBgeometry.h"
#include "../lbsolver/LBgrid.h"
#include "../lbsolver/LBhalfwaybb.h"
#include "../lbsolver/LBfreeSlipCartesian.h"
#include "../lbsolver/LBfreeFlowCartesian.h"
#include "../lbsolver/LBinitiatefield.h"
#include "../lbsolver/LBmacroscopic.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBsnippets.h"
#include "../lbsolver/LButilities.h"
#include "../lbsolver/LBpressurebnd.h"

#include "../io/Input.h"
#include "../io/Output.h"

#include "../lbsolver/LBvtk.h"

#include "../lbsolver/LBbounceback.h"
#include "../lbsolver/LBfreeslipsolid.h"

#include "../lbsolver/LBinletoutlet.h"

#include<algorithm> // std::max


// SET THE LATTICE TYPE
#define LT D3Q19



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
    std::vector<int>  freeSlipRightNodes, freeSlipLeftNodes, freeSlipBottomNodes, freeSlipTopNodes;
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        int val = vtklb.getScalar<int>();
        boundaryMarker[nodeNo] = val;
        if (val == 3) {
            freeSlipLeftNodes.push_back(nodeNo);
        } else if (val == 4) {
            freeSlipRightNodes.push_back(nodeNo);
        } else if (val == 5) {
            freeSlipBottomNodes.push_back(nodeNo);
        } else if (val == 6) {
            freeSlipTopNodes.push_back(nodeNo);
        }
    }


    std::vector<int> normalLeft = {0, -1, 0};
    SolidFreeSlip<LT> freeSlipLeftBnd(normalLeft, freeSlipLeftNodes, nodes, grid);
    std::vector<int> normalRight = {0, 1, 0};
    SolidFreeSlip<LT> freeSlipRightBnd(normalRight, freeSlipRightNodes, nodes, grid);
    std::vector<int> normalBottom = {0, 0, -1};
    SolidFreeSlip<LT> freeSlipBottomBnd(normalBottom, freeSlipBottomNodes, nodes, grid);
    std::vector<int> normalTop = {0, 0, 1};
    SolidFreeSlip<LT> freeSlipTopBnd(normalTop, freeSlipTopNodes, nodes, grid);

    // Bounce back boundary
    std::vector<int> fluidWallNodes = findFluidBndNodes(nodes);
    HalfWayBounceBack<LT> fluidWallBnd(fluidWallNodes, nodes, grid);

    // *****************
    // SETUP BULK NODES
    // *****************
    std::vector<int> bulkNodes;
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        if ( nodes.isFluid(nodeNo) &&       // is a fluid
             nodes.isMyRank(nodeNo) &&      // is on this rank
             (boundaryMarker[nodeNo] < 3) ) { // is marked as a standard fluid or inlet/outlet
            bulkNodes.push_back(nodeNo);
        }
    }

    // *************
    // SET LB VALUES
    // *************
    // Relaxation time
    lbBase_t tau = input["diff"]["tau"];
    // Driver force
    std::valarray<lbBase_t> force = {0.000005, 0.0};


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
    output["diff"].add_variable("vel", vel.get_data(), vel.get_field_index(0, bulkNodes), LT::nD);

    // Print geometry and boundary marker
    outputGeometry("geo", outputDir, myRank, nProcs, nodes, grid, vtklb);
    outputStdVector("boundary", boundaryMarker, outputDir, myRank, nProcs, grid, vtklb);
    // Print smoothed geo
    std::vector<int> smoothgeo;
    vtklb.toAttribute("smoothGeo");
    smoothgeo.push_back(0);
    for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
        smoothgeo.push_back(vtklb.getScalar<int>());
    outputStdVector("smoothgeo", smoothgeo, outputDir, myRank, nProcs, grid, vtklb);
    

    // *********
    // MAIN LOOP
    // *********
    // Number of iterations
    int nIterations = static_cast<int>( input["iterations"]["max"]);
    // Write interval
    int nItrWrite = static_cast<int>( input["iterations"]["write"]);

    // For all time steps
    const std::clock_t beginTime = std::clock();
    for (int i = 0; i <= nIterations; i++) {
        // For all bulk nodes
        for (auto nodeNo: bulkNodes) {
            // Copy of local velocity diestirubtion
            std::valarray<lbBase_t> fNode = f(0, nodeNo);

            // MACROSCOPIC VALUES
            lbBase_t rhoNode = calcRho<LT>(fNode);
            std::valarray<lbBase_t> velNode = calcVel<LT>(fNode, rhoNode, force);
            rho(0, nodeNo) = rhoNode;
            for (int d = 0; d < LT::nD; ++d)
                vel(0, d, nodeNo) = velNode[d];
            // BGK COLLISION TERM
            // Simple SRT
            lbBase_t u2 = LT::dot(velNode, velNode);
            std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
            std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);
            lbBase_t uF = LT::dot(velNode, force);
            std::valarray<lbBase_t> cF = LT::cDotAll(force);
            std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);



            // COLLISION AND PROPAGATION
            for (int q = 0; q < LT::nQ; ++q) {
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = fNode[q]  + omegaBGK[q] + deltaOmegaF[q];
            }
        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(0, f, grid);
        // Free slip
        freeSlipLeftBnd.apply(0, f, grid);
        freeSlipRightBnd.apply(0, f, grid);
        freeSlipBottomBnd.apply(0, f, grid);
        freeSlipTopBnd.apply(0, f, grid);
        // Half way bounce back
        fluidWallBnd.apply(0, f, grid);


        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0) && (i > 0) ) {
            output.write("diff", i);
            if (myRank==0)
                std::cout << "PLOT AT ITERATION : " << i << " ( " << float( std::clock () - beginTime ) /  CLOCKS_PER_SEC << " sec)" << std::endl;
        }
    } // End iterations
 
    MPI_Finalize();

    return 0;
}
