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
#include <algorithm>

#include "../lbsolver/LBbndmpi.h"
//#include "../lbsolver/LBboundary.h"
#include "../lbsolver/LBcollision.h"
//#include "../lbsolver/LBcollision2phase.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBgeometry.h"
#include "../lbsolver/LBgrid.h"
//#include "../lbsolver/LBhalfwaybb.h"
//#include "../lbsolver/LBfreeSlipCartesian.h"
//#include "../lbsolver/LBfreeFlowCartesian.h"
#include "../lbsolver/LBinitiatefield.h"
#include "../lbsolver/LBmacroscopic.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBsnippets.h"
#include "../lbsolver/LButilities.h"
//#include "../lbsolver/LBpressurebnd.h"

#include "../io/Input.h"
#include "../io/Output.h"

#include "../lbsolver/LBvtk.h"

//#include "../src/LBbounceback.h"
//#include "../src/LBfreeslipsolid.h"

//#include "../src/LBinletoutlet.h"

#include "LBturbineforce.h"
#include "LBturbineboundary.h"


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
    // -- Read position data
    VectorField<LT> position(1, grid.size());
    std::vector<std::string> attributeName {"pos_x", "pos_y", "pos_z"};
    for (int n = 0; n < 3; ++n) {
        vtklb.toAttribute(attributeName[n]);
        for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
            position(0, n, nodeNo) = vtklb.getScalar<lbBase_t>();
        }
        std::cout << "READ ATTRIBUTE " << attributeName[n] << " ON RANK " << myRank << std::endl;
    }

    // *************
    // SET LB VALUES
    // *************
    // Relaxation time
    lbBase_t tau = 1.0;
    // Driver force
    std::valarray<lbBase_t> force {0.0, 0.0, 0.0};
    // Rotation speed
    lbBase_t omega_x = 2*3.14159/2000.0;

    // ********
    // SETUP BOUNDARY
    // ********
    // -- Read boundar marker attritute
    std::vector<int> inletNodes, outletNodes, fluidWallNodes;
    std::vector<int> solidWallNodes;
    std::vector<int> bulkNodes;

    vtklb.toAttribute("boundary_marker");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        auto val = vtklb.getScalar<int>();
        switch (val) {
        case 3:  // Outlet
            outletNodes.push_back(nodeNo);
            break;
        case 2:  // Inlet
            inletNodes.push_back(nodeNo);
            break;
        case 1:  // fluid wall
            if (nodes.isFluidBoundary(nodeNo) && nodes.isMyRank(nodeNo)) {
                fluidWallNodes.push_back(nodeNo);
            } else if ( nodes.isFluid(nodeNo) && nodes.isMyRank(nodeNo) ) { // is marked as a standard fluid or inlet/outlet
                bulkNodes.push_back(nodeNo);
            }
            break;
        case 0:  // Bulk node
            if ( nodes.isFluid(nodeNo) && nodes.isMyRank(nodeNo) ) { // is marked as a standard fluid or inlet/outlet
                bulkNodes.push_back(nodeNo);
            }
            if ( nodes.isFluidBoundary(nodeNo) && nodes.isMyRank(nodeNo) ) { // is marked as a standard fluid or inlet/outlet
                solidWallNodes.push_back(nodeNo);
            }
            break;
        default :
            break;
        }
    }

    // -- Inlet boundary
    InletBoundary<LT> inletBnd(omega_x, position, inletNodes, nodes, grid);
    // -- Outlet boundary
    OutletBoundary<LT> outletBnd(omega_x, position, outletNodes, nodes, grid);
    // -- Free Slip Wall
    FreeBoundary<LT> freeBnd(omega_x, position, fluidWallNodes, nodes, grid);
    // -- Solid wall
    HalfWayBounceBack<LT> wallBnd(solidWallNodes, nodes, grid);

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
        std::valarray<lbBase_t> velInit {0, 0, 0};
        std::valarray<lbBase_t> tmpvel = velInert2Rot(omega_x, position(0, nodeNo), velInit);
        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = tmpvel[d];
        if ( fabs(vel(0, 0, nodeNo)) > 0) {
            std::cout << nodeNo << " : " << vel(0, 0, nodeNo) << std::endl;
        }
    }

    // *********
    // LB FIELDS
    // *********
    LbField<LT> f(1, grid.size());  // LBfield
    LbField<LT> fTmp(1, grid.size());  // LBfield
    // initieate values
    for (auto nodeNo: bulkNodes) {
        lbBase_t u2 = LT::dot(vel(0, nodeNo), vel(0, nodeNo));
        std::valarray<lbBase_t> cu = LT::cDotAll(vel(0, nodeNo));
        auto feq = [&](int q, lbBase_t rho) -> lbBase_t {return rho * LT::w[q]*(1.0 + LT::c2Inv*cu[q] + LT::c4Inv0_5*(cu[q]*cu[q] - LT ::c2*u2));};
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = feq(q, rho(0, nodeNo));
        }
    }

    // **********
    // OUTPUT VTK
    // **********
    VectorField<LT> velInert(1, grid.size());
    Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs);
    output.add_file("lb_run");
    output.add_variables({"rho"}, rho);
    output.add_variables({"vel"}, vel);
    // auto node_pos = grid.getNodePos(bulkNodes); // Need a named variable as Outputs constructor takes a reference as input
    // auto global_dimensions = vtklb.getGlobaDimensions();
    // // Setup output file
    // Output output(global_dimensions, outputDir, myRank, nProcs, node_pos);
    // Add density to the output file
    // output.add_file("lb_run");
    // output["lb_run"].add_variable("rho", rho.get_data(), rho.get_field_index(0, bulkNodes), 1);
    // output["lb_run"].add_variable("vel", vel.get_data(), vel.get_field_index(0, bulkNodes), LT::nD);
    // output.write("lb_run", 0);
    // Print geometry and boundary marker
    // outputGeometry("lb_geo", outputDir, myRank, nProcs, nodes, grid, vtklb);
    // Print the flud wall nodes

    std::vector<int> wallMarker(grid.size(), 0);
    int cnt = 0;
    for (auto nodeNo: fluidWallNodes) {
        wallMarker[nodeNo] = -1;
        for (int beta = 0; beta < 18; ++beta)
        {
            int neigNo = grid.neighbor(beta, nodeNo);
            if (nodes.isFluid(neigNo) && !nodes.isFluidBoundary(neigNo))
                wallMarker[nodeNo] = 1;
        }
        
        cnt += 1;
        }
    for (auto nodeNo: inletNodes) {
        wallMarker[nodeNo] = 2;
    }

    for (auto nodeNo: outletNodes) {
        wallMarker[nodeNo] = 3;
    }
    for (auto nodeNo: solidWallNodes) {
        wallMarker[nodeNo] = 4;
    }
    outputStdVector("wall_nodes", wallMarker, outputDir, myRank, nProcs, grid, vtklb);

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
            // -- use old values for the velocitis
            std::valarray<lbBase_t> forceNode = force + rotatingForce(omega_x, rhoNode, position(0, nodeNo), vel(0, nodeNo));
            std::valarray<lbBase_t> velNode = calcVel<LT>(fNode, rhoNode, force);
            rho(0, nodeNo) = rhoNode;
            for (int d = 0; d < LT::nD; ++d)
                vel(0, d, nodeNo) = velNode[d];
            // BGK COLLISION TERM
            // SRT
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
        // Inlet boundary
        inletBnd.apply(tau, force, 1.0, vel, 0, f, grid);
        // Outlet boundary
        outletBnd.apply(tau, force, 1.0, vel, 0, f, grid);
        // Free boundary
        freeBnd.apply(tau, force, 1.0, vel, 0, f, grid);
        // Wall boundary
        wallBnd.apply(0, f, grid);
        // Mpi
        mpiBoundary.communicateLbField(0, f, grid);

        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0) && (i > 0) ) {
            for (auto nodeNo: bulkNodes) {
                velInert(0, nodeNo) = velRot2Inert(omega_x, position(0, nodeNo), vel(0, nodeNo));
            }
            output.write(i);
            if (myRank==0)
                std::cout << "PLOT AT ITERATION : " << i << " ( " << float( std::clock () - beginTime ) /  CLOCKS_PER_SEC << " sec)" << std::endl;
        }
    } // End iterations

    MPI_Finalize();

    return 0;
}
