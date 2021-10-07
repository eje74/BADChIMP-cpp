// //////////////////////////////////////////////
//
//  SUBGRID BOUNDARY CONDTIONS
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

//#include <mpi.h>

#include "../LBSOLVER.h"
#include "../IO.h"

#include "LBsubgridboundary.h"
#include "LBsubgridboundaryd.h"

//  Linear algebra package
#include "../Eigen/Dense"
#include "../Eigen/SVD"

using namespace Eigen;

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
    std::string chimpDir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/";
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


    // SETUP BULK NODES
    std::vector<int> bulkNodes = findBulkNodes(nodes);

    // *************
    // SET LB VALUES
    // *************

    // Number of iterations
    int nIterations = static_cast<int>( input["iterations"]["max"]);
    // Write interval
    int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    // Relaxation time
    lbBase_t tau = input["fluid"]["tau"][0];
    // Vector source
    VectorField<LT> bodyForce(1, 1);
    bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);
    // Driver force
    std::valarray<lbBase_t> force = bodyForce(0, 0);

    //output directory
    std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));
    std::string outDir2 = outputDir+"out"+dirNum;


    // **************
    // GEOMETRY
    // **************
    // -- New boundary condition
    ScalarField qAttribute(1, grid.size());
    VectorField<LT> surfaceNormal(1, grid.size());
    VectorField<LT> surfaceTangent(1, grid.size());

    vtklb.toAttribute("q");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        qAttribute(0, n) = vtklb.getScalarAttribute<double>();
    }
    vtklb.toAttribute("nx");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        surfaceNormal(0, 0, n) = vtklb.getScalarAttribute<double>();
    }
    vtklb.toAttribute("ny");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        surfaceNormal(0, 1, n) = vtklb.getScalarAttribute<double>();
    }
    vtklb.toAttribute("tx");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        surfaceTangent(0, 0, n) = vtklb.getScalarAttribute<double>();
    }
    vtklb.toAttribute("ty");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        surfaceTangent(0, 1, n) = vtklb.getScalarAttribute<double>();
    }


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

    // ******************
    // SETUP BOUNDARY
    // ******************
    OneNodeSubGridBndDyn<LT> fluidWallBnd(findFluidBndNodes(nodes), nodes, grid, qAttribute, surfaceNormal, surfaceTangent, rho, bodyForce, tau);
    // OneNodeSubGridBnd<LT> fluidWallBnd(findFluidBndNodes(nodes), nodes, grid, qAttribute, surfaceNormal, surfaceTangent, rho, bodyForce, tau);
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
    VTK::Output<VTK::voxel, double> output(VTK::BINARY, node_pos, outDir2, myRank, nProcs);
    output.add_file("lb_run");
    //Output output(global_dimensions, outDir2, myRank, nProcs, node_pos);
    // Add density to the output file
    VectorField<D3Q19> velIO(1, grid.size());
    output.add_variable("rho", 1, rho.get_data(), rho.get_field_index(0, bulkNodes));
    output.add_variable("vel", 3, velIO.get_data(), velIO.get_field_index(0, bulkNodes));
    // output.add_variable("rho0", 1, rho.get_data(), rho.get_field_index(0, bulkNodes));
    // output.add_variable("rho1", 1, rho.get_data(), rho.get_field_index(1, bulkNodes));
    // output.add_variable("vel", LT::nD, vel.get_data(), vel.get_field_index(0, bulkNodes));
    output.write(0);

    // Print geometry and boundary marker
    outputGeometry("lb_geo", outDir2, myRank, nProcs, nodes, grid, vtklb);

    // *********
    // MAIN LOOP
    // *********
    // For all time steps
    const std::clock_t beginTime = std::clock();
    for (int i = 0; i <= nIterations; i++) {
        // ***************
        // GLOBAL COUNTERS
        // ***************
        double rhoSumLocal = 0;
        // For all bulk nodes
        for (auto nodeNo: bulkNodes) {
            // Copy of local velocity diestirubtion
            std::valarray<lbBase_t> fNode = f(0, nodeNo);

            // MACROSCOPIC VALUES
            lbBase_t rhoNode = calcRho<LT>(fNode);
            // std::valarray<lbBase_t> velNode = calcVel<LT>(fNode, rhoNode, force);
            auto velNode = calcVel<LT>(fNode, rhoNode, force);
            // velNode = calcVel(fNode, rhoNode, force); // Kan bruke using

            rho(0, nodeNo) = rhoNode;
            rhoSumLocal += rhoNode-1;
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
            fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);

        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(0, f, grid);
        // Half way bounce back
        fluidWallBnd.applyNonSlip(0, f, nodes, grid, qAttribute, surfaceNormal, surfaceTangent, bodyForce, tau);
//        fluidWallBnd.apply(0, f, nodes, grid);
        // *************
        // WRITE TO FILE
        // *************
        double rhoSumGlobal;
        MPI_Allreduce(&rhoSumLocal, &rhoSumGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if ( ((i % nItrWrite) == 0) && (i > 0) ) {
            for (auto nn: bulkNodes) {
                velIO(0, 0, nn) = vel(0, 0, nn);
                velIO(0, 1, nn) = vel(0, 1, nn);
                velIO(0, 2, nn) = 0;
            }
            // output.write("lb_run", i);
            output.write(i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << " ( " << float( std::clock () - beginTime ) /  CLOCKS_PER_SEC << " sec)" << std::endl;
                std::cout << "Error in mass:" << rhoSumGlobal << "  " << std::endl;
                // Setup plot over line
                int nx = 1;
                for (int ny = 1; ny < 19; ++ny ) {
                    std::vector<int> pos { nx, ny};
                    int nodeNo = grid.nodeNo(pos);
//                   std::cout << "vel = " << vel(0, 0, nodeNo) << ", " << vel(0, 1, nodeNo) << std::endl;
                }
            }
        }

    } // End iterations

    MPI_Finalize();

    return 0;
}
