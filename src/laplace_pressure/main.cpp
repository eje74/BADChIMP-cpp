// //////////////////////////////////////////////
//
// BADChIMP laplace_pressure
//
// For documentation see:
//    doc/documentation.pdf
// 
// //////////////////////////////////////////////

#include "../LBSOLVER.h"
#include "../IO.h"
#include "./LBdiffsuion.h"


// SET THE LATTICE TYPE
#define LT D2Q9
#define VTK_CELL VTK::pixel

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
    std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
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
    // Set bulk nodes
    std::vector<int> bulkNodes = findBulkNodes(nodes);

    // *************
    // READ FROM INPUT
    // *************
    // Number of iterations
    int nIterations = static_cast<int>( input["iterations"]["max"]);
    // Write interval
    int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    // Relaxation time
    lbBase_t tau = input["diffusion"]["tau"];

    // *************
    // SET PHYSICS
    // *************
    DiffusionSolver<LT> diffusion(tau, vtklb, nodes, grid);

/*    // TEST
    auto wallBnd{diffusion.getWallBoundary()};
    auto wallNormals{diffusion.getWallNormals()};
    auto wallNeighbors(diffusion.getWallNeighbors());
    std::cout << "WALL NODES (" << wallBnd.size() << ")" << std::endl;

    for (int n = 0; n < wallBnd.size(); ++n) {
        std::cout << " " << wallBnd.nodeNo(n) << "   ";
        std::cout << "(" << wallNormals[n][0] << ", " << wallNormals[n][1] << ") " << wallNeighbors[n]<< std::endl;
    }

    auto wallPressureBnd{diffusion.getWallPressureBoundaryNodes()};
    std::cout << "WALL PRESSURE (" << wallPressureBnd.size() << ")" << std::endl;

    for (int n = 0; n < wallPressureBnd.size(); ++n) {
        std::cout << " " << wallPressureBnd.nodeNo(n) << std::endl;
    }

    auto pressureBnd{diffusion.getPressureBoundaryNodes()};
    std::cout << "PRESSURE (" << pressureBnd.size() << ")" << std::endl;

    for (int n = 0; n < pressureBnd.size(); ++n) {
        std::cout << " " << pressureBnd.nodeNo(n) << std::endl;
    } */

    // ******************
    // MACROSCOPIC FIELDS
    // ******************
    // Density
    ScalarField rho(1, grid.size());
    // Initiate density from file
    vtklb.toAttribute("init_rho");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        rho(0, n) = vtklb.getScalarAttribute<lbBase_t>();
    }

    // ******************
    // SETUP BOUNDARY
    // ******************
    HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

    // *********
    // LB FIELDS
    // *********
    LbField<LT> f(1, grid.size()); 
    LbField<LT> fTmp(1, grid.size());
    // initiate lb distributions
    for (auto nodeNo: bulkNodes) {
        f.set(0, nodeNo) = diffusion.setF(rho(0, nodeNo));
    }

    // **********
    // OUTPUT VTK
    // **********
    VTK::Output<VTK_CELL, double> output(VTK::BINARY, grid.getNodePos(bulkNodes), outputDir, myRank, nProcs);
    output.add_file("lb_run");
    output.add_variable("rho", 1, rho.get_data(), rho.get_field_index(0, bulkNodes));

    // *********
    // MAIN LOOP
    // *********
    for (int i = 0; i <= nIterations; i++) {
        for (auto nodeNo: bulkNodes) {
            // Copy of local velocity diestirubtion
            const std::valarray<lbBase_t> fNode = f(0, nodeNo);

            // Macroscopic values
            const lbBase_t rhoNode = calcRho<LT>(fNode);

            // Save density and velocity for printing
            rho(0, nodeNo) = rhoNode;

            // BGK-collision term
            const std::valarray<lbBase_t> omegaBGK = diffusion.omegaBGK(fNode, rhoNode);

            // Collision and propagation
            fTmp.propagateTo(0, nodeNo, fNode + omegaBGK, grid);

        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(0, f, grid);
        // Half way bounce back
        // bounceBackBnd.apply(f, grid);
        diffusion.applyBoundaryCondition(0, f, grid);

        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0)  ) {
            output.write(i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
            }
            // Check mass conservation
            lbBase_t rhoTot = 0;
            for (auto nodeNo: bulkNodes) {
                rhoTot += rho(0, nodeNo);
            }
            std::cout << "total rho = " << rhoTot << std::endl;

        }
    } // End iterations

    auto pressureForcing = diffusion.getForcing(bulkNodes, f);
    VTK::Output<VTK_CELL, double> outputForce(VTK::BINARY, grid.getNodePos(bulkNodes), outputDir, myRank, nProcs);
    outputForce.add_file("forcing");
    outputForce.add_variable("force", LT::nD, pressureForcing.get_data(), pressureForcing.get_field_index(0, bulkNodes));    
    outputForce.write(0);

    MPI_Finalize();

    return 0;
}
