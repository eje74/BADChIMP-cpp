// //////////////////////////////////////////////
//
// BADChIMP std_case
//
// For documentation see:
//    doc/documentation.pdf
// 
// //////////////////////////////////////////////

#include "../LBSOLVER.h"
#include "../IO.h"

// SET THE LATTICE TYPE
//#define LT D2Q9
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
    //std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
    std::string chimpDir = "/Users/janlv/github/BADCHiMP-cpp/";
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
    // int nIterations = static_cast<int>( input["iterations"]["max"]);
    int nIterations = input["iterations"]["max"];
    // Write interval
    // int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    int nItrWrite = input["iterations"]["write"];
    // Relaxation time
    lbBase_t tau = input["fluid"]["tau"];
    // Body force
    // VectorField<LT> bodyForce(1, 1);
    //bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);
    VectorField<LT> bodyForce(1, 1, input["fluid"]["bodyforce"]);

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
    
    // Velocity
    VectorField<LT> vel(1, grid.size());
    // Initiate velocity
    for (auto nodeNo: bulkNodes) {
        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
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
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = LT::w[q]*rho(0, nodeNo);
        }
    }

    // **********
    // OUTPUT VTK
    // **********
    Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs);
    output.add_file("lb_run");
    output.add_scalar_variables({"rho"}, {rho});
    output.add_vector_variables({"vel"}, {vel});

    // VTK::Output<VTK_CELL, double> output(VTK::BINARY, grid.getNodePos(bulkNodes), outputDir, myRank, nProcs);
    // output.add_file("lb_run");
    // output.add_variable("rho", 1, rho.get_data(), rho.get_field_index(0, bulkNodes));
    // output.add_variable("vel", LT::nD, vel.get_data(), vel.get_field_index(0, bulkNodes));

    // *********
    // MAIN LOOP
    // *********
    for (int i = 0; i <= nIterations; i++) {
        for (auto nodeNo: bulkNodes) {
            // Copy of local velocity diestirubtion
            const std::valarray<lbBase_t> fNode = f(0, nodeNo);

            // Macroscopic values
            const lbBase_t rhoNode = calcRho<LT>(fNode);
            const auto velNode = calcVel<LT>(fNode, rhoNode, bodyForce(0, 0));

            // Save density and velocity for printing
            rho(0, nodeNo) = rhoNode;
            vel.set(0, nodeNo) = velNode;
                
            // BGK-collision term
            const lbBase_t u2 = LT::dot(velNode, velNode);
            const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
            const std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);
            
            // Calculate the Guo-force correction
            const lbBase_t uF = LT::dot(velNode, bodyForce(0, 0));
            const std::valarray<lbBase_t> cF = LT::cDotAll(bodyForce(0, 0));
            const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);

            // Collision and propagation
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
        bounceBackBnd.apply(f, grid);

        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0)  ) {
            output.write(i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
            }
        }

    } // End iterations

    MPI_Finalize();

    return 0;
}
