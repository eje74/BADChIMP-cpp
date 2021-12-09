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
    // Number of fields
    const int nFluidFields = 3;
    // Number of iterations
    int nIterations = static_cast<int>( input["iterations"]["max"]);
    // Write interval
    int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    // Relaxation time
    lbBase_t tau = input["fluid"]["tau"];
    // Body force
    VectorField<LT> bodyForce(1, 1);
    bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);

    // ******************
    // MACROSCOPIC FIELDS
    // ******************
    // Density
    ScalarField rho(nFluidFields, grid.size());
    ScalarField rhoRel(nFluidFields, grid.size());
    // Initiate density from file
    for (int fieldNo=0; fieldNo < rho.num_fields(); fieldNo++) {
        vtklb.toAttribute("init_rho_" + std::to_string(fieldNo));
        for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
            rho(fieldNo, n) = vtklb.getScalarAttribute<lbBase_t>();
        }
    }
    ScalarField rhoTot(1, grid.size());

    // Velocity
    VectorField<LT> vel(nFluidFields, grid.size());
    // Initiate velocity
    for (auto fieldNo=0; fieldNo < vel.num_fields(); ++fieldNo) {
        for (auto nodeNo: bulkNodes) {
            for (int d=0; d < LT::nD; ++d)
                vel(fieldNo, d, nodeNo) = 0.0;
        }
    }
    // ******************
    // SETUP BOUNDARY
    // ******************
    HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

    // *********
    // LB FIELDS
    // *********
    LbField<LT> f(nFluidFields, grid.size()); 
    LbField<LT> fTmp(nFluidFields, grid.size());
    // initiate lb distributions
    for (auto fieldNo=0; fieldNo < f.num_fields(); ++fieldNo) {
        for (auto nodeNo: bulkNodes) {
            for (int q = 0; q < LT::nQ; ++q) {
                f(fieldNo, q, nodeNo) = LT::w[q]*rho(fieldNo, nodeNo);
                fTmp(fieldNo, q, nodeNo) = 0;
            }
        }
    }
    // **********
    // OUTPUT VTK
    // **********
    auto node_pos = grid.getNodePos(bulkNodes); 
    auto global_dimensions = vtklb.getGlobaDimensions();
    Output output(global_dimensions, outputDir, myRank, nProcs, node_pos);
    output.add_file("lb_run");
    VectorField<D3Q19> velIO(1, grid.size());
    for (auto fieldNo=0; fieldNo < nFluidFields; ++fieldNo) {
        output["lb_run"].add_variable("rho" + std::to_string(fieldNo), rho.get_data(), rho.get_field_index(fieldNo, bulkNodes), 1);
    }
    output["lb_run"].add_variable("vel", velIO.get_data(), velIO.get_field_index(0, bulkNodes), 3);
    outputGeometry("lb_geo", outputDir, myRank, nProcs, nodes, grid, vtklb);

    // *********
    // MAIN LOOP
    // *********
    for (int i = 0; i <= nIterations; i++) {
        // Macroscopic values : rho, rhoTot and rhoRel = rho/rhoTot
        for (auto nodeNo: bulkNodes) {
            rhoTot(0, nodeNo) = 0;
            for (auto fieldNo=0; fieldNo < nFluidFields; ++fieldNo) {
                const auto fNode = f(fieldNo, nodeNo);
                rhoTot(0, nodeNo) += rho(fieldNo, nodeNo) = calcRho<LT>(fNode);
            }
            for (auto fieldNo=0; fieldNo < nFluidFields; ++fieldNo) {
                rhoRel(fieldNo, nodeNo) = rho(fieldNo, nodeNo)/rhoTot(0, nodeNo);
            }
        }
        mpiBoundary.communciateScalarField(rhoRel);

        for (auto nodeNo: bulkNodes) {
            // Cacluate gradient
            VectorField<LT> gradNode(nFluidFields, 1);            
            for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {
                gradNode.set(fieldNo, 0) = grad<LT>(rhoRel, fieldNo, nodeNo, grid);
            }

            // Calculate velocity
            // Copy of local velocity diestirubtion
            auto velNode = calcVel<LT>(f(0, nodeNo), rhoTot(0, nodeNo));
            for (auto fieldNo=1; fieldNo<nFluidFields; ++fieldNo)
                velNode += calcVel<LT>(f(fieldNo, nodeNo), rhoTot(0, nodeNo));
            velNode += 0.5*bodyForce(0, 0)/rhoTot(0, nodeNo);

            // Save density and velocity for printing
            vel.set(0, nodeNo) = velNode;
                
            // BGK-collision term
            const lbBase_t u2 = LT::dot(velNode, velNode);
            const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
            // Calculate the Guo-force correction
            const lbBase_t uF = LT::dot(velNode, bodyForce(0, 0));
            const std::valarray<lbBase_t> cF = LT::cDotAll(bodyForce(0, 0));

            for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {
                const auto fNode = f(fieldNo, nodeNo);
                const auto rhoNode = rho(fieldNo, nodeNo);
                const std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);
                
                const std::valarray<lbBase_t> deltaOmegaF = rhoRel(fieldNo, nodeNo)*calcDeltaOmegaF<LT>(tau, cu, uF, cF);

                // Collision and propagation
                fTmp.propagateTo(fieldNo, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);
            }
        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(f, grid);
        // Half way bounce back
        bounceBackBnd.apply(f, grid);

        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0)  ) {
            for (auto nn: bulkNodes) {
                velIO(0, 0, nn) = vel(0, 0, nn);
                velIO(0, 1, nn) = vel(0, 1, nn);
                velIO(0, 2, nn) = 0;
            }
            output.write("lb_run", i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
            }
        }

    } // End iterations

    MPI_Finalize();

    return 0;
}
