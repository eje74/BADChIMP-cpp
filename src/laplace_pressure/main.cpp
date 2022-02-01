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
    // int nIterations = static_cast<int>( input["iterations"]["max"]);
    int nIterations = input["iterations"]["max"];
    // Write interval
    // int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    int nItrWrite = input["iterations"]["write"];
    // Relaxation time
    lbBase_t tau = input["diffusion"]["tau"];

    DiffusionSolver<LT> diffusion(tau, vtklb, nodes, grid);
    // *************
    // FIND NUMBER OF FIELDS 
    // *************
    int numFieldsRank = diffusion.maxBoundaryIndicator()-1;
    int numFields;
    MPI_Allreduce(&numFieldsRank, &numFields, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    // ******************
    // MACROSCOPIC FIELDS
    // ******************
    // Density
    ScalarField rho(numFields, grid.size());
    // Initiate density from file
    vtklb.toAttribute("init_rho");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        auto rhoVal = vtklb.getScalarAttribute<lbBase_t>();
        for (int i=0; i < numFields; ++i)
            rho(i, n) = rhoVal;
    }

    // *********
    // LB FIELDS
    // *********
    LbField<LT> f(numFields, grid.size()); 
    LbField<LT> fTmp(numFields, grid.size());
    // initiate lb distributions
    for (auto nodeNo: bulkNodes) {
        for (int i=0; i < numFields; ++i)
            f.set(i, nodeNo) = diffusion.setF(rho(0, nodeNo));
    }

    // **********
    // OUTPUT VTK
    // **********
    Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs);
    output.add_file("lb_run");
    output.add_variable("rho", rho);


    // *********
    // MAIN LOOP
    // *********
    for (int i = 0; i <= nIterations; i++) {
        for (int fieldNo = 0; fieldNo < numFields; ++fieldNo) {
            for (auto nodeNo: bulkNodes) {
                const std::valarray<lbBase_t> fHat = diffusion.collision(f(fieldNo, nodeNo), rho(fieldNo, nodeNo));
                fTmp.propagateTo(fieldNo, nodeNo, fHat, grid);
            } // End nodes
        }
        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************

        // Mpi send-recive
        mpiBoundary.communicateLbField(f, grid);
        // Bondary condition

        diffusion.applyBoundaryCondition(f, grid);

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


    //----------------------------------------------------------------------------------- Write forces and pressures to file
    ScalarField psiTot(1, grid.size());
    VectorField<LT> jTot(1, grid.size());

    for (auto &nodeNo: bulkNodes) {
        psiTot(0, nodeNo) = 1;
        jTot.set(0, nodeNo) = 0; 
    } 

    for (int fieldNum=0; fieldNum < numFields; ++fieldNum) {
        VectorField<LT> jVec = diffusion.getForcing(fieldNum, f, bulkNodes);
        ScalarField psi(1, grid.size());
        for (auto & nodeNo: bulkNodes) {
            psi(0, nodeNo) = rho(fieldNum, nodeNo);
            psiTot(0, nodeNo) -= psi(0, nodeNo);
            jTot.set(0, nodeNo) -=  jVec(0, nodeNo);
        }
        std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(fieldNum);
        jVec.writeToFile(mpiDir + filename);
        psi.writeToFile(mpiDir + filename);
    }
    std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(numFields);
    jTot.writeToFile(mpiDir + filename);
    psiTot.writeToFile(mpiDir + filename);


    VectorField<LT> jRead(numFields + 1, grid.size());
    ScalarField psiRead(numFields + 1, grid.size());
    for (int fieldNum=0; fieldNum < (numFields+1); ++fieldNum) {
        VectorField<LT> jTmp(1, grid.size());
        ScalarField psiTmp(1, grid.size());
        std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(fieldNum);
        jTmp.readFromFile(mpiDir + filename);
        psiTmp.readFromFile(mpiDir + filename);
        for (auto & nodeNo: bulkNodes) {
            jRead.set(fieldNum, nodeNo) = jTmp(0, nodeNo);
            psiRead(fieldNum, nodeNo) = psiTmp(0, nodeNo);
        }
    }
    Output<LT> outputForce(grid, bulkNodes, outputDir, myRank, nProcs);
    outputForce.add_file("forcing");
    outputForce.add_variables({"force", "pressure"}, {jRead, psiRead});
    outputForce.write(0); 
    // VTK::Output<VTK_CELL, double> outputForce(VTK::BINARY, grid.getNodePos(bulkNodes), outputDir, myRank, nProcs);
    // outputForce.add_file("forcing");
    // for (int fieldNum=0; fieldNum < (numFields+1); ++fieldNum) {
    //     outputForce.add_variable("force" + std::to_string(fieldNum), LT::nD, jRead.get_data(), jRead.get_field_index(fieldNum, bulkNodes));    
    //     outputForce.add_variable("pressure" + std::to_string(fieldNum), 1, psiRead.get_data(), psiRead.get_field_index(fieldNum, bulkNodes));
    // }
    // outputForce.write(0); 

    MPI_Finalize();

    return 0;
}
