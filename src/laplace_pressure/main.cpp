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
#include "./LBdiffusion.h"

#include <random>

// SET THE LATTICE TYPE
#define LT D2Q9
//#define LT D3Q19
template <typename Lattice, typename T=double, int FMT=VTK::BINARY, typename CELL=VTK::voxel>
using Output = LBOutputUnstructured<Lattice, T, FMT, CELL>;

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
    // std::string chimpDir = "./";
    std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";
    std::string laplacePressureDir = mpiDir + "uturn/";

    // ***********************
    // SETUP GRID AND GEOMETRY
    // ***********************
    Input input(inputDir + "input.dat");
    LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
    Grid<LT> grid(vtklb);
    Nodes<LT> nodes(vtklb, grid);
    BndMpi<LT> mpiBoundary(vtklb, nodes, grid);

    // ---------------------------------------------------------------------------------- Read signed distance
    // need to be defined before pressure.
    auto sd = readSignedDistance("signed_distance", vtklb, nodes, grid);

    // ---------------------------------------------------------------------------------- Add pressure boundary nodes
    int maxBoundaryIndicator = addPressureBoundary("pressure_boundary", vtklb, nodes, grid);

    std::cout << "maxBoundaryIndicator = " << maxBoundaryIndicator << std::endl;

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

    // ---------------------------------------------------------------------------------- Setup pressure boundary
    std::vector<int> bulkNodes = findBulkNodes<LT>(nodes);
    laplaceBoundary<LT> bnd(bulkNodes, sd, nodes, grid);

    // setup the pressure boundary information (normal vector, distance to boundary)

    std::vector<int> pressureNodes;
    std::vector<lbBase_t> pressureNodesS;
    std::vector<std::valarray<lbBase_t>> pressureNodesNorm;

    for (const auto &n : bulkNodes)
    {
        if (nodes.getTag(n) == 1)
        {
            pressureNodes.push_back(n);
            pressureNodesS.push_back(0.0);
            pressureNodesNorm.push_back((std::initializer_list<lbBase_t>){1, 0});
        }
        else if (nodes.getTag(n) == 2)
        {
            pressureNodes.push_back(n);
            pressureNodesS.push_back(0.0);
            pressureNodesNorm.push_back((std::initializer_list<lbBase_t>){1, 0});
        }
    }

    bnd.addPressureNodeData(pressureNodes, pressureNodesS, pressureNodesNorm);

    // Set bulk nodes
    const int numFields = maxBoundaryIndicator - 1;

    // std::vector<int> bulkNodes = diffusion.findBulkNodes(vtklb, nodes);
    ScalarField rho(numFields, grid.size());
    ScalarField rhoOld(1, grid.size());

    // Initiate density from file
    // vtklb.toAttribute("init_rho");
    for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
    {
        // auto rhoVal = vtklb.getScalarAttribute<lbBase_t>();
        for (int fieldNum = 0; fieldNum < numFields; ++fieldNum)
        {
            rho(fieldNum, n) = 0.5;
        }
    }

    for (const auto &nodeNo : bulkNodes)
        rhoOld(0, nodeNo) = rho(0, nodeNo);

    VectorField<LT> jVecOut(numFields, grid.size());

    // *********
    // LB FIELDS
    // *********
    LbField<LT> f(numFields, grid.size());
    LbField<LT> fTmp(numFields, grid.size());
    // initiate lb distributions
    for (auto nodeNo : bulkNodes)
    {
        std::valarray<lbBase_t> wVec(LT::w, LT::nQ);
        for (int fieldNum = 0; fieldNum < numFields; ++fieldNum)
        {
            f.set(fieldNum, nodeNo) = wVec * rho(fieldNum, nodeNo);
        }
    }


    ScalarField tagsField(1, grid.size());
    ScalarField nodeTypeField(1, grid.size());

    ScalarField applyBnd(1, grid.size());

    for (auto nodeNo: bulkNodes) {
      tagsField(0, nodeNo) = nodes.getTag(nodeNo);
      nodeTypeField(0, nodeNo) = nodes.getType(nodeNo);
    }

    // **********
    // OUTPUT VTK
    // **********

    Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs);
    output.add_file("lb_run_laplace");

    output.add_scalar_variables({"rho", "sd", "tags", "nodeType", "applyBnd"}, {rho, sd, tagsField, nodeTypeField, applyBnd});
    output.add_vector_variables({"PressureForceField"}, {jVecOut});
    // output.write(0);

    // *********
    // MAIN LOOP
    // *********
    // ----------------------------------------------------------------------------------- Total number of nodes
    lbBase_t sizeLocal = bulkNodes.size();
    lbBase_t sizeGlobal;
    MPI_Allreduce(&sizeLocal, &sizeGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i <= nIterations; i++)
    {
        const lbBase_t omega = 1.0 / tau;
        const std::valarray<lbBase_t> wVec(LT::w, LT::nQ);
        for (int fieldNum = 0; fieldNum < numFields; ++fieldNum)
        {
            for (auto nodeNo : bulkNodes)
            {
                const std::valarray<lbBase_t> fNode = f(fieldNum, nodeNo);
                const lbBase_t rhoNode = LT::qSum(fNode);
                rho(fieldNum, nodeNo) = rhoNode;
                const std::valarray<lbBase_t> fEq = wVec * rhoNode;
                const std::valarray<lbBase_t> fHat = fNode - omega * (fNode - fEq);
                for (int q = 0; q < LT::nQ; ++q)
                    fTmp(fieldNum, q, grid.neighbor(q, nodeNo)) = fHat[q];
            } // End nodes
        }

        // Swap data_ from fTmp to f;
        f.swapData(fTmp); // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************

        // Mpi send-recive
        mpiBoundary.communicateLbField(f, grid);
        // Bondary condition
        for (int fieldNum=0; fieldNum < numFields; ++fieldNum) {
	  bnd.apply(fieldNum, f, tau, nodes, grid, applyBnd);
        }

        // *************
        // WRITE TO FILE
        // *************
        if (((i % nItrWrite) == 0))
        {
            output.write(i);
            if (myRank == 0)
            {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
            }
            // Local changes
            lbBase_t deltaRho = 0;

            // Check mass conservation
            lbBase_t rhoTot = 0;
            for (auto nodeNo : bulkNodes)
            {
                const auto rhoNode = rho(0, nodeNo);
                rhoTot += rhoNode;
                deltaRho += std::abs(rhoNode - rhoOld(0, nodeNo));
                rhoOld(0, nodeNo) = rhoNode;
            }
            std::cout << "total rho = " << rhoTot << std::endl;
            lbBase_t deltaRhoGlobal;
            MPI_Allreduce(&deltaRho, &deltaRhoGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (myRank == 0) {
                std::cout << "Relative change = " << deltaRhoGlobal/sizeGlobal/nItrWrite << std::endl;
            }

            //----------------------------------------------------------------------------------- Write forces and pressures to file
            ScalarField psiTot(1, grid.size());
            VectorField<LT> jTot(1, grid.size());

            for (auto &nodeNo : bulkNodes)
            {
                psiTot(0, nodeNo) = 1;
                jTot.set(0, nodeNo) = 0;
            }

            for (int fieldNum = 0; fieldNum < numFields; ++fieldNum)
            {
                VectorField<LT> jVec = calcLaplaceForcing(fieldNum, f, tau, bulkNodes);
		//if (fieldNum == 0) {
		for (auto & nodeNo: bulkNodes) {
		  jVecOut.set(fieldNum, nodeNo) = jVec(0, nodeNo);
		}
		  //}
                ScalarField psi(1, grid.size());
                for (auto &nodeNo : bulkNodes)
                {
                    psi(0, nodeNo) = rho(fieldNum, nodeNo);
                    psiTot(0, nodeNo) -= psi(0, nodeNo);
                    jTot.set(0, nodeNo) -= jVec(0, nodeNo);
                }
                std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(fieldNum);
                jVec.writeToFile(laplacePressureDir + filename);
                psi.writeToFile(laplacePressureDir + filename);
            }

            std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(numFields);
            jTot.writeToFile(laplacePressureDir + filename);
            psiTot.writeToFile(laplacePressureDir + filename);

            VectorField<LT> jRead(numFields + 1, grid.size());
            ScalarField psiRead(numFields + 1, grid.size());
            for (int fieldNum = 0; fieldNum < (numFields + 1); ++fieldNum)
            {
                VectorField<LT> jTmp(1, grid.size());
                ScalarField psiTmp(1, grid.size());
                std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(fieldNum);
                jTmp.readFromFile(laplacePressureDir + filename);
                psiTmp.readFromFile(laplacePressureDir + filename);
                for (auto &nodeNo : bulkNodes)
                {
                    jRead.set(fieldNum, nodeNo) = jTmp(0, nodeNo);
                    psiRead(fieldNum, nodeNo) = psiTmp(0, nodeNo);
                }
            }
        }
    } // End iterations

    //----------------------------------------------------------------------------------- Write forces and pressures to file
    ScalarField psiTot(1, grid.size());
    VectorField<LT> jTot(1, grid.size());

    for (auto &nodeNo : bulkNodes)
    {
        psiTot(0, nodeNo) = 1;
        jTot.set(0, nodeNo) = 0;
    }

    for (int fieldNum = 0; fieldNum < numFields; ++fieldNum)
    {
        VectorField<LT> jVec = calcLaplaceForcing(fieldNum, f, tau, bulkNodes);
        ScalarField psi(1, grid.size());
        for (auto &nodeNo : bulkNodes)
        {
            psi(0, nodeNo) = rho(fieldNum, nodeNo);
            psiTot(0, nodeNo) -= psi(0, nodeNo);
            jTot.set(0, nodeNo) -= jVec(0, nodeNo);
        }
        std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(fieldNum);
        jVec.writeToFile(laplacePressureDir + filename);
        psi.writeToFile(laplacePressureDir + filename);
    }
    std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(numFields);
    jTot.writeToFile(laplacePressureDir + filename);
    psiTot.writeToFile(laplacePressureDir + filename);

    VectorField<LT> jRead(numFields + 1, grid.size());
    ScalarField psiRead(numFields + 1, grid.size());
    for (int fieldNum = 0; fieldNum < (numFields + 1); ++fieldNum)
    {
        VectorField<LT> jTmp(1, grid.size());
        ScalarField psiTmp(1, grid.size());
        std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(fieldNum);
        jTmp.readFromFile(laplacePressureDir + filename);
        psiTmp.readFromFile(laplacePressureDir + filename);
        for (auto &nodeNo : bulkNodes)
        {
            jRead.set(fieldNum, nodeNo) = jTmp(0, nodeNo);
            psiRead(fieldNum, nodeNo) = psiTmp(0, nodeNo);
        }
    }

    /*
    Output<LT> outputForce(grid, bulkNodes, outputDir, myRank, nProcs);
    outputForce.add_file("forcing");
    outputForce.add_scalar_variables({"pressure"}, {psiRead});
    outputForce.add_vector_variables({"force"}, {jRead});
    outputForce.write(0);
    */
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
