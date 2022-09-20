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
    std::string chimpDir = "./";   
    //std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
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

    // ---------------------------------------------------------------------------------- Read signed distance
    // need to be defined before pressure.
    auto sd = readSignedDistance("signed_distance", vtklb, nodes, grid);

    // ---------------------------------------------------------------------------------- Add pressure boundary nodes
    int maxBoundaryIndicator = addPressureBoundary("pressure_boundary", vtklb, nodes, grid);
    
    std::cout <<"maxBoundaryIndicator = " << maxBoundaryIndicator << std::endl;

    
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

    //                                    Output directory number
    //------------------------------------------------------------------------------------- Output directory number
    std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));  
    std::string outputDir2 = outputDir + "/out" + dirNum;
    

    // ---------------------------------------------------------------------------------- Setup pressure boundary
    std::vector<int> bulkNodes = findBulkNodes<LT>(nodes);
    laplaceBoundary<LT> bnd(bulkNodes, sd, nodes, grid);

    
    
    // setup the pressure boundary information (normal vector, distance to boundary)
    
    vtklb.toAttribute("boundary_distance");
    auto numPNodes = vtklb.numSubsetEntries();

    std::cout << "Rank " << myRank << ": " << numPNodes << " pressure nodes." << std::endl;
    
    std::vector<int> pressureNodes(numPNodes);
    std::vector<lbBase_t> pressureNodesS(numPNodes);
    //std::vector<std::valarray<lbBase_t>> pressureNodesNorm;

    std::vector<lbBase_t> pressureNodeNormX(numPNodes);
    std::vector<lbBase_t> pressureNodeNormY(numPNodes);
    std::vector<lbBase_t> pressureNodeNormZ(numPNodes);
    
    std::valarray<lbBase_t> TmpArr(0.0, LT::nD);
    std::vector<std::valarray<lbBase_t>> pressureNodesNorm(numPNodes, TmpArr);
    

    for (int n=0; n<numPNodes; ++n) {
      const auto ent = vtklb.getSubsetAttribure<lbBase_t>();
      //std::cout << "S, Node number = " << ent.nodeNo << "     value = " << ent.val << std::endl;
      pressureNodes[n] = ent.nodeNo;
      pressureNodesS[n] = ent.val;
    }
    std::cout << "RANK " << myRank <<" READ FROM INPUT FILE: Pressure nodes' s distance" << std::endl;
    
    vtklb.toAttribute("boundary_normal_x");
    if (vtklb.numSubsetEntries() != numPNodes)
      std::cout << "ERROR: NUMBER OF boundary_normal_x ENTRIES DOES NOT MATCH NUMBER OF PRESSURE NODES!" << std::endl;
    for (int n=0; n<vtklb.numSubsetEntries(); ++n) {
      const auto ent = vtklb.getSubsetAttribure<lbBase_t>();
      //std::cout << "X, Node number = " << ent.nodeNo << "     value = " << ent.val << std::endl;
      //pressureNodeNormX[n] = ent.val;
      pressureNodesNorm[n][0] = ent.val;
    }
    std::cout << "RANK " << myRank <<" READ FROM INPUT FILE: Pressure nodes' normal x-component" << std::endl;
    
    if(LT::nD>=2){
      vtklb.toAttribute("boundary_normal_y");
      if (vtklb.numSubsetEntries() != numPNodes)
      std::cout << "ERROR: NUMBER OF boundary_normal_y ENTRIES DOES NOT MATCH NUMBER OF PRESSURE NODES!" << std::endl;
      for (int n=0; n<vtklb.numSubsetEntries(); ++n) {
	const auto ent = vtklb.getSubsetAttribure<lbBase_t>();
	//std::cout << "Y, Node number = " << ent.nodeNo << "     value = " << ent.val << std::endl;
	pressureNodeNormY[n] = ent.val;
	pressureNodesNorm[n][1] = ent.val;
      }
      std::cout << "RANK " << myRank <<" READ FROM INPUT FILE: Pressure nodes' normal y-component" << std::endl;
    }
    if(LT::nD==3){
      vtklb.toAttribute("boundary_normal_z");
      if (vtklb.numSubsetEntries() != numPNodes)
      std::cout << "ERROR: NUMBER OF boundary_normal_z ENTRIES DOES NOT MATCH NUMBER OF PRESSURE NODES!" << std::endl;
      for (int n=0; n<vtklb.numSubsetEntries(); ++n) {
	const auto ent = vtklb.getSubsetAttribure<lbBase_t>();
	//std::cout << "Z, Node number = " << ent.nodeNo << "     value = " << ent.val << std::endl;
	pressureNodeNormZ[n] = ent.val;
	pressureNodesNorm[n][2] = ent.val;
      }
      std::cout << "RANK " << myRank <<" READ FROM INPUT FILE: Pressure nodes' normal z-component" << std::endl;
    }
    
  
    
    /*
    std::vector<int> pressureNodes;
    std::vector<lbBase_t> pressureNodesS;
    std::vector<std::valarray<lbBase_t>> pressureNodesNorm;  

    for (const auto & n: bulkNodes) {
        // norms.set(0, n) = normalFromSignedDistance(n, sd, grid);
      if (nodes.getTag(n) == 1 
      //|| nodes.getTag(n) == 3
      ) {
            pressureNodes.push_back(n);
            //pressureNodesS.push_back(0.5);
	    pressureNodesS.push_back(0.0);
            pressureNodesNorm.push_back((std::initializer_list<lbBase_t>){1, 0, 0});
        }
        if (nodes.getTag(n) == 2) {
            pressureNodes.push_back(n);
            //pressureNodesS.push_back(0.5);
	    pressureNodesS.push_back(0.0);
            pressureNodesNorm.push_back((std::initializer_list<lbBase_t>){-1, 0, 0});
        }
    }
    */
    bnd.addPressureNodeData(pressureNodes, pressureNodesS, pressureNodesNorm);

    
    // Set bulk nodes
    const int numFields = maxBoundaryIndicator - 1;
    
    
    //std::vector<int> bulkNodes = diffusion.findBulkNodes(vtklb, nodes);
    ScalarField rho(numFields, grid.size());
    std::cout <<" test 2 rank: " << myRank << std::endl;
    ScalarField rhoOld(1, grid.size());
    std::cout <<" test 3 rank: " << myRank << std::endl;
    
    
    // Initiate density from file
    // vtklb.toAttribute("init_rho");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        // auto rhoVal = vtklb.getScalarAttribute<lbBase_t>();
        for (int fieldNum=0; fieldNum < numFields; ++fieldNum) {
	        rho(fieldNum, n) = 0.5;
        }
    }

    
    
    for (const auto & nodeNo: bulkNodes)
        rhoOld(0, nodeNo) = rho(0, nodeNo);

    VectorField<LT> jVecOut(1, grid.size());

    
    
    // *********
    // LB FIELDS
    // *********
    LbField<LT> f(numFields, grid.size()); 
    LbField<LT> fTmp(numFields, grid.size());
    // initiate lb distributions
    for (auto nodeNo: bulkNodes) {
        std::valarray<lbBase_t> wVec(LT::w, LT::nQ);
        for (int fieldNum=0; fieldNum < numFields; ++fieldNum) {
            f.set(fieldNum, nodeNo) = wVec*rho(fieldNum, nodeNo);
        }
    }


    // **********
    // OUTPUT VTK
    // **********
    Output<LT> output(grid, bulkNodes, outputDir2, myRank, nProcs);
    output.add_file("lb_run_laplace");
    output.add_scalar_variables({"rho", "sd"}, {rho, sd});
    output.add_vector_variables({"PressureForceField"}, {jVecOut});
    //output.write(0);

    // *********
    // MAIN LOOP
    // *********
    // ----------------------------------------------------------------------------------- Total number of nodes
    lbBase_t sizeLocal = bulkNodes.size();
    lbBase_t sizeGlobal;
    MPI_Allreduce(&sizeLocal, &sizeGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i <= nIterations; i++) {
        const lbBase_t omega = 1.0/tau;
        const std::valarray<lbBase_t> wVec(LT::w, LT::nQ);
        for (int fieldNum=0; fieldNum < numFields; ++fieldNum) {
            for (auto nodeNo: bulkNodes) {
                const std::valarray<lbBase_t> fNode = f(fieldNum, nodeNo);
                const lbBase_t rhoNode = LT::qSum(fNode);
                rho(fieldNum, nodeNo) = rhoNode;
                const std::valarray<lbBase_t> fEq = wVec*rhoNode;
                const std::valarray<lbBase_t> fHat =  fNode - omega*(fNode - fEq);
                for (int q = 0; q < LT::nQ; ++q)
                    fTmp(fieldNum, q, grid.neighbor(q, nodeNo)) = fHat[q];
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
        for (int fieldNum=0; fieldNum < numFields; ++fieldNum) {
            bnd.apply(fieldNum, f, tau, nodes, grid);
        }

        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0)  ) {
            output.write(i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
            }
            // Local changes
            lbBase_t deltaRho = 0;

            // Check mass conservation
            lbBase_t rhoTot = 0;
            for (auto nodeNo: bulkNodes) {
                const auto rhoNode = rho(0, nodeNo);
                rhoTot += rhoNode;
                deltaRho += std::abs(rhoNode - rhoOld(0, nodeNo));
                rhoOld(0, nodeNo) = rhoNode;
            }
            std::cout << "total rho = " << rhoTot << std::endl;
            lbBase_t deltaRhoGlobal;
            MPI_Allreduce(&deltaRho, &deltaRhoGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (myRank == 0) {
                std::cout << "Relative change = " << deltaRhoGlobal/sizeGlobal << std::endl;
            }
           

            //----------------------------------------------------------------------------------- Write forces and pressures to file
            ScalarField psiTot(1, grid.size());
            VectorField<LT> jTot(1, grid.size());

            for (auto &nodeNo: bulkNodes) {
                psiTot(0, nodeNo) = 1;
                jTot.set(0, nodeNo) = 0; 
            } 

            for (int fieldNum=0; fieldNum < numFields; ++fieldNum) {
                VectorField<LT> jVec = calcLaplaceForcing(fieldNum, f, tau, bulkNodes);
		if (fieldNum == 0) {
		  for (auto & nodeNo: bulkNodes) {
		    jVecOut.set(0, nodeNo) = jVec(0, nodeNo);
		  }
		}
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
        }
    } // End iterations


    /*
    //----------------------------------------------------------------------------------- Write forces and pressures to file
    ScalarField psiTot(1, grid.size());
    VectorField<LT> jTot(1, grid.size());

    for (auto &nodeNo: bulkNodes) {
        psiTot(0, nodeNo) = 1;
        jTot.set(0, nodeNo) = 0; 
    } 

    for (int fieldNum=0; fieldNum < numFields; ++fieldNum) {
        VectorField<LT> jVec = calcLaplaceForcing(fieldNum, f, tau, bulkNodes);
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
    */



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
