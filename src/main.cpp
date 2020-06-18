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

#include "Input.h"
#include "Output.h"

#include<algorithm> // std::max

// SET THE LATTICE TYPE
#define LT D3Q19


int main()
{
  /*
    std::cout << "Test c2q" << std::endl;
    for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
            for (int k = -1; k <= 1; ++k)
            {
                std::vector<int> v = {i, j, k};
                std::cout << "q = " <<  LT::c2q(v) << "  vec = " << v << std::endl;
            }
  */
    // SETUP MPI
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // TEST PRINTOUT
    std::cout << "Begin test Two phase new" << std::endl;

    // SETUP THE INPUT AND OUTPUT PATHS
    //std::string chimpDir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/";
    std::string chimpDir = "./";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

    std::cout << "before read files" << std::endl;
    // READ BOUNDARY FILES WITH MPI INFORMATION
    MpiFile<LT> rankFile(mpiDir + "rank_" + std::to_string(myRank) + "_rank.mpi");
    MpiFile<LT> localFile(mpiDir + "rank_" + std::to_string(myRank) + "_local_labels.mpi");
    MpiFile<LT> globalFile(mpiDir + "rank_" + std::to_string(myRank) + "_global_labels.mpi");

    std::cout << "Sizes = " << rankFile.size() << " " << localFile.size() << " " << globalFile.size() << std::endl;

    // SETUP GRID
    std::cout << "grid" << std::endl;
    auto grid  = Grid<LT>::makeObject(localFile, rankFile, myRank);
    std::cout << "Process "<<myRank<<": grid.size = " << grid.size() << std::endl;

    // SETUP NODE
    std::cout << "nodes" << std::endl;
    auto nodes = Nodes<LT>::makeObject(localFile, rankFile, myRank, grid);

    // SETUP MPI BOUNDARY
    std::cout << "mpi boundary" << std::endl;
    BndMpi<LT> mpiBoundary(myRank);
    mpiBoundary.setup(localFile, globalFile, rankFile, nodes,  grid);

    std::vector<int> fluidBndNodes = findFluidBndNodes(nodes);
    std::vector<int> fluidBndNodesSouth;
    std::vector<int> fluidBndNodesNorth;
    
    
    for (auto nodeNo: fluidBndNodes){ 
      if (grid.pos(nodeNo, 1)==2)  fluidBndNodesSouth.push_back(nodeNo);
      if (grid.pos(nodeNo, 1)==rankFile.dim_global(1)-2)  fluidBndNodesNorth.push_back(nodeNo);
    }

    for (auto nodeNo: fluidBndNodesSouth) {
        if (!nodes.isFluid(nodeNo)) {
            std::cout << "node = " << nodeNo << " is not fluid " << std::endl;
            exit(1);
        }
    }

 
    /*
    std::cout<<std::endl;
    for (int i=0; i< freeSlipFluidBndNodesNorth.size(); i++) {
      std::cout<<"freeSlipFluidBndNodesNorth["<<i<<"] = "<<freeSlipFluidBndNodesNorth[i]<<" "<<std::endl;
    }
    
    for (int nodeNo=0; nodeNo<grid.size(); nodeNo++){
      if(grid.pos(nodeNo, 2)==0)
	std::cout<<"NodeNo = "<<nodeNo<<std::endl;
    }
    */
    
    // SETUP BOUNCE BACK BOUNDARY (fluid boundary)
    std::cout << "bbBnd" << std::endl;
    //HalfWayBounceBack<LT> bbBnd(fluidBndNodes, nodes, grid); // = makeFluidBoundary<HalfWayBounceBack>(nodes, grid);
    HalfWayBounceBack<LT> bbBndSth(fluidBndNodesSouth, nodes, grid); 
    HalfWayBounceBack<LT> bbBndNrth(fluidBndNodesNorth, nodes, grid);

    
    
    // SETUP FREE SLIP BOUNDARY (fluid boundary)
    
    std::cout << "freeSlipBndSouth" << std::endl;
    std::vector<int> normVecY = {0, 1, 0};
    //FreeSlipCartesian<LT> frSlpBndSth(normVecY, fluidBndNodesSouth, nodes, grid);
    FreeFlowCartesian<LT> frFlwBndSth(fluidBndNodesSouth, nodes, grid);
    
    std::cout << "freeSlipBndNorth" << std::endl;
    std::vector<int> normVecMinY = {0, -1, 0};
    FreeSlipCartesian<LT> frSlpBndNrth(normVecMinY, fluidBndNodesNorth, nodes, grid);
    

    // SETUP SOLID BOUNDARY
    std::vector<int> solidBnd = findSolidBndNodes(nodes);

    // SETUP BULK NODES
    std::vector<int> bulkNodes = findBulkNodes(nodes);


    std::cout<<"c(0) = "<<LT::c(0)<<std::endl;
    std::valarray<int> c_va(LT::c(0).data(), LT::nD);
    std::cout<<"c_va(0) = "<<c_va[0]<<std::endl;
    
    //---------------------END OF SETUP OF BOUNDARY CONDITION AND MPI---------------------

    // READ INPUT FILE
    Input input(inputDir + "input.dat");

    // Vector source
    VectorField<LT> bodyForce(1, 1);
    bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);

    int nIterations = static_cast<int>( input["iterations"]["max"]);

    lbBase_t tau0 = input["fluid"]["tau"][0];
    lbBase_t tau1 = input["fluid"]["tau"][1];
    
    
    lbBase_t sigma = input["fluid"]["sigma"];
    lbBase_t beta = input["fluid"]["beta"];

    lbBase_t tauDiff0 = input["diff-field"]["tauDiff"][0];
    lbBase_t tauDiff1 = input["diff-field"]["tauDiff"][1];
    lbBase_t tauDiff2 = input["diff-field"]["tauDiff"][2];

    lbBase_t H = input["Partial-Misc"]["H"];
    lbBase_t kinConst = input["Partial-Misc"]["kinConst"];

    lbBase_t CSmagorinsky = input["LES"]["CSmag"];;

    std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));

    std::string outDir2 = outputDir+"out"+dirNum;

    //std::string outDir2 = <std::string>(input["out"]["directory"][0]); //virker ikke!

    // SET DERIVED VARAIBLES
    lbBase_t nu0Inv = 1.0 / (LT::c2 * (tau0 - 0.5));
    lbBase_t nu1Inv = 1.0 / (LT::c2 * (tau1 - 0.5));

    lbBase_t DwInv = 1.0 / (LT::c2 * (tauDiff0 - 0.5));
    lbBase_t DdwInv = 1.0 / (LT::c2 * (tauDiff1 - 0.5));
    lbBase_t DsInv = DwInv; //1.0 / (LT::c2 * (tauDiff0 - 0.5));
    lbBase_t DoilInv = 1.0 / (LT::c2 * (tauDiff2 - 0.5));
    lbBase_t Dw = (LT::c2 * (tauDiff0 - 0.5));
    lbBase_t Ddw = (LT::c2 * (tauDiff1 - 0.5));
    lbBase_t Ds = Dw;//(LT::c2 * (tauDiff0 - 0.5));
    lbBase_t Doil = (LT::c2 * (tauDiff2 - 0.5));


    // SETUP LB FIELDS
    LbField<LT> f(1, grid.size());  // LBfield
    LbField<LT> fTmp(1, grid.size());  // LBfield

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(1, grid.size()); // LBfield
    VectorField<LT> vel(1, grid.size()); // LBfield
    ScalarField eff_nu(1, grid.size()); // LBfield
    ScalarField qSrc(1, grid.size()); // LBfield
    VectorField<LT> forceTot(1, grid.size()); // LBfield
    /*
    ScalarField velXKernel(1, grid.size()); // LBfield
    ScalarField velYKernel(1, grid.size()); // LBfield
    ScalarField velZKernel(1, grid.size()); // LBfield
    */
    
    // FILL MACROSCOPIC FIELDS
    //   Fluid densities and velocity

    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 1.0;
        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
	vel(0,1,nodeNo) = -0.0;
	
	//vel(0, 1, nodeNo) = 0.03;
	int ymax = rankFile.dim_global(1)-2;
	int y = grid.pos(nodeNo, 1);
	
	//if(y <= 5 && y > 2)
	/*
	if(y == 5)  
	  qSrc(0, nodeNo) = -0.02;
	  else*/
	if(y == (ymax-2))
	  qSrc(0, nodeNo) = 0.02;
	
    }

    lbBase_t r0 = rankFile.dim_global(0)/4;
    lbBase_t theta0 = /*1.5*3.1425;*/0.5*3.1425;//0;
    lbBase_t Tperiod = 10000.;
    lbBase_t startTime = 4*Tperiod;
    lbBase_t angVel = 2*3.1415/Tperiod;
    lbBase_t R = 10;
    lbBase_t epsilon = 0.;//5.0;//5;//3;
    lbBase_t meanSphereRho=1.;
    lbBase_t meanSphereRhoGlobal=1.;

    //Fixed outlet pressure -----------------------------------------------
    lbBase_t pHat = 0.985*LT::c2;
    
    std::vector<int> rotCenterPos = {rankFile.dim_global(0)/2, (int)(rankFile.dim_global(1)-1.4*(r0+R+epsilon)), rankFile.dim_global(2)/2};
    //std::vector<int> centerPos = {rankFile.dim_global(0)/2, (int)(0.67*rankFile.dim_global(1)), 0};  

    
    //---------------------END OF INPUT TO INITIALIZATION OF FIELDS---------------------


    // INITIATE LB FIELDS
    // -- fluid field
    initiateLbField(0, 0, 0, bulkNodes, rho, vel, f);  // LBinitiatefield

    // JLV
    //---------------------SETUP OUTPUT---------------------
    std::vector<std::vector<int>> node_pos;
    node_pos.reserve(bulkNodes.size());
    for (const auto& node:bulkNodes) {
        node_pos.push_back(grid.pos(node));
    }

    Output output(globalFile.dim_global(), outDir2, myRank, nProcs-1, node_pos);
    output.add_file("fluid");
    output["fluid"].add_variable("rho", rho.get_data(), rho.get_field_index(0, bulkNodes), 1);
    output["fluid"].add_variable("vel", vel.get_data(), vel.get_field_index(0, bulkNodes), LT::nD);
    output["fluid"].add_variable("eff_nu", eff_nu.get_data(), eff_nu.get_field_index(0, bulkNodes), 1);
    output["fluid"].add_variable("qSrc", qSrc.get_data(), qSrc.get_field_index(0, bulkNodes), 1);
    output["fluid"].add_variable("force", forceTot.get_data(), forceTot.get_field_index(0, bulkNodes), LT::nD);
    
    output.write("fluid", 0);

    //---------------------END OF SETUP OUTPUT---------------------


    // -----------------MAIN LOOP------------------
    /* Comments to main loop:
     * Calculation of cu is kept outside of calcOmega and calcDeltaOmega
     *  to avoid recalculation  of the dot-product.
     * It is therefore natural to give all scalar products as inputs to the
     *  calcOmega and calcDeltaOmega.
     * Note that we have a oneIteration function in 'LBiteration.h'. This seemed
     *   to affect the speed an gave a slow down of 10-15 %.
     */



    for (int i = 0; i <= nIterations; i++) {

      /*
      for (auto nodeNo : bulkNodes) {
	// UPDATE grad vel kernel
	velXKernel(0,nodeNo) = vel(0,0,nodeNo);
	velYKernel(0,nodeNo) = vel(0,1,nodeNo);
	velZKernel(0,nodeNo) = vel(0,2,nodeNo);	
      }
      
      for (auto nodeNo : solidBnd){
	// UPDATE grad vel kernel
	velXKernel(0,nodeNo) = 0.0;
	velYKernel(0,nodeNo) = 0.0;
	velZKernel(0,nodeNo) = 0.0;
      }

      mpiBoundary.communciateScalarField(velXKernel);
      mpiBoundary.communciateScalarField(velYKernel);
      mpiBoundary.communciateScalarField(velZKernel);
      */

      //----------------------------------Rotating sphere---------------------------------------
      lbBase_t angVelLoc = angVel;
      int phaseTLoc = i;

      phaseTLoc = i - (int)startTime;
      if(i<startTime){
	angVelLoc = 0;
	phaseTLoc = 0;
      }
	        
      std::valarray<lbBase_t> sphereCenterLoc = sphereCenter(rotCenterPos, r0, theta0, angVelLoc, phaseTLoc);
      //End-------------------------------Rotating sphere---------------------------------------
         
        for (auto nodeNo : bulkNodes) {
            // UPDATE MACROSCOPIC DENSITIES
            lbBase_t rhoTotNode = rho(0, nodeNo) = calcRho<LT>(f(0,nodeNo));  // LBmacroscopic
        }  // End for all bulk nodes



	//----------------------------------Rotating sphere---------------------------------------
	int numSphereNodes = 0;
	int numSphereNodesGlobal;
	meanSphereRho = 0.0;
	for (auto nodeNo : bulkNodes) {
	  std::valarray<lbBase_t> fTot = f(0, nodeNo);
	  std::vector<int> pos = grid.pos(nodeNo);

	  lbBase_t fromSphereCenterSq = 0;
	  for(int i = 0; i < LT::nD; ++i){
	    fromSphereCenterSq += (pos[i]-sphereCenterLoc[i])*(pos[i]-sphereCenterLoc[i]);
	  }
	  
	  	  
	  if (fromSphereCenterSq <= R*R){
	    meanSphereRho += rho(0, nodeNo);
	    numSphereNodes++;
	  }
        }  // End for all bulk nodes
	MPI_Allreduce(&meanSphereRho, &meanSphereRhoGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&numSphereNodes, &numSphereNodesGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	meanSphereRhoGlobal/=numSphereNodesGlobal;
	//End-------------------------------Rotating sphere---------------------------------------
	
	
        for (auto nodeNo: bulkNodes) {

            // Set the local total lb distribution
            std::valarray<lbBase_t> fTot = f(0, nodeNo);
	    
	    std::vector<int> pos = grid.pos(nodeNo);
	    std::valarray<lbBase_t> forceNode = bodyForce(0, 0);
	    
            // -- total density
            lbBase_t rhoTotNode = rho(0, nodeNo);
	    lbBase_t qSrcNode = qSrc(0, nodeNo);

	    //Fixed outlet pressure -----------------------------------------------
	    
	    //if(pos[1] < 4 && pos[1]>1){
	    
	    if(pos[1]<=5 && qSrcNode!=0){
	      qSrcNode = 2*(pHat*LT::c2Inv - LT::qSum(fTot));
	      qSrc(0, nodeNo)=qSrcNode;
	    }

	    /*
	    std::valarray<lbBase_t> gradVelX = grad(velXKernel, 0, nodeNo, grid);
	    std::valarray<lbBase_t> gradVelY = grad(velYKernel, 0, nodeNo, grid);
	    std::valarray<lbBase_t> gradVelZ = grad(velZKernel, 0, nodeNo, grid);
	    
	    
	    std::valarray<lbBase_t> gradVelTensor(LT::nD*LT::nD);
	    int it=0;
	    for(int i = 0; i < LT::nD; ++i){
	      gradVelTensor[it]=gradVelX[i];
	      it++;
	      gradVelTensor[it]=gradVelY[i];
	      it++;
	      gradVelTensor[it]=gradVelZ[i];
	      it++;
	    }
	    std::valarray<lbBase_t> gradVelTensorT(LT::nD*LT::nD);
	    it=0;
	    for(int i = 0; i < LT::nD; ++i){
	      gradVelTensorT[it]=gradVelX[i];
	      it++;
	    }
	    for(int i = 0; i < LT::nD; ++i){
	      gradVelTensorT[it]=gradVelY[i];
	      it++;
	    }
	    for(int i = 0; i < LT::nD; ++i){
	      gradVelTensorT[it]=gradVelZ[i];
	      it++;
	    }

	    std::valarray<lbBase_t> WALETensor(LT::nD*LT::nD);
	    std::valarray<lbBase_t> g2 = LT::matrixMultiplication(gradVelTensor, gradVelTensor);
	    std::valarray<lbBase_t> gT2 = LT::matrixMultiplication(gradVelTensorT, gradVelTensorT);

	    
	    WALETensor = 0.5*(g2 + gT2);
	    lbBase_t trace1 = LT::traceOfMatrix(gradVelTensor)/LT::nD;
	    WALETensor[0] += trace1;
	    WALETensor[4] += trace1;
	    WALETensor[8] += trace1;

	    //lbBase_t WALETensorAbs = sqrt(2*LT::contractionRank2(WALETensor, WALETensor));
	    lbBase_t WALETensorAbs = LT::contractionRank2(WALETensor, WALETensor);
	    //std::cout<<WALETensorAbs<<std::endl;
	    */
	    
	    //----------------------------------Rotating sphere---------------------------------------
	    std::valarray<lbBase_t> rotPointForceNode = rotatingPointForce<LT>(fTot, pos, rotCenterPos, sphereCenterLoc, r0, theta0, angVelLoc, R, epsilon, phaseTLoc);
	    qSrcNode+=qSrcConstSphere<LT>(fTot, pos, sphereCenterLoc, R, meanSphereRhoGlobal, epsilon);
	    //End-------------------------------Rotating sphere---------------------------------------
	    // -- force
            
	    forceNode += rotPointForceNode;

	    /*
	    if(pos[1] == 1 && pos[0]>0 && pos[0]< rankFile.dim_global(0)-1){
	      int qYDir = 1;
	      int neighNodeNo = grid.neighbor(qYDir,nodeNo);
	      //std::cout<<neighNodeNo<<std::endl;
	      //forceNode[0] = 2*(vel(0,0,neighNodeNo)*LT::qSum(fTot)- LT::qSumC(fTot)[0]);
	      forceNode[0] = 2*(- LT::qSumC(fTot)[0]);
	      //fTot = f(0, neighNodeNo);
	    }
	    */
	    
	    forceTot.set(0,nodeNo) = forceNode;
	    
            // -- velocity
            std::valarray<lbBase_t> velNode = calcVel<LT>(fTot, rhoTotNode, forceNode);  // LBmacroscopic
            vel.set(0, nodeNo) = velNode;
	    
	    
	    //qSrc(0, nodeNo)=qSrcNode;
	    
	   
	    rhoTotNode = rho(0, nodeNo)+=0.5*qSrcNode;

	     
	    //std::valarray<lbBase_t> Stilde = calcShearRateTilde<LT>(fTot, rhoTotNode, velNode, forceNode, qSrcNode); // LBmacroscopic

	    lbBase_t tau = tau0;

	   	    
	    

	    // CALCULATE SMAGORINSKY CONSTANT AND LES EFFECTIVE TAU

	    std::valarray<lbBase_t> fTotSym(LT::nQ);
	    for (int q = 0; q < LT::nQ; ++q){
	      fTotSym[q] = 0.5*(fTot[q]+fTot[LT::reverseDirection(q)]);
	    } 
	    //std::valarray<lbBase_t> StildeLowTri = calcShearRateTildeLowTri<LT>(fTot, rhoTotNode, velNode, forceNode, qSrcNode); // LBmacroscopic
	    std::valarray<lbBase_t> StildeLowTri = calcShearRateTildeLowTri<LT>(fTotSym, rhoTotNode, velNode, forceNode, qSrcNode); // LBmacroscopic
	    
	    lbBase_t StildeAbs= sqrt(2*LT::contractionLowTri(StildeLowTri, StildeLowTri));

	    std::valarray<lbBase_t> sphereShellUnitNormalNode = sphereShellUnitNormal<LT>(pos, sphereCenterLoc);
	    std::valarray<lbBase_t> tangentialVelNormNode = tangentialUnitVector<LT>(velNode, sphereShellUnitNormalNode);

	    //std::valarray<lbBase_t> ELowTri = 1/(rhoTotNode*tau0)*LT::c2Inv*calcStrainRateTildeLowTri<LT>(fTotSym, rhoTotNode, velNode, forceNode, qSrcNode);
	    std::valarray<lbBase_t> ELowTri = 1/(rhoTotNode*tau0)*LT::c2Inv*StildeLowTri;

	    lbBase_t  wallVelGrad = LT::dot(LT::contractionLowTriVec(ELowTri,sphereShellUnitNormalNode),tangentialVelNormNode);
	    wallVelGrad=sqrt(wallVelGrad*wallVelGrad);
	    
	    lbBase_t friction_vel = sqrt(LT::c2*(tau0-0.5)*wallVelGrad);
	    
	   
            lbBase_t CSmagorinskyEff= CSmagorinsky;

	    //----------------------------------Rotating sphere---------------------------------------
	    lbBase_t fromSphereCenterSq = 0;
	    for(int i = 0; i < LT::nD; ++i){
	      fromSphereCenterSq += (pos[i]-sphereCenterLoc[i])*(pos[i]-sphereCenterLoc[i]);
	    }
	    	    
	    if (fromSphereCenterSq < R*R)
	      CSmagorinskyEff = 0.0;
	    else{
	      
	      lbBase_t yPlus = (sqrt(fromSphereCenterSq)-R)*friction_vel/(LT::c2*(tau0-0.5));
	      //lbBase_t yPlus = sqrt(fromSphereCenterSq)-R;
	      lbBase_t APlus = 25;
	      CSmagorinskyEff = CSmagorinsky*(1-exp(-yPlus/APlus));
	      
	      //CSmagorinskyEff = CSmagorinsky;//0.5;
	    }
	    //End-------------------------------Rotating sphere---------------------------------------

	    tau = 0.5*(tau0+sqrt(tau0*tau0+2*CSmagorinskyEff*CSmagorinskyEff*LT::c4Inv*StildeAbs/rhoTotNode));

	    //tau = tau0+ CSmagorinskyEff*CSmagorinskyEff*LT::c2Inv*WALETensorAbs;
	    
	    lbBase_t tauAntiNode=tau;

	    
	    
	    //ONLY FOR OUTPUT PURPOSES
	    eff_nu(0, nodeNo) = LT::c2*(tau - 0.5);
	    
	    
	    // CALCULATE BGK COLLISION TERM
            lbBase_t uu = LT::dot(velNode, velNode);  // Square of the velocity
            std::valarray<lbBase_t> cu = LT::cDotAll(velNode);  // velocity dotted with lattice vectors
            //std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fTot, tau, rhoTotNode, uu, cu);  // LBcollision
	    std::valarray<lbBase_t> omegaBGK = calcOmegaBGKTRT<LT>(fTot, tau, tauAntiNode, rhoTotNode, uu, cu);  // LBcollision
  
            // CALCULATE FORCE CORRECTION TERM
            lbBase_t  uF = LT::dot(velNode, forceNode);
            std::valarray<lbBase_t>  cF = LT::cDotAll(forceNode);
            //std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);  // LBcollision
	    std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaFTRT<LT>(tau, tauAntiNode, cu, uF, cF);  // LBcollision
	    
	    // CALCULATE MASS SOURCE CORRECTION TERM
            //std::valarray<lbBase_t> deltaOmegaQ = calcDeltaOmegaQ<LT>(tau, cu, uu, qSrcNode);  // LBcollision
	    std::valarray<lbBase_t> deltaOmegaQ = calcDeltaOmegaQTRT<LT>(tau, tauAntiNode, cu, uu, qSrcNode);  // LBcollision

            // COLLISION AND PROPAGATION


            //-----------------------
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be

                fTmp(0, q,  grid.neighbor(q, nodeNo)) =    fTot[q]            + omegaBGK[q]           + deltaOmegaF[q] + deltaOmegaQ[q];
                //-----------------------
            }
        } // End nodes

        // PRINT

        if ( (i % static_cast<int>(input["iterations"]["write"])) == 0) {

            // JLV
            output.write("fluid", i);
            // JLV
            if (myRank==0)
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
        }

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield


        // MPI Boundary
        // Hente verdier fra ghost
        // Sette i bulk
        mpiBoundary.communicateLbField(0, f, grid);


        // BOUNDARY CONDITIONS
        //bbBnd.apply(0, f, grid);  // LBboundary
	//bbBndSth.apply(0, f, grid);  // LBboundary
	//bbBndNrth.apply(0, f, grid);  // LBboundary
	//frSlpBndSth.apply(0, f, grid);  // LBboundary
	frSlpBndNrth.apply(0, f, grid);  // LBboundary
	frFlwBndSth.apply(0, f, grid);  // LBboundary

	
    } // End iterations
    // -----------------END MAIN LOOP------------------

    //mpi.end();

    MPI_Finalize();

    return 0;
}
