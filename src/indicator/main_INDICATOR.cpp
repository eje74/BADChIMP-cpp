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
// TO RUN PROGRAM: type "mpirun -np <#procs> bdchmp" in command
// line in main directory
//
// //////////////////////////////////////////////

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "../LBSOLVER.h"
#include "../IO.h"

#include<algorithm> // std::max

// SET THE LATTICE TYPE
//#define LT D2Q9
#define LT D3Q19
template <typename T>
using Output = LBOutputUnstructured<LT, T, VTK::BINARY, VTK::voxel>;
int main()
{
    // SETUP MPI
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

    // SETUP BOUNCE BACK BOUNDARY (fluid boundary)
    std::cout << "bbBnd" << std::endl;
    HalfWayBounceBack<LT> bbBnd(findBulkNodes(nodes), nodes, grid); // = makeFluidBoundary<HalfWayBounceBack>(nodes, grid);

    // SETUP SOLID BOUNDARY
    std::vector<int> solidBnd = findSolidBndNodes(nodes);

    // SETUP BULK NODES
    std::vector<int> bulkNodes = findBulkNodes(nodes);

    // ****************
    // SETUP MASS SOURCES
    // ****************
    // Source markers:
    
    std::vector<int> sourceMarker(grid.size());
    //   read from file
    vtklb.toAttribute("source");
    std::vector<int> constDensNodes, sourceNodes;  
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        int val = vtklb.getScalarAttribute<int>();
        sourceMarker[nodeNo] = val;
        if (val == 1) { // sink/source
	  		sourceNodes.push_back(nodeNo);
        } else if (val == 2) { // Const pressure
	  		constDensNodes.push_back(nodeNo);
		}
    }
    
    // Vector source
    VectorField<LT> bodyForce(1, 1);
    bodyForce.set(0, 0) = input["fluid"]["bodyforce"];

    int nIterations = input["iterations"]["max"];

    lbBase_t tau0 = input["fluid"]["tau"][0];
    lbBase_t tau1 = input["fluid"]["tau"][1];
    lbBase_t sigma = input["fluid"]["sigma"];
    lbBase_t beta = input["fluid"]["beta"];

    lbBase_t tauDiff0 = input["diff-field"]["tauDiff"][0];
    lbBase_t tauDiff1 = input["diff-field"]["tauDiff"][1];
    lbBase_t tauDiff2 = input["diff-field"]["tauDiff"][2];

    lbBase_t H = input["Partial-Misc"]["H"];
    lbBase_t kinConst = input["Partial-Misc"]["kinConst"];
    
	std::string dirNum = input["out"]["directoryNum"];
    std::string outDir2 = "output/out" + dirNum;
    
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

    lbBase_t testX = 0;
    for (int q = 0; q < LT::nQNonZero_; ++q) {
      testX += LT::w[q] * LT::c(q,0)*LT::c(q,0) /  LT::cNorm[q];
    }
    std::cout<<"Recolor Constant = "<<testX<<std::endl;
    lbBase_t betaPrime = testX*LT::c2Inv*beta;
    

    // Scalar source
    ScalarField Q(2, grid.size());
    for (int n = 0; n < Q.size(); ++n) {
        Q(0, n) = 0.0;
        Q(1, n) = 0.0;
    }

    
    
    // SETUP LB FIELDS
    LbField<LT> f(1, grid.size());  // LBfield
    LbField<LT> fTmp(1, grid.size());  // LBfield
    LbField<LT> gInd(1, grid.size()); // LB indicator field
    LbField<LT> gIndTmp(1, grid.size()); // LB indicator field

    LbField<LT> g(3, grid.size()); // LBfield diffusion
    LbField<LT> gTmp(3, grid.size()); // LBfield diffusion
   
    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(1, grid.size()); // LBfield
    VectorField<LT> vel(1, grid.size()); // LBfield
    ScalarField cField(1, grid.size()); // LBfield
    
    ScalarField cgNormField(1, grid.size()); // LBfield
    ScalarField cgFieldX(1, grid.size()); // LBfield
    ScalarField cgFieldY(1, grid.size()); // LBfield
    ScalarField cgFieldZ(1, grid.size()); // LBfield
    ScalarField cgFieldNormX(1, grid.size()); // LBfield
    ScalarField cgFieldNormY(1, grid.size()); // LBfield
    ScalarField cgFieldNormZ(1, grid.size()); // LBfield
    ScalarField kappaField(2, grid.size()); // LBfield
    ScalarField divGradColorField(1, grid.size()); // LBfield
    ScalarField divGradPhiBField(1, grid.size()); // LBfield
   
    ScalarField diffWField(1, grid.size()); // LBfield

    
    ScalarField R(1, grid.size()); // LBfield
    ScalarField waterChPot(1, grid.size()); // LBfield
    ScalarField indField(1, grid.size()); // LB indicator field
    VectorField<LT> gradients(5, grid.size()); // LBfield
    VectorField<LT> divGradF(1, grid.size()); // LBfield
    VectorField<LT> gradAbsFField(1, grid.size()); // LBfield
    //TEST
    VectorField<LT> diffgradients(4, grid.size()); // LBfield

    ScalarField rhoDiff(3, grid.size()); // LBfield

    //TEST
    ScalarField phiDiff0(1, grid.size()); // LBfield
    ScalarField phiDiff1(1, grid.size()); // LBfield
    ScalarField phiDiff2(1, grid.size()); // LBfield
    ScalarField phiIndicator(1, grid.size()); // LBfield
    
    // FILL MACROSCOPIC FIELDS
    //   Fluid densities and velocity
    std::srand(8549388);

    // ****************
    // SETUP MASS DENSITIES
    // ****************
    // Density markers:
    // read from file
    vtklb.toAttribute("rho0");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        float val = vtklb.getScalarAttribute<float>();
        rho(0, nodeNo) = val;
	for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }
    vtklb.toAttribute("phiInd");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        float val = vtklb.getScalarAttribute<float>();
        indField(0, nodeNo) = val*rho(0, nodeNo);
	rhoDiff(0, nodeNo) = 0.0;
	rhoDiff(1, nodeNo) = 0.0;
	rhoDiff(2, nodeNo) = 0.0;
    }
    /*
    vtklb.toAttribute("phiDiff0");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        float val = vtklb.getScalarAttribute<float>();
        rhoDiff(0, nodeNo) = val*rho(0, nodeNo);	
    }
    vtklb.toAttribute("phiDiff1");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        float val = vtklb.getScalarAttribute<float>();
        rhoDiff(1, nodeNo) = val*rho(0, nodeNo);
    }
    vtklb.toAttribute("phiDiff2");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        float val = vtklb.getScalarAttribute<float>();
        rhoDiff(2, nodeNo) = val*rho(0, nodeNo);
	//rhoDiff(2, nodeNo) = rho(0, nodeNo) - indField(0, nodeNo) - rhoDiff(0, nodeNo) - rhoDiff(1, nodeNo);
    }
    */
    
    auto global_dimensions = vtklb.getGlobaDimensions();

    //   Solid boundary (Wettability)
      
    vtklb.toAttribute("wettability");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
      float val = vtklb.getScalarAttribute<float>();
      if(nodes.isSolidBoundary(nodeNo)){
	indField(0, nodeNo) = val;
	rhoDiff(0, nodeNo) = 0.0; // LB Diff field
	rhoDiff(1, nodeNo) = 0.0; // LB Diff field
	rhoDiff(2, nodeNo) = 1-val; // LB Diff field
      }
    }

    int numNodes = bulkNodes.size();
    int numNodesGlobal;
    MPI_Allreduce(&numNodes, &numNodesGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    
    //---------------------END OF INPUT TO INITIALIZATION OF FIELDS---------------------

    
    // INITIATE LB FIELDS
    // -- total fluid
    initiateLbField(0, 0, 0, bulkNodes, rho, vel, f);  // LBinitiatefield
    // -- indicator field
    initiateLbField(0, 0, 0, bulkNodes, indField, vel, gInd);  // LBinitiatefield
    // -- diffusive fields
    initiateLbField(0, 0, 0, bulkNodes, rhoDiff, vel, g);  // LBinitiatefield
    initiateLbField(1, 1, 0, bulkNodes, rhoDiff, vel, g);  // LBinitiatefield
    initiateLbField(2, 2, 0, bulkNodes, rhoDiff, vel, g);  // LBinitiatefield
    
    
    // **********
    // OUTPUT VTK
    // **********
    // Setup output file
    Output<double> output(grid, bulkNodes, outDir2, myRank, nProcs);
    output.add_file("fluid");
    output.add_variables({ {{"rho"},rho}, {{"vel"},vel}, {{"water_fluid"},indField}, 
	                       {{"salt_in_water","water_in_oil","salt_in_oil"},rhoDiff},
					       {{"waterGrad","colorGrad","flWGrad","dWGrad"}, gradients}, 
						   {{"kappa"}, kappaField}, {{"divgradPhiB"}, divGradPhiBField}, 
						   {{"divgradCol"}, divGradColorField}, {{"WchP"}, waterChPot},
						   {{"R_dw"}, R}, {{"gradAbsF"},gradAbsFField} });
	// output["fluid"].add_variable("rho", rho.get_data(), rho.get_field_index(0, bulkNodes), 1);
    // output["fluid"].add_variable("vel", vel.get_data(), vel.get_field_index(0, bulkNodes), LT::nD);
    // output["fluid"].add_variable("water_fluid", indField.get_data(), indField.get_field_index(0, bulkNodes), 1);
    // output["fluid"].add_variable("salt_in_water", rhoDiff.get_data(), rhoDiff.get_field_index(0, bulkNodes), 1);
    // output["fluid"].add_variable("water_in_oil", rhoDiff.get_data(), rhoDiff.get_field_index(1, bulkNodes), 1);
    // output["fluid"].add_variable("salt_in_oil", rhoDiff.get_data(), rhoDiff.get_field_index(2, bulkNodes), 1);
    // output["fluid"].add_variable("waterGrad", gradients.get_data(), gradients.get_field_index(0, bulkNodes), LT::nD);
    // output["fluid"].add_variable("colorGrad", gradients.get_data(), gradients.get_field_index(1, bulkNodes), LT::nD);
    // output["fluid"].add_variable("flWGrad", gradients.get_data(), gradients.get_field_index(2, bulkNodes), LT::nD);
    // output["fluid"].add_variable("dWGrad", gradients.get_data(), gradients.get_field_index(3, bulkNodes), LT::nD);
    // output["fluid"].add_variable("kappa", kappaField.get_data(), kappaField.get_field_index(0, bulkNodes), 1);
    // output["fluid"].add_variable("divgradPhiB", divGradPhiBField.get_data(), divGradPhiBField.get_field_index(0, bulkNodes), 1);
    // output["fluid"].add_variable("divgradCol", divGradColorField.get_data(), divGradColorField.get_field_index(0, bulkNodes), 1);
    // output["fluid"].add_variable("WchP", waterChPot.get_data(), waterChPot.get_field_index(0, bulkNodes), 1);
    // output["fluid"].add_variable("R_dw", R.get_data(), R.get_field_index(0, bulkNodes), 1);
    // output["fluid"].add_variable("gradAbsF", gradAbsFField.get_data(), gradAbsFField.get_field_index(0, bulkNodes), 1);
         
    output.write(0);
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

        for (auto nodeNo : bulkNodes) {
            // UPDATE MACROSCOPIC DENSITIES
            // Calculate rho for each phase
	    lbBase_t rhoTotNode = rho(0, nodeNo) = calcRho<LT>(f(0,nodeNo));  // LBmacroscopic
	    lbBase_t indicator0Node = indField(0, nodeNo) = calcRho<LT>(gInd(0,nodeNo));  // LBmacroscopic
	    lbBase_t rhoDiff0Node = rhoDiff(0, nodeNo) = calcRho<LT>(g(0,nodeNo));  // LBmacroscopic
	    lbBase_t rhoDiff1Node = rhoDiff(1, nodeNo) = calcRho<LT>(g(1,nodeNo));  // LBmacroscopic
	    lbBase_t rhoDiff2Node = rhoDiff(2, nodeNo) = calcRho<LT>(g(2,nodeNo));  // LBmacroscopic

	    	    
            // Calculate color gradient kernel
	    const lbBase_t rhoType0Node = indicator0Node+rhoDiff0Node;
	    const lbBase_t rhoType1Node = rhoTotNode-rhoType0Node-rhoDiff1Node;
	    const lbBase_t cType0Node = rhoType0Node/(rhoTotNode);
	    const lbBase_t cType1Node = 1-cType0Node;
	    const lbBase_t rhoColorNode = rhoType0Node + rhoType1Node;
	    
	    //cField(0, nodeNo) = (rhoType0Node-rhoType1Node)/rhoTotNode;
	    cField(0, nodeNo) = (rhoType0Node-rhoType1Node)/rhoColorNode;
	    diffWField(0, nodeNo) = rhoDiff1Node;
	    phiDiff0(0, nodeNo) = rhoDiff0Node/rhoTotNode;
	    phiDiff1(0, nodeNo) = rhoDiff1Node/rhoTotNode;
	    phiDiff2(0, nodeNo) = rhoDiff2Node/rhoTotNode;
	    phiIndicator(0,nodeNo) = indicator0Node/rhoTotNode;
	    
	    // Calculate water gradient kernel
	    //waterChPot(0, nodeNo) = (indicator0Node+rhoDiff1Node/H);
	    //waterChPot(0, nodeNo) = indicator0Node/rhoTotNode+rhoDiff1Node/H/rhoTotNode;
	    //waterChPot(0, nodeNo) = (indicator0Node+rhoDiff1Node/H);
	    //waterChPot(0, nodeNo) = (indicator0Node-rhoDiff1Node/H);//*(indicator0Node-rhoDiff1Node/H);

	    //waterChPot(0, nodeNo) = (indicator0Node*indicator0Node-rhoDiff1Node*rhoDiff1Node/H);

	    //waterChPot(0, nodeNo) = (indicator0Node-rhoDiff1Node/H)*(indicator0Node-rhoDiff1Node/H)*(indicator0Node-rhoDiff1Node/H)-(indicator0Node-rhoDiff1Node/H);
	    waterChPot(0, nodeNo) = indicator0Node+rhoDiff1Node/H/rhoTotNode ;

	    

	    
        }  // End for all bulk nodes

        // Need to set the density when when calculating the densities
        // Predefine two regions:
        //   1) Pressure is constant
        //   2) Constant sink
        // Strategy - When mass is removed, then remove both oil and water
        //          - If mass is added at the top then then only add oil.

	/*Phase sources Q at the Constant-Density (Pressure) nodes are calculated*/
	
	for (auto nodeNo: constDensNodes) {
            lbBase_t rho0Node = rho(0, nodeNo);
	    lbBase_t indicator0Node = indField(0, nodeNo);
            //setConstDensity(Q(0, nodeNo), Q(1, nodeNo), rho0Node, rho1Node, 1.0);
	    lbBase_t rhoConst = 1.0;
	    Q(0, nodeNo) = 2*(rhoConst - rho0Node);
	    R(0, nodeNo) = indicator0Node/rho0Node*Q(0, nodeNo);
	    
	    rho0Node = rho(0, nodeNo)+=0.5*Q(0, nodeNo);
	    indicator0Node = indField(0, nodeNo)+= 0.5*R(0, nodeNo);
            // Calculate color gradient kernel
            cField(0, nodeNo) = 0;//(rho0Node - rho1Node)/(rho0Node + rho1Node);
        }
	
	/*Phase sources Q at the Constant-rate nodes are calculated*/
	
	for (auto nodeNo: sourceNodes) {
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t indicator0Node = indField(0, nodeNo);
	    //setConstSource(Q(0, nodeNo), R(0, nodeNo), rho0Node, rho1Node, 1e-4);
	    
	    Q(0, nodeNo) = 1e-3;
	    R(0, nodeNo) = 1e-3;
	    rho0Node = rho(0, nodeNo)+=0.5*Q(0, nodeNo);
	    indicator0Node = indField(0, nodeNo)+= 0.5*R(0, nodeNo);
            // Calculate color gradient kernel
            cField(0, nodeNo) = 0;//(rho0Node - rho1Node)/(rho0Node + rho1Node);
        }
	

        for (auto nodeNo: solidBnd) {
	  const lbBase_t rhoTotNode = rho(0, nodeNo);
	  const lbBase_t indicator0Node = indField(0, nodeNo);
	  const lbBase_t rhoDiff0Node = rhoDiff(0, nodeNo);
	  const lbBase_t rhoDiff1Node = rhoDiff(1, nodeNo);
	  const lbBase_t rhoDiff2Node = rhoDiff(2, nodeNo);

	  const lbBase_t rhoType0Node = indicator0Node+rhoDiff0Node;
	  const lbBase_t rhoType1Node = rhoTotNode-rhoType0Node;
	  
	  cField(0, nodeNo) = 2*indicator0Node-1;

	  diffWField(0, nodeNo) = rhoDiff1Node;
	 
	  phiDiff0(0, nodeNo) = rhoDiff0Node/rhoTotNode;
	  phiDiff1(0, nodeNo) = rhoDiff1Node/rhoTotNode;
	  phiDiff2(0, nodeNo) = rhoDiff2Node/rhoTotNode;
	  phiIndicator(0,nodeNo) = indicator0Node/rhoTotNode;
	  
	  waterChPot(0, nodeNo) = indicator0Node;

	}

        //  MPI: COMMUNCATE SCALAR 'cField'
        mpiBoundary.communciateScalarField(cField);
	
	//  MPI: COMMUNCATE SCALAR 'waterChPot'
        mpiBoundary.communciateScalarField(waterChPot);

	//  MPI: COMMUNCATE SCALAR 'indField'
        mpiBoundary.communciateScalarField(indField);

	//  MPI: COMMUNCATE SCALAR 'rhoDiff'
        mpiBoundary.communciateScalarField(diffWField);
	mpiBoundary.communciateScalarField(phiDiff0);
	mpiBoundary.communciateScalarField(phiDiff1);
	mpiBoundary.communciateScalarField(phiDiff2);
	mpiBoundary.communciateScalarField(phiIndicator);

	//mpiBoundary.communciateScalarField(rhoDiff);

	// Curvature kappa calculation
	// color gradient Norm
	//------------------------------------------
	for (auto nodeNo: bulkNodes) {
	  // -- calculate the color gradient
	  std::valarray<lbBase_t> colorGradNode = grad(cField, 0, nodeNo, grid);
	  gradients.set(1, nodeNo) = colorGradNode;
	  cgNormField(0, nodeNo) = vecNorm<LT>(colorGradNode);
	  //nablaSquared c
	  lbBase_t divGradcNode = divGrad(cField, 0, nodeNo, grid);
	  divGradColorField(0, nodeNo) = divGradcNode;
	  
	}
	mpiBoundary.communciateScalarField(divGradColorField);
	
	for (auto nodeNo: bulkNodes) {
	  std::valarray<lbBase_t> colorGradNode = gradients(1, nodeNo);
	  const lbBase_t rhoTotNode = rho(0, nodeNo);

	  
	  
	  diffgradients.set(0, nodeNo) = grad(phiDiff0, 0, nodeNo, grid);
	  diffgradients.set(1, nodeNo) = grad(phiDiff1, 0, nodeNo, grid);
	  diffgradients.set(2, nodeNo) = grad(phiDiff2, 0, nodeNo, grid);
	  diffgradients.set(3, nodeNo) = grad(phiIndicator, 0, nodeNo, grid);
	  
	  

	  //nablaSquared phiB
	  lbBase_t divGradphiBNode = divGrad(phiIndicator, 0, nodeNo, grid)+divGrad(phiDiff0, 0, nodeNo, grid);
	  divGradPhiBField(0, nodeNo) = divGradphiBNode;
	  
	  
	  //END nablaSquared c

	  
	  gradients.set(2, nodeNo) = grad(indField, 0, nodeNo, grid);
	  gradients.set(3, nodeNo) = grad(diffWField, 0, nodeNo, grid);
	 
	  
	  // Calculate 1st kappa kernel
	  
	  cgFieldX(0, nodeNo) = colorGradNode[0];
	  cgFieldY(0, nodeNo) = colorGradNode[1];
	  cgFieldZ(0, nodeNo) = colorGradNode[2];

	  lbBase_t CGNorm = cgNormField(0, nodeNo);
	  cgFieldNormX(0, nodeNo) = colorGradNode[0]/(CGNorm + (CGNorm < lbBaseEps));
	  cgFieldNormY(0, nodeNo) = colorGradNode[1]/(CGNorm + (CGNorm < lbBaseEps));
	  cgFieldNormZ(0, nodeNo) = colorGradNode[2]/(CGNorm + (CGNorm < lbBaseEps));

	  

	  
	  
	  
	}
	
	for (auto nodeNo: solidBnd) {
	  
	  // Calculate 1st kappa kernel
	  cgNormField(0, nodeNo) = 0.0;
	  
	  cgFieldNormX(0, nodeNo) = 0.0;
	  cgFieldNormY(0, nodeNo) = 0.0;
	  cgFieldNormZ(0, nodeNo) = 0.0;

	  cgFieldX(0, nodeNo) = 0.0;
	  cgFieldY(0, nodeNo) = 0.0;
	  cgFieldZ(0, nodeNo) = 0.0;
	}
	
	mpiBoundary.communciateScalarField(cgNormField);
	mpiBoundary.communciateScalarField(cgFieldNormX);
	mpiBoundary.communciateScalarField(cgFieldNormY);
	mpiBoundary.communciateScalarField(cgFieldNormZ);

	
	
	
	//------------------------------------------
	for (auto nodeNo: bulkNodes) {
	  // -- 1st kappa calculation --------------------------
	  //std::valarray<lbBase_t> gradColorGradNormNode = grad(cgNormField, 0, nodeNo, grid);
	  std::valarray<lbBase_t> colorGradNode = gradients(1, nodeNo);
	  lbBase_t CGNorm = cgNormField( 0, nodeNo);

	  
	  std::valarray<lbBase_t> gradColorGradNormXNode = grad(cgFieldNormX, 0, nodeNo, grid);
	  std::valarray<lbBase_t> gradColorGradNormYNode = grad(cgFieldNormY, 0, nodeNo, grid);
	  std::valarray<lbBase_t> gradColorGradNormZNode = grad(cgFieldNormZ, 0, nodeNo, grid);
	  
	  kappaField(0, nodeNo) = -(gradColorGradNormXNode[0]+gradColorGradNormYNode[1]+gradColorGradNormZNode[2]);
	  
	  
	  //uten denne får man generering av bevegelsesmengde
	  if(sqrt(kappaField(0, nodeNo)*kappaField(0, nodeNo))>0.5 || CGNorm<1e-4)
	    //if(CGNorm<1e-4)
	    kappaField(0, nodeNo) = 0;

	  // -- Alternative kappa calculation --------------------------
	  lbBase_t divGradcNode = divGrad(cField, 0, nodeNo, grid);
	  std::valarray<lbBase_t> gradAbsColorGradNode = grad(cgNormField, 0, nodeNo, grid);

	  kappaField(1, nodeNo) = (LT::dot(colorGradNode/(CGNorm + (CGNorm < lbBaseEps)),gradAbsColorGradNode)-divGradcNode)/(CGNorm + (CGNorm < lbBaseEps));
	  if(CGNorm<lbBaseEps){
	    //kappaField(0, nodeNo) = 0;
	    kappaField(1, nodeNo) = 0;
	  }
	  
	  //if(sqrt(kappaField(0, nodeNo)*kappaField(0, nodeNo))>1)
	  //  kappaField(0, nodeNo) = 1;

	  if(sqrt(kappaField(1, nodeNo)*kappaField(1, nodeNo))>1)
	    kappaField(1, nodeNo) = 1;
	}

	//mpiBoundary.communciateScalarField(kappaField);

	for (auto nodeNo: bulkNodes) {
	  // -- calculate the color gradient
	  std::valarray<lbBase_t> kappaGradNode = grad(kappaField, 0, nodeNo, grid);
	  gradients.set(4, nodeNo) = kappaGradNode;
	}
		
	/*
	for (auto nodeNo: bulkNodes) {
	  
	  lbBase_t indicator0Node = indField(0, nodeNo);
	  lbBase_t rhoDiff0Node = rhoDiff(0, nodeNo);
	  lbBase_t rhoDiff1Node = rhoDiff(1, nodeNo);
	  lbBase_t rhoTotNode = rho(0, nodeNo);
	  lbBase_t CGNorm = cgNormField( 0, nodeNo);
	  lbBase_t cIndNode = indicator0Node/rhoTotNode;
	  lbBase_t c1Node = rhoDiff1Node/rhoTotNode;

	  const lbBase_t rhoType0Node = indicator0Node+rhoDiff0Node;
	  const lbBase_t rhoType1Node = rhoTotNode-rhoType0Node;
	    
	  const lbBase_t cType0Node = rhoType0Node/rhoTotNode;
	  const lbBase_t cType1Node = rhoType1Node/rhoTotNode;
	  	  
	  waterChPot(0, nodeNo) += sigma*divGradColorField(0, nodeNo)*(indicator0Node-rhoDiff1Node/H);
	}
       
	mpiBoundary.communciateScalarField(waterChPot);
	*/


	for (auto nodeNo: bulkNodes) {
	  // -- calculate the water gradient
	  std::valarray<lbBase_t> waterChPotGradNode = grad(waterChPot, 0, nodeNo, grid);
	  
	  lbBase_t rhoDiff0Node = rhoDiff(0, nodeNo);
	  lbBase_t rhoDiff1Node = rhoDiff(1, nodeNo);
	  lbBase_t rhoDiff2Node = rhoDiff(2, nodeNo);
	  lbBase_t indicator0Node = indField(0, nodeNo);
	  
	  lbBase_t rhoTotNode = rho(0, nodeNo);

	  lbBase_t cIndNode = indicator0Node/rhoTotNode;
	  lbBase_t c0Node = rhoDiff0Node/rhoTotNode;
	  lbBase_t c1Node = rhoDiff1Node/rhoTotNode;
	  lbBase_t c2Node = rhoDiff2Node/rhoTotNode;

	  const lbBase_t rhoType0Node = indicator0Node+rhoDiff0Node;
	  const lbBase_t rhoType1Node = rhoTotNode-rhoType0Node-rhoDiff1Node;
	  const lbBase_t rhoColorNode = rhoType0Node + rhoType1Node;
	  const lbBase_t cType0Node = rhoType0Node/rhoColorNode;
	  const lbBase_t cType1Node = 1-cType0Node;
	  
	  
	  /*
	  const lbBase_t rhoType0Node = indicator0Node+rhoDiff0Node;
	  const lbBase_t rhoType1Node = rhoTotNode-rhoDiff1Node-rhoType0Node;
	  const lbBase_t cType0Node = rhoType0Node/(rhoTotNode-rhoDiff1Node);
	  */
	  
	  std::valarray<lbBase_t> colorGradNode = gradients(1, nodeNo);
	  lbBase_t CGNorm = cgNormField( 0, nodeNo);

	  //waterChPotGradNode += cIndNode*LT::c2Inv*0.5*sigma*kappaField(0, nodeNo)*colorGradNode;
	  
	  //waterGradNode -= LT::c2*beta*waterChPot(0, nodeNo)*(1-waterChPot(0, nodeNo))*colorGradNode/(CGNorm + (CGNorm < lbBaseEps));
	  gradients.set(0, nodeNo) = waterChPotGradNode;
	  
	  
	  
	  lbBase_t R0Node = 0.0;
	  lbBase_t R1Node = 0.0;                                    //<------------------------------ Diffusive source set to zero
	  //R1Node = -kinConst*LT::dot(waterChPotGradNode,colorGradNode/(CGNorm + (CGNorm < lbBaseEps)));
	    
	  
	  
	  //R1Node = kinConst*LT::dot(waterChPotGradNode,colorGradNode*0.5);

	  //R1Node = -1e-1*kinConst*waterChPot(0, nodeNo)*CGNorm*0.5;
	  
	  //R1Node = kinConst*LT::dot(waterChPotGradNode,colorGradNode/(CGNorm + (CGNorm < lbBaseEps)));
	  
	  //R1Node = 1e-3*kinConst*waterChPot(0, nodeNo)*CGNorm*0.5*kappaField(0, nodeNo)/(sqrt(kappaField(0, nodeNo)*kappaField(0, nodeNo))+(sqrt(kappaField(0, nodeNo)*kappaField(0, nodeNo))<lbBaseEps));
	  /*
	  if(cType0Node>0.99 && i>=1000){
	    lbBase_t rhoWaterNode= indicator0Node + rhoDiff1Node;	    
	    //R1Node = kinConst*(H*rhoWaterNode-c1Node);
	    R1Node = 2*(rhoTotNode*H*rhoWaterNode - rhoDiff(1, nodeNo));
	  }
	  */

	  /*
	  lbBase_t tmp_thresh = 0.75;
	  if(i>=1000 && cType0Node>tmp_thresh && cType0Node<(1-1e-4)){
	    lbBase_t rhoWaterNode= indicator0Node + rhoDiff1Node;	    
	    //R1Node = kinConst*(H*rhoWaterNode-c1Node);
	    //lbBase_t c1tilde = H*(rhoWaterNode+LT::c2Inv*sigma*kappaField(0, nodeNo)*(1-cType0Node*cType0Node));
	    lbBase_t c1tilde = H*(1+LT::c2Inv*sigma*kappaField(0, nodeNo));
	    R1Node = 2*(rhoTotNode*c1tilde - rhoDiff(1, nodeNo));
	  }
	  */

	  /*
	  //if(cType0Node>0.99 && cType0Node<0.9999 && cType1Node<=0.01 && cType1Node>=0.0001 && i>=1000){
	  if(cType0Node>0.85 && cType0Node<0.999 && i>=1000){  
	    lbBase_t rhoWaterNode= indicator0Node + rhoDiff1Node;
	    //lbBase_t c1tilde = H*(rhoWaterNode/cType0Node+LT::c2Inv*sigma*kappaField(0, nodeNo)*(1-cType0Node));
	    lbBase_t c1tilde = H*rhoWaterNode/cType0Node;
	    R1Node = 2*(rhoTotNode*c1tilde - rhoDiff(1, nodeNo));
	  }
	  else if(cType0Node>=0.999){
	    lbBase_t c1tilde = 0;
	    R1Node = 2*(rhoTotNode*c1tilde - rhoDiff(1, nodeNo));
	  }
	  */
	  
	  R(0,nodeNo)+=R1Node;
	  lbBase_t R2Node = 0.0; 
	  lbBase_t RIndNode =-R1Node;

	  if(rhoDiff1Node+R1Node<0.0 || indicator0Node+RIndNode<0.0){
	    R1Node = 0.0;
	    RIndNode = 0.0;
	  }
	  
	  
	  rhoDiff(0, nodeNo) = rhoDiff0Node += 0.5*R0Node;
	  rhoDiff(1, nodeNo) = rhoDiff1Node += 0.5*R1Node;
	  rhoDiff(2, nodeNo) = rhoDiff2Node += 0.5*R2Node;
	  indicator0Node = indField(0, nodeNo)+= 0.5*RIndNode;
	  
	}


	//----------------------------------Set flux---------------------------------------
	int cartDir = 0;
	bodyForce(0,0,cartDir)=0.0; //calcFluxForceCartDir<LT>(0, f, bulkNodes, cartDir, input["fluid"]["momx"], numNodesGlobal);
	
	//----------------------------------Set Capillary Number---------------------------------------
	/*
	//lbBase_t meanfcX=0.0;
	lbBase_t meancType0fcX=0.0;
	lbBase_t meancType1fcX=0.0;
	//lbBase_t meanfcXGlobal;
	lbBase_t meancType0fcXGlobal;
	lbBase_t meancType1fcXGlobal;

	lbBase_t meancType0=0.0;
	lbBase_t meancType1=0.0;
	lbBase_t meancType0Global;
	lbBase_t meancType1Global;
	
	for (auto nodeNo : bulkNodes) {
	  std::valarray<lbBase_t> fTot = f(0, nodeNo);
	  std::valarray<lbBase_t> sumfc = LT::qSumC(fTot);

	  lbBase_t indicator0Node = indField(0, nodeNo);
	  lbBase_t rhoDiff0Node = rhoDiff(0, nodeNo); 
	  lbBase_t rhoTotNode = rho(0, nodeNo);   
	 
	  const lbBase_t rhoType0Node = indicator0Node+rhoDiff0Node;
	  const lbBase_t rhoType1Node = rhoTotNode-rhoType0Node;
	  
	  const lbBase_t cType0Node = rhoType0Node/rhoTotNode;
	  const lbBase_t cType1Node = rhoType1Node/rhoTotNode;

	  //meanfcX+=sumfc[0];
	  meancType0fcX+=cType0Node*sumfc[0];
	  meancType1fcX+=cType1Node*sumfc[0];

	  meancType0+=cType0Node;
	  meancType1+=cType1Node;
	}
	
	//MPI_Allreduce(&meanfcX, &meanfcXGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&meancType0fcX, &meancType0fcXGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&meancType1fcX, &meancType1fcXGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&meancType0, &meancType0Global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&meancType1, &meancType1Global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	

	//meanfcXGlobal/=numNodesGlobal;
	meancType0fcXGlobal/=numNodesGlobal;
	meancType1fcXGlobal/=numNodesGlobal;

	meancType0Global/=numNodesGlobal;
	meancType1Global/=numNodesGlobal;

	lbBase_t CapNo=input["fluid"]["Ca"];
	
	lbBase_t fluxForceX=2*(sigma*CapNo - meancType0fcXGlobal/nu0Inv-meancType1fcXGlobal/nu1Inv)/(meancType0Global/nu0Inv+meancType1Global/nu1Inv);
	
	bodyForce(0,0,0)=fluxForceX;
	//----------------------------------end Flux---------------------------------------
	*/
	//----------------------------------Set Capillary Number---------------------------------------
	/*
	int cartDir = 0;
	lbBase_t CapNumb=input["fluid"]["Ca"];

	bodyForce(0,0,cartDir)= calcCapNumbForceCartDir(0, f, rho, bulkNodes, cartDir, sigma*CapNumb,  1./nu0Inv, 1./nu1Inv, numNodesGlobal);
	*/
	//----------------------------------end Flux---------------------------------------
	
        for (auto nodeNo: bulkNodes) {

            // Set the local total lb distribution
            // Use std:valarray insatead of auto as valarray uses expression templates that do not work well with auto
	    std::valarray<lbBase_t> fTot = f(0, nodeNo);
	    
	    lbBase_t indicator0Node = indField(0, nodeNo);
	    
	    // -- total density
	    lbBase_t rhoTotNode = rho(0, nodeNo);
	    
	    lbBase_t rhoDiff0Node = rhoDiff(0, nodeNo);
	    lbBase_t rhoDiff1Node = rhoDiff(1, nodeNo);
	    lbBase_t rhoDiff2Node = rhoDiff(2, nodeNo);

	    std::valarray<lbBase_t> colorGradNode = gradients(1, nodeNo);
	    lbBase_t CGNorm = cgNormField( 0, nodeNo);//vecNorm<LT>(colorGradNode);

	    const lbBase_t c0Node = rhoDiff0Node/rhoTotNode;
	    const lbBase_t c1Node = rhoDiff1Node/rhoTotNode;
	    const lbBase_t c2Node = rhoDiff2Node/rhoTotNode;
	    const lbBase_t cIndNode = indicator0Node/rhoTotNode;

	    /*
	    const lbBase_t rhoType0Node = indicator0Node+rhoDiff0Node;
	    const lbBase_t rhoType1Node = rhoTotNode-rhoType0Node;
	    
	    const lbBase_t cType0Node = rhoType0Node/rhoTotNode;
	    const lbBase_t cType1Node = rhoType1Node/rhoTotNode;
	    */

	    
	    const lbBase_t rhoType0Node = indicator0Node+rhoDiff0Node;
	    const lbBase_t rhoType1Node = rhoTotNode-rhoType0Node-rhoDiff1Node;
	    const lbBase_t rhoColorNode = rhoType0Node + rhoType1Node;
	    const lbBase_t cType0Node = rhoType0Node/rhoColorNode;
	    const lbBase_t cType1Node = 1-cType0Node;
	    
	    
	    /*
	    const lbBase_t rhoType0Node = indicator0Node+rhoDiff0Node;
	    const lbBase_t rhoType1Node = rhoTotNode-rhoDiff1Node-rhoType0Node;
	    const lbBase_t cType0Node = rhoType0Node/(rhoTotNode-rhoDiff1Node);
	    const lbBase_t cType1Node = 1-cType0Node;
	    */
	    // -----------------------------------------------
	    
            // -- force
            //std::valarray<lbBase_t> forceNode = setForceGravity(rhoTotNode*indicator0Node, rhoTotNode*(1-indicator0Node), bodyForce, 0);
	    //std::valarray<lbBase_t> IFTforceNode = 0*0.5*sigma*kappaField(0, nodeNo)*colorGradNode;
	    //std::valarray<lbBase_t> IFTforceNode = sigma*kappaField(0, nodeNo)*beta*cType0Node*cType1Node*colorGradNode/(CGNorm + (CGNorm < lbBaseEps));
	    
	    /*
	    lbBase_t ksiPhaseField = 4/beta;
	    lbBase_t betaPhaseField = 12*sigma/ksiPhaseField;
	    lbBase_t kappaPhaseField = 1.5*sigma*ksiPhaseField;
	    std::valarray<lbBase_t> IFTforceNode = (4*betaPhaseField*cType0Node*cType1Node*(cType0Node-0.5)
						    -kappaPhaseField*divGradPhiBField(0, nodeNo))*0.5*colorGradNode;
	    */
	    //std::valarray<lbBase_t> IFTforceNode = cType0Node*gradients(0, nodeNo);
	    //std::valarray<lbBase_t> IFTforceNode =waterChPot(0, nodeNo)*0.5*colorGradNode;
	    
	    //std::valarray<lbBase_t> kappaGradNode = gradients.set(4, nodeNo);
	    //IFTforceNode -= colorGradNode*(1-4*beta*cType0Node*cType1Node/(CGNorm + (CGNorm < lbBaseEps)));
	    //std::valarray<lbBase_t> IFTforceNode = 2/beta*sigma*kappaField(0, nodeNo)*CGNorm*colorGradNode;
	    //std::valarray<lbBase_t> forceNode = IFTforceNode;
	    std::valarray<lbBase_t> forceNode = bodyForce(0,0);
	    //std::valarray<lbBase_t> forceNodeDiff = 0.*IFTforceNode;//0.0*colorGradNode;
	    //lbBase_t absKappa = sqrt(kappaField(0, nodeNo)*kappaField(0, nodeNo));
	    //lbBase_t pGradNorm = 0.5*sigma*absKappa*CGNorm;
	    //std::valarray<lbBase_t> gradP_alt = 0.5*sigma*absKappa*colorGradNode;
            // -- velocity
            std::valarray<lbBase_t> velNode = calcVel<LT>(fTot, rhoTotNode, forceNode);  // LBmacroscopic
            vel.set(0, nodeNo) = velNode;
	    //std::valarray<lbBase_t> divGradFNode = divGradF(0, nodeNo);
	    
	    

	    const lbBase_t q0Node = Q(0, nodeNo);
	    const lbBase_t q1Node = 0.0;// Q(1, nodeNo);
	    
	    
	    
	    
            // CALCULATE BGK COLLISION TERM
            // Mean collision time /rho_tot/\nu_tot = \sum_s \rho_s/\nu_s
	    lbBase_t tau = LT::c2Inv * rhoTotNode / ((indicator0Node)*nu0Inv + (rhoTotNode-indicator0Node-rhoDiff0Node)*nu1Inv) + 0.5;

	    

	    lbBase_t tauD_eff = LT::c2Inv / ((cIndNode+c0Node)*DwInv + (1-cIndNode-c0Node)*DdwInv) + 0.5;
	    
	    
	    
            lbBase_t uu = LT::dot(velNode, velNode);  // Square of the velocity
            std::valarray<lbBase_t> cu = LT::cDotAll(velNode);  // velocity dotted with lattice vectors
            std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fTot, tau, rhoTotNode, uu, cu);  // LBcollision

	    
	    std::valarray<lbBase_t> omegaBGKInd   = calcOmegaBGK<LT>(gInd(0, nodeNo), tauD_eff, indicator0Node, uu, cu);  // LBcollisionIndicator
	    std::valarray<lbBase_t> omegaBGKDiff0 = calcOmegaBGK<LT>(g(0,nodeNo),     tauD_eff, rhoDiff0Node,   uu, cu);  // LBcollision
	    std::valarray<lbBase_t> omegaBGKDiff1 = calcOmegaBGK<LT>(g(1,nodeNo),     tauD_eff, rhoDiff1Node,   uu, cu);  // LBcollision
	    std::valarray<lbBase_t> omegaBGKDiff2 = calcOmegaBGK<LT>(g(2,nodeNo),     tauD_eff, rhoDiff2Node,   uu, cu);  // LBcollision
	    
	    

	    // CALCULATE BGK COLLISION TERM FOR ORIGINAL COLOR GRADIENT
            // Mean collision time /rho_tot/\nu_tot = \sum_s \rho_s/\nu_s
	    
	    
            // CALCULATE FORCE CORRECTION TERM
            lbBase_t  uF = LT::dot(velNode, forceNode);
            std::valarray<lbBase_t>  cF = LT::cDotAll(forceNode);
	    //std::valarray<lbBase_t>  cgradP_alt = LT::cDotAll(gradP_alt);
	    //std::valarray<lbBase_t>  cDivGradF = LT::cDotAll(divGradFNode);
            std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);  // LBcollision

	    /*
	    uF = LT::dot(velNode, 1.5*forceNode);
	    cF = LT::cDotAll(1.5*forceNode);
	    */
	    
	    
	    std::valarray<lbBase_t> deltaOmegaFDiffInd = calcDeltaOmegaFDiff<LT>(tauD_eff, cIndNode, cu, uF, cF);  // LBcollision
	    std::valarray<lbBase_t> deltaOmegaFDiff0   = calcDeltaOmegaFDiff<LT>(tauD_eff, c0Node,   cu, uF, cF);  // LBcollision
	    std::valarray<lbBase_t> deltaOmegaFDiff1   = calcDeltaOmegaFDiff<LT>(tauD_eff, c1Node,   cu, -uF, -cF);  // LBcollision
	    std::valarray<lbBase_t> deltaOmegaFDiff2   = calcDeltaOmegaFDiff<LT>(tauD_eff, c2Node,   cu, -uF, -cF);  // LBcollision

            // CALCULATE MASS SOURCE CORRECTION TERM
            std::valarray<lbBase_t> deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tau, cu, uu, q0Node);
            //std::valarray<lbBase_t> deltaOmegaQ1 = calcDeltaOmegaQ<LT>(tau, cu, uu, q1Node);

	    lbBase_t  uCG = LT::dot(velNode, colorGradNode);
	    
	    // CALCULATE DIFFUSIVE MASS SOURCE CORRECTION TERM

	    std::valarray<lbBase_t> deltaOmegaR1   = calcDeltaOmegaR<LT>(tauD_eff, cu, 0*R(0, nodeNo));
	    //std::valarray<lbBase_t> deltaOmegaRInd = calcDeltaOmegaR<LT>(tauD_eff, cu, -R(0, nodeNo));
	    std::valarray<lbBase_t> deltaOmegaRInd = calcDeltaOmegaR<LT>(tauD_eff, cu, R(0, nodeNo));
	    
	    
	    
            // CALCULATE SURFACE TENSION PERTURBATION

	    // -- calculate the normalized color gradient
	    
	    //std::valarray<lbBase_t> colorGradNode = gradients.set(1, nodeNo);

	    std::valarray<lbBase_t> cCG = LT::cDotAll(colorGradNode);
	    
            colorGradNode *= 1.0/(CGNorm + (CGNorm < lbBaseEps));

	    
	    

            std::valarray<lbBase_t> cCGNorm = LT::cDotAll(colorGradNode);
	 
	    
	    //lbBase_t CGNormF = LT::dot(colorGradNode, IFTforceNode);

	    //lbBase_t divGradRhoNode = 3*divGradPField(0, nodeNo);
	    
	    lbBase_t uCGNorm = LT::dot(colorGradNode, velNode);

	    /*
            std::valarray<lbBase_t> deltaOmegaST = calcDeltaOmegaST<LT>(tau, sigma, CGNorm, cCGNorm);
	    */
	    
	    
	    std::valarray<lbBase_t> deltaOmegaRCInd   = calcDeltaOmegaRC<LT>(beta, indicator0Node, cType1Node, 1, cCGNorm);
	   
	    std::valarray<lbBase_t> deltaOmegaRCDiff0 = calcDeltaOmegaRC<LT>(beta, rhoDiff0Node, cType1Node, 1, cCGNorm);
	    
	    std::valarray<lbBase_t> deltaOmegaRCDiff1 = calcDeltaOmegaRC<LT>(beta, 0*rhoDiff1Node, (1-cType1Node), 1, -cCGNorm);
	    
	    std::valarray<lbBase_t> deltaOmegaRCDiff2 = calcDeltaOmegaRC<LT>(beta, rhoDiff2Node,  cType0Node, 1, -cCGNorm);
	    
	    
	    /*
	    std::valarray<lbBase_t> deltaOmegaRCInd   = calcDeltaOmegaRC2<LT>(beta*(1-0.5/tauD_eff), indicator0Node, cType1Node, 1, cCGNorm);
	   
	    std::valarray<lbBase_t> deltaOmegaRCDiff0 = calcDeltaOmegaRC2<LT>(beta*(1-0.5/tauD_eff), rhoDiff0Node, cType1Node, 1, cCGNorm);
	    
	    std::valarray<lbBase_t> deltaOmegaRCDiff1 = calcDeltaOmegaRC2<LT>(beta*(1-0.5/tauD_eff), 0*rhoDiff1Node, (1-cType1Node), 1, -cCGNorm);
	    
	    std::valarray<lbBase_t> deltaOmegaRCDiff2 = calcDeltaOmegaRC2<LT>(beta*(1-0.5/tauD_eff), rhoDiff2Node,  cType0Node, 1, -cCGNorm);
	    */
	    
	    
	    // CALCULATE SURFACE TENSION PERTURBATION FOR ORIGINAL COLOR GRADIENT
            // -- calculate the normalized color gradient
	    
	    std::valarray<lbBase_t> deltaOmegaST = calcDeltaOmegaST<LT>(tau, sigma, CGNorm, cCGNorm);
	    //std::valarray<lbBase_t> deltaOmegaST = calcDeltaOmegaST2<LT>(tau, sigma, cu, CGNorm, cCGNorm);
	    /*
	    std::valarray<lbBase_t> capTensorLowTri = tau*LT::qSumCCLowTri(deltaOmegaST);
	    VectorField<LT> Tgradphi(4, 1);
	    

	    for (int p = 0; p < 4; ++p) {
	      Tgradphi(p,0)[0] = capTensorLowTri[0]*diffgradients(p, nodeNo)[0]
		+capTensorLowTri[1]*diffgradients(p, nodeNo)[1] + capTensorLowTri[3]*diffgradients(p, nodeNo)[2];
 	      Tgradphi(p,0)[1] = capTensorLowTri[1]*diffgradients(p, nodeNo)[0]
		+capTensorLowTri[2]*diffgradients(p, nodeNo)[1] + capTensorLowTri[4]*diffgradients(p, nodeNo)[2];
	      Tgradphi(p,0)[2] = capTensorLowTri[3]*diffgradients(p, nodeNo)[0]
		+capTensorLowTri[4]*diffgradients(p, nodeNo)[1] + capTensorLowTri[5]*diffgradients(p, nodeNo)[2];
	    }

	    std::valarray<lbBase_t>  cTgradphi0 = LT::cDotAll(Tgradphi(0,0));
	    std::valarray<lbBase_t>  cTgradphi1 = LT::cDotAll(Tgradphi(1,0));
	    std::valarray<lbBase_t>  cTgradphi2 = LT::cDotAll(Tgradphi(2,0));
	    std::valarray<lbBase_t>  cTgradphi3 = LT::cDotAll(Tgradphi(3,0));
	    
	    std::valarray<lbBase_t> deltaOmegaSTDiffCorr0= calcDeltaOmegaSTDiffCorr<LT>(tauD_eff, cTgradphi0);
	    std::valarray<lbBase_t> deltaOmegaSTDiffCorr1= calcDeltaOmegaSTDiffCorr<LT>(tauD_eff, cTgradphi1);
	    std::valarray<lbBase_t> deltaOmegaSTDiffCorr2= calcDeltaOmegaSTDiffCorr<LT>(tauD_eff, cTgradphi2);
	    std::valarray<lbBase_t> deltaOmegaSTDiffCorrInd= calcDeltaOmegaSTDiffCorr<LT>(tauD_eff, cTgradphi3);
	    
	    lbBase_t  uGrad0 = LT::dot(velNode, rhoTotNode*diffgradients(0, nodeNo));
            std::valarray<lbBase_t>  cGrad0 = LT::cDotAll(rhoTotNode*diffgradients(0, nodeNo));
	    lbBase_t  uGrad1 = LT::dot(velNode, rhoTotNode*diffgradients(1, nodeNo));
            std::valarray<lbBase_t>  cGrad1 = LT::cDotAll(rhoTotNode*diffgradients(1, nodeNo));
	    lbBase_t  uGrad2 = LT::dot(velNode, rhoTotNode*diffgradients(2, nodeNo));
            std::valarray<lbBase_t>  cGrad2 = LT::cDotAll(rhoTotNode*diffgradients(2, nodeNo));
	    lbBase_t  uGrad3 = LT::dot(velNode, rhoTotNode*diffgradients(3, nodeNo));
            std::valarray<lbBase_t>  cGrad3 = LT::cDotAll(rhoTotNode*diffgradients(3, nodeNo));
	    */
	    
	    // CALCULATE DIFFUSIVE GRAD CORRECTION TERM
	    /*
	    std::valarray<lbBase_t> deltaOmegaGradDiff0   = calcDeltaOmegaGradDiff<LT>(tauD_eff, cu, uGrad0, cGrad0);	    
	    std::valarray<lbBase_t> deltaOmegaGradInd = calcDeltaOmegaGradDiff<LT>(tauD_eff, cu, uGrad3, cGrad3);
	    std::valarray<lbBase_t> deltaOmegaGradDiff1   = calcDeltaOmegaGradDiff<LT>(tauD_eff, cu, uGrad1, cGrad1);
	    std::valarray<lbBase_t> deltaOmegaGradDiff2 = calcDeltaOmegaGradDiff<LT>(tauD_eff, cu, uGrad2, cGrad2);
	    */
	    

	    
            // COLLISION AND PROPAGATION
	  
       
	    //-----------------------
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be
	      
	      fTmp(0, q,  grid.neighbor(q, nodeNo)) =    fTot[q]            + omegaBGK[q]      + deltaOmegaF[q]  +  deltaOmegaST[q];
                
	      gIndTmp(0, q,  grid.neighbor(q, nodeNo)) = gInd(0, q, nodeNo) + omegaBGKInd[q]   + deltaOmegaFDiffInd[q] /*+ cIndNode*deltaOmegaF[q]*/ + deltaOmegaRCInd[q]   + deltaOmegaRInd[q] + cIndNode*(tau/tauD_eff)*deltaOmegaST[q]/* + deltaOmegaSTDiffCorrInd[q]*/;
	      gTmp(0, q,  grid.neighbor(q, nodeNo)) =    g(0, q, nodeNo)    + omegaBGKDiff0[q] + deltaOmegaFDiff0[q]   /*+ c0Node*deltaOmegaF[q]*/ + deltaOmegaRCDiff0[q]                     + c0Node*(tau/tauD_eff)*deltaOmegaST[q]/* + deltaOmegaSTDiffCorr0[q]*/;
	      gTmp(1, q,  grid.neighbor(q, nodeNo)) =    g(1, q, nodeNo)    + omegaBGKDiff1[q] + deltaOmegaFDiff1[q]   /*+ c1Node*deltaOmegaF[q]*/ + deltaOmegaRCDiff1[q]   + deltaOmegaR1[q] + c1Node*(tau/tauD_eff)*deltaOmegaST[q]/* + deltaOmegaSTDiffCorr1[q]*/;
	      gTmp(2, q,  grid.neighbor(q, nodeNo)) =    g(2, q, nodeNo)    + omegaBGKDiff2[q] + deltaOmegaFDiff2[q]   /*+ c2Node*deltaOmegaF[q]*/ + deltaOmegaRCDiff2[q]                     + c2Node*(tau/tauD_eff)*deltaOmegaST[q]/* + deltaOmegaSTDiffCorr2[q]*/;
	      /*
	      if (gTmp(0, q,  grid.neighbor(q, nodeNo)) <0){
		std::cout<<"Negative distribution g="<< gTmp(0, q,  grid.neighbor(q, nodeNo))<<std::endl;
		}
	      */
	      //-----------------------
            }
        } // End nodes

        // PRINT

        if ( (i % input["iterations"]["write"]) == 0) {
	  
	  
	  //-------------------------------------------------------------
	  /*
	  std::string tmpName(outDir);
          tmpName += std::to_string(myRank) + "_" + std::to_string(i);
          tmpName += ".dat";
          std::ofstream ofs;
          ofs.open(tmpName);
          if (!ofs) {
            std::cout << "Error: could not open file: " << tmpName << std::endl;
            return 1;
          }

          for (auto nodeNo: bulkNodes) {
            ofs << std::setprecision(23) << grid.pos(nodeNo, 0) << " " << grid.pos(nodeNo, 1) << " " << grid.pos(nodeNo, 2) << " " << rho(0, nodeNo) << " " << rhoDiff(0, nodeNo)
                            << " " << vel(0, 0, nodeNo) << " " << vel(0, 1, nodeNo) << " " << vel(0, 2, nodeNo) << " " << nodes.getType(nodeNo) << std::endl;
          }
          ofs.close();
	  */
	  //-------------------------------------------------------------

	  
          // JLV
          output.write(i);
          // JLV
	  if (myRank==0)
            std::cout << "PLOT AT ITERATION : " << i << std::endl;
        }

	if ( (i % 100) == 0 && myRank==0){
	  std::ofstream ofs;
	  std::string tmpName(outDir2+"/force.dat");
	  ofs.open(tmpName, std::ios::app);
	  if (!ofs) {
	    std::cout << "Error: could not open file: " << tmpName << std::endl;
	    return 1;
	  }
	  ofs << i << " " << std::setprecision(23) << bodyForce(0, 0, 0) << std::endl;
	  
	  ofs.close();
	}

	
        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield
	g.swapData(gTmp);  // LBfield
	gInd.swapData(gIndTmp);  // LBfield
	
	
        // MPI Boundary
        // Hente verdier hente fra ghost
        // Sette i bulk
        mpiBoundary.communicateLbField(0, f, grid);
	/*
	mpiBoundary.communicateLbField(0, g, grid);
        mpiBoundary.communicateLbField(1, g, grid);
	mpiBoundary.communicateLbField(2, g, grid);
	*/	
	mpiBoundary.communicateLbField(g, grid);
	mpiBoundary.communicateLbField(0, gInd, grid);

	
	
        // BOUNDARY CONDITIONS
        bbBnd.apply(0, f, grid);  // LBboundary
	bbBnd.apply(0, g, grid);  // LBboundary
        bbBnd.apply(1, g, grid);  // LBboundary
	bbBnd.apply(2, g, grid);  // LBboundary
	bbBnd.apply(0, gInd, grid); // LBboundary

	
	
    } // End iterations
    // -----------------END MAIN LOOP------------------

    //mpi.end();

    MPI_Finalize();

    return 0;
}

