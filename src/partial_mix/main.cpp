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
#include "./LBpartial_mixHelp.h"

// SET THE LATTICE TYPE
#define LT D2Q9
#define VTK_CELL VTK::pixel
//#define LT D3Q19
//#define VTK_CELL VTK::voxel

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
    // std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
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
    // Set solid boundary 
    std::vector<int> solidBoundaryNodes = findSolidBndNodes(nodes);

    // ***************
    // READ FROM INPUT
    // ***************
    // Number of iterations
    int nIterations = static_cast<int>( input["iterations"]["max"]);
    // Write interval
    int nItrWrite = static_cast<int>( input["iterations"]["write"]);

    // Fluid Flow
    // **********
    // Number of fluid fields
    const int nFluidFields = input["fluid"]["numfluids"];
    std::cout << "nFluidFields = " << nFluidFields << std::endl;
    // Relaxation time
    // lbBase_t tau = 1;
    // Kinematic viscosities
    std::valarray<lbBase_t> kin_visc = inputAsValarray<lbBase_t>(input["fluid"]["viscosity"]);
    std::valarray<lbBase_t> kin_visc_inv(nFluidFields);
    for (int n = 0; n < nFluidFields; ++n){
      kin_visc_inv[n] = 1./kin_visc[n];
    }
    // Body force
    VectorField<LT> bodyForce(1, 1);
    bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);
    std::cout << std::endl;
    //Parameters for compressibility
    std::valarray<lbBase_t> alpha = LT::w0*inputAsValarray<lbBase_t>(input["fluid"]["alpha"]);
    std::cout << "alpha =";
    for (int n = 0; n < nFluidFields; ++n)
        std::cout << " " << alpha[n];
    std::cout << std::endl;
    std::valarray<lbBase_t> Gamma0 = inputAsValarray<lbBase_t>(input["fluid"]["alpha"]);
    std::valarray<lbBase_t> GammaNonZero = (1-LT::w0*Gamma0)/(1-LT::w0);    
    // Phase separation: fluids
    std::valarray<lbBase_t> beta = inputAsValarray<lbBase_t>(input["fluid"]["phaseseparation"]["beta"]);
    std::cout << "beta = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
        for (int j=0; j < nFluidFields; ++j) {
            std::cout << " " << beta[i*nFluidFields + j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    // Surface tension
    std::valarray<lbBase_t> sigma = inputAsValarray<lbBase_t>(input["fluid"]["phaseseparation"]["sigma"]);
    std::cout << "sigma = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
        for (int j=0; j < nFluidFields; ++j) {
            std::cout << " " << sigma[i*nFluidFields + j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Diffusion
    // *********
    // Number of diffusive fields
    const int nDiffFields = input["diffusion"]["numfields"];
    std::cout << "nDiffFields = " << nDiffFields << std::endl;
    // Diffusion coefficients
    std::valarray<lbBase_t> diff_coef = inputAsValarray<lbBase_t>(input["diffusion"]["fluidinteraction"]["diffcoef"]);
    std::valarray<lbBase_t> diff_coef_inv(nDiffFields*nFluidFields);
    std::cout << "Diff. coef = " << std::endl;
    for (int i=0; i < nDiffFields; ++i) {
      std::cout << "diff.field[" << i <<"]: ";
        for (int j=0; j < nFluidFields; ++j) {
	    diff_coef_inv[i*nFluidFields + j] = 1./diff_coef[i*nFluidFields + j];
            std::cout << " " << diff_coef[i*nFluidFields + j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    std::valarray<lbBase_t> betaDiff = inputAsValarray<lbBase_t>(input["diffusion"]["fluidinteraction"]["beta"]);
    std::cout << "betaDiff = " << std::endl;
    for (int k=0; k < nDiffFields; ++k) {
      std::cout << "diff.field[" << k <<"]: "<< std::endl;
      for (int i=0; i < nFluidFields; ++i) {
        for (int j=0; j < nFluidFields; ++j) {
	  std::cout << " " << betaDiff[k*nFluidFields*nFluidFields+i*nFluidFields + j];
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));
    
    std::string outputDir2 = outputDir + "/out" + dirNum;
    //END READ FROM INPUT
    
    
    // ******************
    // MACROSCOPIC FIELDS
    // ******************
    // Fluid Flow
    // **********
    // Density
    ScalarField rho(nFluidFields, grid.size());
    ScalarField rhoRel(nFluidFields, grid.size());
    ScalarField rhoTot(1, grid.size());

    // Initiate density from file
    setScalarAttribute(rho, "init_rho_", vtklb);


    /*
    // Wall wettability
    setScalarAttributeWall(rhoRel, "init_rho_wall_", vtklb, nodes);

    // -- scale wettability
    normelizeScalarField(rhoRel, solidBoundaryNodes);
    */
    
    // Velocity
    VectorField<LT> vel(1, grid.size());
    // Initiate velocity
    
    for (auto nodeNo: bulkNodes) {
      vel.set(0, nodeNo) = 0;
    }
    

    //Phase separation
    //Color gradients for fluid pair
    //Absolute values of color gradients for fluid pair

    ScalarField FNorm(nFluidFields*(nFluidFields-1)/2, grid.size());
    VectorField<LT> F(nFluidFields*(nFluidFields-1)/2, grid.size());
    VectorField<LT> ForceField(1, grid.size());
   
    VectorField<LT> unitNormal(nFluidFields*(nFluidFields-1)/2, grid.size());

    //Curvature of interface for fluid pair
    //(lower triangular matrix)
    ScalarField kappa(nFluidFields*(nFluidFields-1)/2, grid.size());
    
    for (auto fieldNo=0; fieldNo < F.num_fields(); ++fieldNo) {
        for (auto nodeNo: bulkNodes) {
	  F.set(fieldNo, nodeNo) = 0;
	  unitNormal.set(fieldNo, nodeNo) = 0;
	  kappa(fieldNo, nodeNo) = 0;
        }
    }
    for (auto fieldNo=0; fieldNo < ForceField.num_fields(); ++fieldNo) {
        for (auto nodeNo: bulkNodes) {
	  ForceField.set(fieldNo, nodeNo) = 0;
        }
    }
    
    
    // Advection-Diffusion 
    // *******************
    ScalarField rhoD(nDiffFields, grid.size());
    // Initiate diffusive fields
   
    for (auto nodeNo: bulkNodes) {
        for (int fieldNo=0; fieldNo < rhoD.num_fields(); ++fieldNo) {
	  rhoD(fieldNo, nodeNo) = 0.01*rho(0, nodeNo);
	}

    }
    
    // ****************** 
    // SETUP BOUNDARY
    // ******************
    HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

    // *********
    // LB FIELDS
    // *********
    // Fluid Fields
    // ************
    LbField<LT> f(nFluidFields, grid.size()); 
    LbField<LT> fTmp(nFluidFields, grid.size());
    // initiate lb distributions
    for (int fieldNo=0; fieldNo < f.num_fields(); ++fieldNo) {
        for (auto nodeNo: bulkNodes) {
            for (int q = 0; q < LT::nQ; ++q) {
                f(fieldNo, q, nodeNo) = LT::w[q]*rho(fieldNo, nodeNo);
                fTmp(fieldNo, q, nodeNo) = 0;
            }
        }
    }

    LbField<LT> fTot(1, grid.size()); 
    LbField<LT> fTotTmp(1, grid.size());
    // initiate lb distributions
    
    for (auto nodeNo: bulkNodes) {
      rhoTot(0, nodeNo) = 0.0;
      for (int fieldNo=0; fieldNo < f.num_fields(); ++fieldNo) {
	rhoTot(0, nodeNo) += rho(fieldNo, nodeNo);
      }
      for (int q = 0; q < LT::nQ; ++q) {
	fTot(0, q, nodeNo) = LT::w[q]*rhoTot(0, nodeNo);
	fTotTmp(0, q, nodeNo) = 0;
      }
    }
    

    
    // Diffusive Fields
    //*****************
    LbField<LT> g(nDiffFields, grid.size()); 
    LbField<LT> gTmp(nDiffFields, grid.size());
    // initiate lb distributions
    for (int fieldNo=0; fieldNo < g.num_fields(); ++fieldNo) {
        for (auto nodeNo: bulkNodes) {
            for (int q = 0; q < LT::nQ; ++q) {	      
	        g(fieldNo, q, nodeNo) = LT::w[q]*rhoD(fieldNo, nodeNo);
                gTmp(fieldNo, q, nodeNo) = 0;
            }
        }
    }
    
    
    // **********
    // OUTPUT VTK
    // **********
    VTK::Output<VTK_CELL, double> output(VTK::BINARY, grid.getNodePos(bulkNodes), outputDir2, myRank, nProcs);
    output.add_file("lb_run");

    output.add_variable("rhoTot", 1, rhoTot.get_data(), rhoTot.get_field_index(0, bulkNodes));
 
    output.add_variable("vel", LT::nD, vel.get_data(), vel.get_field_index(0, bulkNodes));
    
    for (int fieldNo=0; fieldNo < nFluidFields; ++fieldNo) {
        output.add_variable("rho" + std::to_string(fieldNo), 1, rho.get_data(), rho.get_field_index(fieldNo, bulkNodes));
    }
    

    for (int fieldNo=0; fieldNo < nDiffFields; ++fieldNo) {
        output.add_variable("rhoD" + std::to_string(fieldNo), 1, rhoD.get_data(), rhoD.get_field_index(fieldNo, bulkNodes));
    }
    
    for (int cnt=0; cnt < nFluidFields*(nFluidFields-1)/2; ++cnt) {
        output.add_variable("kappa_" + std::to_string(cnt), 1, kappa.get_data(), kappa.get_field_index(cnt, bulkNodes));
	output.add_variable("FNorm_" + std::to_string(cnt), 1, FNorm.get_data(), FNorm.get_field_index(cnt, bulkNodes));
	//output.add_variable("F" + std::to_string(cnt), LT::nD, F.get_data(), F.get_field_index(cnt, bulkNodes));
    }
    output.add_variable("F", LT::nD, F.get_data(), F.get_field_index(0, bulkNodes));
    output.add_variable("unitNormal", LT::nD, unitNormal.get_data(), unitNormal.get_field_index(0, bulkNodes));
    output.add_variable("forceField", LT::nD, ForceField.get_data(), ForceField.get_field_index(0, bulkNodes));

    
    // *********
    // MAIN LOOP
    // *********
    // Just ad hoc helper fields (Should be set outside the loop structure)
    std::valarray<lbBase_t> cNormInv(LT::nQ);
    std::valarray<lbBase_t> wAll(LT::nQ);
    for (int q=0; q < LT::nQNonZero_; ++q) {
        cNormInv[q] = 1.0/LT::cNorm[q];
        wAll[q] = LT::w[q];
    }
    cNormInv[LT::nQNonZero_] = 0;
    wAll[LT::nQNonZero_] = LT::w[LT::nQNonZero_];
    
    for (int i = 0; i <= nIterations; i++) {
        // Macroscopic values : rho, rhoTot and rhoRel = rho/rhoTot
        calcDensityFields(rho, rhoRel, rhoTot, rhoD, bulkNodes, f, fTot, g);
      
      
        mpiBoundary.communciateScalarField(rhoRel);

	for (auto nodeNo: bulkNodes) {
	  // Cacluate gradient
	  ScalarField rhoRelNode(1, nFluidFields);
	  for (int fieldNo = 0; fieldNo < nFluidFields; ++fieldNo){
	    rhoRelNode(0, fieldNo) = rhoRel(fieldNo, nodeNo);
	  }
	  const CGAttributes<LT> cgat(nFluidFields, nodeNo, cNormInv, Gamma0, rhoRelNode, rhoRel, grid);

	  for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	    FNorm(cnt, nodeNo)= cgat.FNorm_(0, cnt);
	    F.set(cnt, nodeNo)= cgat.F_(0, cnt);
	    unitNormal.set(cnt, nodeNo)= F(cnt, nodeNo)/(FNorm(cnt, nodeNo)+(FNorm(cnt, nodeNo)<lbBaseEps));
	  }
	}

	mpiBoundary.communciateVectorField_TEST(unitNormal);
	//mpiBoundary.communciateVectorField_TEST(F);
	//mpiBoundary.communciateScalarField(FNorm);
	
	for (auto nodeNo: bulkNodes) {

	  
	    // Cacluate gradient
            ScalarField rhoRelNode(1, nFluidFields);
            for (int fieldNo = 0; fieldNo < nFluidFields; ++fieldNo){
                rhoRelNode(0, fieldNo) = rhoRel(fieldNo, nodeNo);
	    }
            //const CGAttributes<LT> cgat(nFluidFields, nodeNo, cNormInv, Gamma0, rhoRelNode, rhoRel, grid);

	    
	    VectorField<LT> IFTforceNode(1,1);
	    IFTforceNode.set(0 ,0) = 0;
	    ForceField.set(0, nodeNo)=0.0;
	    
	    int cnt = 0;
	    for (int fieldNo_k=0; fieldNo_k<nFluidFields; ++fieldNo_k) {
	      for (int fieldNo_l = 0; fieldNo_l < fieldNo_k; ++fieldNo_l) {
		const int sigmaBeta_ind = fieldNo_k*nFluidFields + fieldNo_l;
		
		kappa(cnt, nodeNo) = - div_test<LT>(unitNormal, cnt, nodeNo, grid);
	
		if (/*sqrt(kappa(cnt, nodeNo)*kappa(cnt, nodeNo))>0.5 ||*/ FNorm(cnt, nodeNo) < 1e-3)
		  kappa(cnt, nodeNo) = 0.0;
		
		IFTforceNode.set(0 ,0) += 0.5*sigma[sigmaBeta_ind]*kappa(cnt, nodeNo)*F(cnt, nodeNo); 
		
		cnt++;
	      }
	    }
	    

	    /*
	    kappa(0, nodeNo) = - div_test2<LT>(unitNormal, 0, nodeNo, grid);
	    if (sqrt(kappa(0, nodeNo)*kappa(0, nodeNo))>=0.5 || FNorm(0, nodeNo) < 1e-4)
		  kappa(0, nodeNo) = 0.0;
	    
	    //IFTforceNode.set(0 ,0) += 0.5*sigma[1]*kappa(0, nodeNo)*F(0, nodeNo); 
	    IFTforceNode.set(0 ,0) += 0.5*0.01*kappa(0, nodeNo)*F(0, nodeNo); 
	    */
	    
	    ForceField.set(0, nodeNo)+=IFTforceNode(0, 0)+bodyForce(0, 0);
	    /*
            std::valarray<lbBase_t> feqTotRel0Node(LT::nQ);
            for (int q=0; q < LT::nQNonZero_; ++q) {
                feqTotRel0Node[q]= wAll[q]* cgat.GammaNonZeroTotNode_;
            }
            feqTotRel0Node[LT::nQNonZero_] = wAll[LT::nQNonZero_]*cgat.Gamma0TotNode_;
	    */
            
            // Calculate local fluid tau
	    lbBase_t visc_inv = rhoRelNode(0, 0)*kin_visc_inv[0];
            for (int fieldNo=1; fieldNo<nFluidFields; ++fieldNo){
		visc_inv += rhoRelNode(0, fieldNo)*kin_visc_inv[fieldNo];
	    }
	    const lbBase_t tauFlNode = LT::c2Inv/visc_inv + 0.5;
	    
	    // Calculate velocity
            // Copy of local velocity distribution
	    auto velNode = calcVel<LT>(fTot(0, nodeNo), rhoTot(0, nodeNo), ForceField(0, nodeNo));
	    /*
	    for (int fieldNo=1; fieldNo<nFluidFields; ++fieldNo){
	      velNode += calcVel<LT>(f(fieldNo, nodeNo), rhoTot(0, nodeNo));
	    }
	    */
	    
           
	    
	    
            // Save density and velocity for printing
            vel.set(0, nodeNo) = velNode;
                
            // BGK-collision term
            const auto u2 = LT::dot(velNode, velNode);
            const auto cu = LT::cDotAll(velNode);
            // Calculate the Guo-force correction
            const auto uF = LT::dot(velNode, ForceField(0, nodeNo));
            const auto cF = LT::cDotAll(ForceField(0, nodeNo));

	    /*
            LbField<LT> fTotNode(1,1);
            fTotNode.set(0,0) = 0;
            */
	    LbField<LT> deltaOmegaRC(1, nFluidFields);

	    LbField<LT> deltaOmegaST(1,1);
	    deltaOmegaST.set(0 ,0) = 0;
	    
	    const auto fTotNode = fTot(0, nodeNo);
	    const auto rhoTotNode = rhoTot(0, nodeNo);
	    const auto fToteqNode = calcfeq<LT>(rhoTotNode, u2, cu);
	    const auto omegaBGKTot = calcOmegaBGK_TEST<LT>(fTotNode, fToteqNode, tauFlNode);
	    
	    const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tauFlNode, cu, uF, cF);

	    
	    fTot.set(0, nodeNo) = fTotNode + deltaOmegaF + omegaBGKTot;
	    
            for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {	        		
	        const auto fNode = f(fieldNo, nodeNo);
	        const auto rhoNode = rho(fieldNo, nodeNo);
              
	        const auto feqNode = calcfeq<LT>(rhoNode, u2, cu);
	        const auto omegaBGK = calcOmegaBGK_TEST<LT>(fNode, feqNode, 1/*tauFlNode*/);

	        

		LbField<LT> deltaOmegaFDiff(1,1);
		deltaOmegaFDiff.set(0 ,0) = calcDeltaOmegaFDiff<LT>(1, rhoRel(fieldNo, nodeNo), cu, uF, cF);  // LBcollision

		
                // Recoloring step
                int field_k_ind = (fieldNo*(fieldNo-1))/2;
                deltaOmegaRC.set(0, fieldNo) = 0;
               
                for (int field_l = 0; field_l < fieldNo; ++field_l) {
                    const int F_ind = field_k_ind + field_l;
		    const int sigmaBeta_ind = fieldNo*nFluidFields + field_l;
                    //deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, sigma[sigmaBeta_ind], FNorm(F_ind, nodeNo), cgat.cDotFRC_(0, F_ind)/cgat.FNorm_(0, F_ind));
		    //deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, sigma[sigmaBeta_ind], FNorm(F_ind, nodeNo), LT::cDotAll(F(F_ind, nodeNo))/(FNorm(F_ind, nodeNo)+(FNorm(F_ind, nodeNo)<lbBaseEps)));
                    deltaOmegaRC.set(0, fieldNo) += beta[sigmaBeta_ind]*rhoRel(field_l, nodeNo)*LT::cDotAll(F(F_ind, nodeNo))*cNormInv/(FNorm(F_ind, nodeNo)+(FNorm(F_ind, nodeNo)<lbBaseEps));
		    //deltaOmegaRC.set(0, fieldNo) += calcDeltaOmegaRC<LT>(beta[sigmaBeta_ind], 1, rhoRel(field_l, nodeNo), 1, LT::cDotAll(F(F_ind, nodeNo)/(FNorm(F_ind, nodeNo)+(FNorm(F_ind, nodeNo)<lbBaseEps)));
		    
                }
                for (int field_l = fieldNo + 1; field_l < nFluidFields; ++field_l) {
                    const int field_k_ind = (field_l*(field_l-1))/2;
                    const int F_ind =  field_k_ind + fieldNo;
		    const int sigmaBeta_ind = fieldNo*nFluidFields + field_l;
                    //deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, sigma[sigmaBeta_ind], FNorm(F_ind, nodeNo), -cgat.cDotFRC_(0, F_ind)/cgat.FNorm_(0, F_ind));
		    //deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, sigma[sigmaBeta_ind], FNorm(F_ind, nodeNo), -LT::cDotAll(F(F_ind, nodeNo))/(FNorm(F_ind, nodeNo)+(FNorm(F_ind, nodeNo)<lbBaseEps)));
                    deltaOmegaRC.set(0, fieldNo) -= beta[sigmaBeta_ind]*rhoRel(field_l, nodeNo)*LT::cDotAll(F(F_ind, nodeNo))*cNormInv/(FNorm(F_ind, nodeNo)+(FNorm(F_ind, nodeNo)<lbBaseEps));
                }

                deltaOmegaRC.set(0, fieldNo) *= wAll*rhoNode;
                //deltaOmegaRC.set(0, fieldNo) *= rhoNode*feqTotRel0Node;

		//fTot.set(0, nodeNo) += deltaOmegaST(0, 0);
            
		

		
                // Calculate lb field
                f.set(fieldNo, nodeNo) = fNode + omegaBGK + deltaOmegaRC(0, fieldNo) + deltaOmegaFDiff(0 ,0)/*+ deltaOmegaST(0, 0)*/;
		// Collision and propagation
		fTmp.propagateTo(fieldNo, nodeNo, f(fieldNo, nodeNo), grid);
		
		//fTot.set(0, nodeNo) +=  deltaOmegaST(0, 0);
            }
	    // Collision and propagation
	    /*
	    for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {
	      fTmp.propagateTo(fieldNo, nodeNo, rhoRel(fieldNo, nodeNo)*fTotNode(0, 0) + omegaRC(0, fieldNo), grid);
            }
	    */
	    
	    // Collision and propagation
	    //fTot.set(0, nodeNo) += deltaOmegaF + omegaBGKTot;
	    fTotTmp.propagateTo(0, nodeNo, fTot(0, nodeNo) + deltaOmegaST(0, 0), grid);
	    
            
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // DIFFUSION
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	   
	    
	    
	    //lbBase_t tauDiff_aveNode = LT::c2Inv/diff_ave_inv + 0.5;
	    for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	      
	      lbBase_t diff_ave_inv = 0.0;
	      for (int fluidNo=0; fluidNo<nFluidFields; ++fluidNo){
		diff_ave_inv += rhoRel(fluidNo, nodeNo)*diff_coef_inv[fieldNo*nFluidFields+fluidNo];
	      }
	      lbBase_t tauDiff_aveNode = LT::c2Inv/diff_ave_inv + 0.5;

	      const auto gNode = g(fieldNo, nodeNo);
	      const lbBase_t rhoDNode = calcRho<LT>(gNode);
	      // Save density and velocity for printing
	      rhoD(fieldNo, nodeNo) = rhoDNode;

	      const auto geqNode = calcfeq<LT>(rhoDNode, u2, cu);
	      const auto omegaBGK = calcOmegaBGK_TEST<LT>(gNode, geqNode, tauDiff_aveNode);
	      g.set(fieldNo, nodeNo) = gNode + omegaBGK;
	    }

	    //------------------------------------------------------
	    //Interaction between diffusive fields and fluid fields
	    //------------------------------------------------------
	    
	    LbField<LT> deltaOmegaDI(1, nDiffFields);
	    for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	      deltaOmegaDI.set(0, fieldNo) = 0;
	    }

	    //lbBase_t W, W_1, W_2, W1, W2;
	    int diffPhaseInd;
	    int solventPhaseInd;
	    LbField<LT> cosPhiTmp(1, 1);
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    //Stored interface normals are the lower triangular part of the interface normal matrix
	    //and point from phase of lower phase index toward phase of higher phase index. e.g., 0->1, 0->2, 1->2 etc.
	    //Instead of changing sign on omegaDI contribution, change sign on potential W to obtain wanted interaction result.
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    /*
	    diffPhaseInd = 10;
	    //soluble in phase 0
	    W = rhoRel(0, nodeNo) - 1;   
	    cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1));
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);   
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 11;
	    //not soluble in phase 1
	    solventPhaseInd = 0;
	    W = -rhoRel(1, nodeNo);
	    cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) - cgat.FNorm_(0,2)*cgat.cosPhi_(0, 2))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,2));
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);   
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 12;
	    //soluble in phase 0
	    solventPhaseInd = 0;
	    W = rhoRel(0, nodeNo) - 1;   
	    cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1));
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);   
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    
	    diffPhaseInd = 0;
	    //not soluble in phase 1
	    solventPhaseInd = 0;
	    W = -rhoRel(1, nodeNo);
	    cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) - cgat.FNorm_(0,2)*cgat.cosPhi_(0, 2))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,2));
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);   
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 1;
	    solventPhaseInd = 1;
	    //surfactant 0-1-interfaces, while soluble in phase 2 
	    W1 = rhoRel(0, nodeNo) - rhoRel(1, nodeNo);
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 0]*W1*cgat.cosPhi_(0, 0);
	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 2;
	    solventPhaseInd = 1;
	    //surfactant 0-1-interfaces
	    W1 = rhoRel(0, nodeNo); // At interface 0-1, not soluble in phase 0 (positive sign since phase 0 is lowest phase in interface) 
	    W1+= -rhoRel(1, nodeNo); // At interface 0-1, not soluble in phase 1
	    W2 = -rhoRel(2, nodeNo); // At interfaces 0-2 and 1-2, not soluble in phase 2 
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 0]*W1*cgat.cosPhi_(0, 0);
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 2]*W2*cgat.cosPhi_(0, 2);
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 2]*W2*cgat.cosPhi_(0, 1);
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 3;
	    //surfactant 1-2-interfaces
	    W1 = rhoRel(1, nodeNo) - rhoRel(2, nodeNo);
	    W2 = rhoRel(0, nodeNo);
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 1*nFluidFields + 2]*W1*cgat.cosPhi_(0, 2);
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 1]*W2*cgat.cosPhi_(0, 0);
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 2]*W2*cgat.cosPhi_(0, 1);
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 4;
	    //surfactant 0-2-interfaces
	    W1 = rhoRel(0, nodeNo); // At interface 0-2, not soluble in phase 0 (positive sign since phase 0 is lowest phase in interface) 
	    W1+= - rhoRel(2, nodeNo); // At interface 0-2, not soluble in phase 2
	    W2 = - rhoRel(1, nodeNo); // At interface 0-1, not soluble in phase 1
	    
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 2]*W1*cgat.cosPhi_(0, 1);
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 1]*W2*cgat.cosPhi_(0, 0);
	    deltaOmegaDI.set(0, diffPhaseInd) -= betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 1*nFluidFields + 2]*W2*cgat.cosPhi_(0, 2);
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 5;
	    solventPhaseInd = 1;
	    //soluble in phase 0 & 1 (or, i.e., not soluble phase 2) 
	    //W = rhoRel(0, nodeNo) - 1;
	   
	    W_1 = -rhoRel(2, nodeNo); // At interface 0-2, not soluble in phase 2
	    W_2 = -rhoRel(2, nodeNo); // At interface 1-2, not soluble in phase 2
	    //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 1]*W_1*cgat.cosPhi_(0, 0);
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 2]*W_1*cgat.cosPhi_(0, 1);
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 1*nFluidFields + 2]*W_2*cgat.cosPhi_(0, 2);
	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 6;
	    solventPhaseInd = 0;
	    //not soluble in phase 0 
	    W = -rhoRel(1, nodeNo);
	    cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) - cgat.FNorm_(0,2)*cgat.cosPhi_(0, 2))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,2));
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);
	    //W_1 = -rhoRel(1, nodeNo);
	    //W_2 = -rhoRel(2, nodeNo);
	    //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cgat.cosPhi_(0, 0);
	    //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 1);	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 7;
	    solventPhaseInd = 0;
	    //W = -(rhoRel(0, nodeNo)-1);
	    //surfactant midpoint phase 0 interfaces
	    W = rhoRel(0, nodeNo) - 0.5;
	    cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1));
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);
	    //W_1 = -rhoRel(1, nodeNo);
	    //W_2 = -rhoRel(2, nodeNo);
	    //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cgat.cosPhi_(0, 0);
	    //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 1);	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 8;
	    solventPhaseInd = 0;
	    //surfactant on phase 0 side of interfaces
	    W = rhoRel(0, nodeNo) - 0.75;
	    cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1));
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);
	    //W_1 = -rhoRel(1, nodeNo);
	    //W_2 = -rhoRel(2, nodeNo);
	    //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cgat.cosPhi_(0, 0);
	    //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 1);	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 9;
	    solventPhaseInd = 0;
	    //surfactant on phase 0 side of interfaces
	    W = rhoRel(0, nodeNo) - 0.75;
	    cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1));
	    deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);
	    //W_1 = -rhoRel(1, nodeNo);
	    //W_2 = -rhoRel(2, nodeNo);
	    //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cgat.cosPhi_(0, 0);
	    //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 1);	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    
	    */	    
	    //------------------------------------------------------
	    //END Interaction between diffusive fields and fluid fields
	    //------------------------------------------------------
	    
	    
	    for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	      deltaOmegaDI.set(0, fieldNo) *= wAll*rhoD(fieldNo, nodeNo);
	      gTmp.propagateTo(fieldNo, nodeNo, g(fieldNo, nodeNo) + deltaOmegaDI(0, fieldNo) , grid);
	    }
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // END DIFFUSION
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	    

        } // End nodes
        // Swap data_ from fTmp to f;
	fTot.swapData(fTotTmp);  // LBfield
        f.swapData(fTmp);  // LBfield
	g.swapData(gTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
	mpiBoundary.communicateLbField(fTot, grid);
        // Half way bounce back
        bounceBackBnd.apply(fTot, grid);

	
        mpiBoundary.communicateLbField(f, grid);
        // Half way bounce back
        bounceBackBnd.apply(f, grid);

	mpiBoundary.communicateLbField(g, grid);
        // Half way bounce back
        bounceBackBnd.apply(g, grid);

        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0)  ) {
	  /*
	    for (auto nn: bulkNodes) {
                velIO(0, 0, nn) = vel(0, 0, nn);
                velIO(0, 1, nn) = vel(0, 1, nn);
                velIO(0, 2, nn) = 0; 

            }
            output.write("lb_run", i);
	  */
	    output.write(i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
            }
        }

    } // End iterations

    MPI_Finalize();

    return 0;
}
