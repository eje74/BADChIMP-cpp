// //////////////////////////////////////////////
//
// BADChIMP partial_mix2
//
// For documentation see:
//    doc/documentation.pdf
// 
// //////////////////////////////////////////////

#include "../LBSOLVER.h"
#include "../IO.h"
#include "./LBpartial_mixHelp.h"

//                                SET THE LATTICE TYPE
//------------------------------------------------------------------------------------- SET THE LATTICE TYPE
#define LT D2Q9

#define VTK_CELL VTK::pixel
//#define LT D3Q19
//#define VTK_CELL VTK::voxel

int main()
{
  //=====================================================================================
  //
  //                                      SETUP MPI
  //
  //=====================================================================================
  MPI_Init(NULL, NULL);
  int nProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  
  
  //=====================================================================================
  //
  //                         SETUP THE INPUT AND OUTPUT PATHS
  //
  //=====================================================================================
  std::string chimpDir = "./";
  // std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
  std::string mpiDir = chimpDir + "input/mpi/";
  std::string inputDir = chimpDir + "input/";
  std::string outputDir = chimpDir + "output/";
  
  
  //=====================================================================================
  //
  //                               SETUP GRID AND GEOMETRY
  //
  //=====================================================================================
  Input input(inputDir + "input.dat");
  LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
  Grid<LT> grid(vtklb);
  Nodes<LT> nodes(vtklb, grid);
  BndMpi<LT> mpiBoundary(vtklb, nodes, grid);
  // Set bulk nodes
  std::vector<int> bulkNodes = findBulkNodes(nodes);
  // Set solid boundary 
  std::vector<int> solidBoundaryNodes = findSolidBndNodes(nodes);
  

  //=====================================================================================
  //
  //                                  READ FROM INPUT
  //
  //=====================================================================================

  //                                Number of iterations
  //------------------------------------------------------------------------------------- Number of iterations
  int nIterations = input["iterations"]["max"];
  
  //                                  Write interval
  //------------------------------------------------------------------------------------- Write interval  
  int nItrWrite = input["iterations"]["write"];
  
  //                                    Fluid Flow
  //===================================================================================== Fluid Flow

  //                              Number of fluid fields
  //------------------------------------------------------------------------------------- Number of fluid fields
  const int nFluidFields = input["fluid"]["numfluids"];
 
  
  //                               Kinematic viscosities
  //------------------------------------------------------------------------------------- Kinematic viscosities
  auto kin_visc = input["fluid"]["viscosity"];
  std::valarray<lbBase_t> kin_visc_inv(nFluidFields);
  for (int n = 0; n < nFluidFields; ++n){
    kin_visc_inv[n] = 1./kin_visc[n];
  }

  //                                     Body force
  //------------------------------------------------------------------------------------- Body force
  VectorField<LT> bodyForce(1, 1, input["fluid"]["bodyforce"]);
  
  //                          Parameters for compressibility
  //------------------------------------------------------------------------------------- Parameters for compressibility
  std::valarray<lbBase_t> alpha = LT::w0*inputAsValarray<lbBase_t>(input["fluid"]["alpha"]);
  std::valarray<lbBase_t> Gamma0 = inputAsValarray<lbBase_t>(input["fluid"]["alpha"]);
  std::valarray<lbBase_t> GammaNonZero = (1-LT::w0*Gamma0)/(1-LT::w0);
  
  //                                 Phase separation: fluids
  //------------------------------------------------------------------------------------- Phase separation: fluids
  std::valarray<lbBase_t> beta = inputAsValarray<lbBase_t>(input["fluid"]["phaseseparation"]["beta"]);
  

  
  //                                        Surface tension
  //------------------------------------------------------------------------------------- Surface tension
  std::valarray<lbBase_t> sigma = inputAsValarray<lbBase_t>(input["fluid"]["phaseseparation"]["sigma"]);
  

  
  //                                           Diffusion
  //===================================================================================== Diffusion

  // Number of diffusive fields
  //------------------------------------------------------------------------------------- Number of diffusive fields
  const int nDiffFields = input["diffusion"]["numfields"];
  
  // Diffusion coefficients
  //------------------------------------------------------------------------------------- Diffusion coefficients
  std::valarray<lbBase_t> diff_coef = inputAsValarray<lbBase_t>(input["diffusion"]["fluidinteraction"]["diffcoef"]);
  std::valarray<lbBase_t> diff_coef_inv(nDiffFields*nFluidFields);
  for (int i=0; i < nDiffFields; ++i) {
    for (int j=0; j < nFluidFields; ++j) {
      diff_coef_inv[i*nFluidFields + j] = 1./diff_coef[i*nFluidFields + j];
    }
  }

  //                                 Phase separation: Diffusive fields
  //------------------------------------------------------------------------------------- Phase separation: Diffusive fields
  std::valarray<lbBase_t> betaDiff = inputAsValarray<lbBase_t>(input["diffusion"]["fluidinteraction"]["beta"]);
  
  
  //                                      Partial Miscibility
  //===================================================================================== Partial Miscibility

  //                                       Henry's constant H
  //------------------------------------------------------------------------------------- Henry's constant H
  const lbBase_t H = input["Partial-Misc"]["H"];

  //                                   Fixed Concentration Factor
  //------------------------------------------------------------------------------------- Fixed Concentration Factor
  const lbBase_t FixedConcFactor = input["Partial-Misc"]["FixedConcFactor"];

  //                                   Initial Fixed Salt Concentration 
  //------------------------------------------------------------------------------------- Initial Fixed Salt Concentration 
  const lbBase_t InitialSaltConc = input["Partial-Misc"]["InitialSaltConc"];

  //                                   Kinetic Reaction Constant 
  //------------------------------------------------------------------------------------- Kinetic Reaction Constant
  const lbBase_t kinConst = input["Partial-Misc"]["kinConst"];
  
  //                                    Output directory number
  //------------------------------------------------------------------------------------- Output directory number
  std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));  
  std::string outputDir2 = outputDir + "/out" + dirNum;
  
  
  //                                Recoloration Constant 2k
  //------------------------------------------------------------------------------------- Recoloration Constant 2k
  lbBase_t kx2_test = 0;
  for (int q = 0; q < LT::nQNonZero_; ++q) {
    kx2_test += LT::w[q] * LT::c(q,0)*LT::c(q,0) /  LT::cNorm[q];
  }
  //------------------------------------------------------------------------------------- END READ FROM INPUT



  
  //=====================================================================================
  //
  //                               WRITE PARAMETERS TO SCREEN
  //
  //=====================================================================================

  if (myRank==0) {

    //                              Number of fluid fields
    //------------------------------------------------------------------------------------- Number of fluid fields
    std::cout << "nFluidFields = " << nFluidFields << std::endl;
    std::cout << std::endl;

    // kinematic viscosity
    //------------------------------------------------------------------------------------- kinematic viscosity
    std::cout << "kin. visc = " << std::endl;
    for (int n = 0; n < nFluidFields; ++n){
      std::cout << " " <<  kin_visc[n];
    }
    std::cout << std::endl;

    
    //                          Parameters for compressibility
    //------------------------------------------------------------------------------------- Parameters for compressibility
    std::cout << "alpha ="<< std::endl;
    for (int n = 0; n < nFluidFields; ++n)
      std::cout << " " << alpha[n];
    std::cout << std::endl;

    //                                 Phase separation: fluids
    //------------------------------------------------------------------------------------- Phase separation: fluids
    std::cout << "beta = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
      for (int j=0; j < nFluidFields; ++j) {
	std::cout << " " << beta[i*nFluidFields + j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    //                                     Surface tension
    //------------------------------------------------------------------------------------- Surface tension
    std::cout << "sigma = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
      for (int j=0; j < nFluidFields; ++j) {
	std::cout << " " << sigma[i*nFluidFields + j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    // Number of diffusive fields
    //------------------------------------------------------------------------------------- Number of diffusive fields
    std::cout << "nDiffFields = " << nDiffFields << std::endl;

    // Diffusion coefficients
    //------------------------------------------------------------------------------------- Diffusion coefficients
    std::cout << "Diff. coef = " << std::endl;
    for (int i=0; i < nDiffFields; ++i) {
      std::cout << "diff.field[" << i <<"]: ";
      for (int j=0; j < nFluidFields; ++j) {
	std::cout << " " << diff_coef[i*nFluidFields + j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    //                                 Phase separation: Diffusive fields
    //------------------------------------------------------------------------------------- Phase separation: Diffusive fields
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

    //                                Recoloration Constant 2k
    //------------------------------------------------------------------------------------- Recoloration Constant 2k
    std::cout<<"Recoloration Constant 2k = "<< (kx2_test) <<std::endl;
  }
  //------------------------------------------------------------------------------------- END WRITE PARAMETERS TO SCREEN


  
  //=====================================================================================
  //
  //                                  DEFINE RHEOLOGY
  //
  //=====================================================================================
  Newtonian<LT> newtonian(kin_visc_inv[0]*LT::c2Inv + 0.5);

    
  //=====================================================================================
  //
  //                                MACROSCOPIC FIELDS
  //
  //=====================================================================================

  //                                     Fluid Flow
  //===================================================================================== Fluid Flow
  //                                      Density
  //------------------------------------------------------------------------------------- Density
  ScalarField rho(nFluidFields, grid.size());
  ScalarField rhoRel(nFluidFields, grid.size());
  ScalarField phi(nFluidFields, grid.size());
  ScalarField rhoTot(1, grid.size());
  
  //                                      Fracture Geometry
  //------------------------------------------------------------------------------------- Fracture Geometry
  ScalarField height(1, grid.size());
  ScalarField tmpGradHeight(2, grid.size());
  VectorField<LT> gradHeight(1, grid.size());

  //                                      Sources
  //------------------------------------------------------------------------------------- Sources
  ScalarField Rfield(nFluidFields, grid.size());
  ScalarField Q(1, grid.size());
  
  //                              Initiate density from file
  //------------------------------------------------------------------------------------- Initiate density from file
  setScalarAttribute(rho, "init_rho_", vtklb);
  
  /*
    for (auto nodeNo: bulkNodes) {
    rho(0, nodeNo) += rho(1, nodeNo);
    rho(1, nodeNo) = 0.0;
    
    }
  */ 

  
  
  //                              Initiate frac geo from file
  //------------------------------------------------------------------------------------- Initiate frac geo from file
  setScalarAttribute(height, "frac_height_", vtklb);
  setScalarAttribute(tmpGradHeight, "gradH_", vtklb);
  for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {	
    gradHeight(0, 0, nodeNo) = tmpGradHeight(0, nodeNo);
    gradHeight(0, 1, nodeNo) = tmpGradHeight(1, nodeNo);
  }

  VectorField<LT> ForceField2D(5, grid.size());
  
    
  //                                    Wall wettability
  //------------------------------------------------------------------------------------- Wall wettability
  setScalarAttributeWall(rhoRel, "init_rho_wall_", vtklb, nodes);
  
  //                                   Scale wettability
  //------------------------------------------------------------------------------------- scale wettability
  normelizeScalarField(rhoRel, solidBoundaryNodes);
  
    
    
    
  //                                       Velocity
  //------------------------------------------------------------------------------------- Velocity
  VectorField<LT> vel(1, grid.size());
  //                                   Initiate velocity
  //------------------------------------------------------------------------------------- Initiate velocity
  for (auto nodeNo: bulkNodes) {
    vel.set(0, nodeNo) = 0;
  }
    

  //                                    Phase separation
  //                             Color gradients for fluid pair
  //                     Absolute values of color gradients for fluid pair
  //------------------------------------------------------------------------------------- Phase separation
  ScalarField FNorm(nFluidFields*(nFluidFields-1)/2, grid.size());
  VectorField<LT> F(nFluidFields*(nFluidFields-1)/2, grid.size());
  VectorField<LT> ForceField(1, grid.size());
  VectorField<LT> J_RC_prev(nFluidFields, grid.size());
  VectorField<LT> J_DI_prev(nDiffFields, grid.size());
  
  VectorField<LT> unitNormal(nFluidFields*(nFluidFields-1)/2, grid.size());


  
  //                         Curvature of interface for fluid pair
  //                                 (lower triangular matrix)
  //------------------------------------------------------------------------------------- Curvature of interface for fluid pair
  ScalarField kappa(nFluidFields*(nFluidFields-1)/2, grid.size());
  ScalarField kappa2(nFluidFields*(nFluidFields-1)/2, grid.size());

  ScalarField cosAng(nFluidFields*(nFluidFields-1)/2, grid.size());
  
  ScalarField TCap(LT::nD*(LT::nD+1)/2, grid.size());
  
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
  for (auto fieldNo=0; fieldNo < J_DI_prev.num_fields(); ++fieldNo) {
    for (auto nodeNo: bulkNodes) {
      J_DI_prev.set(fieldNo, nodeNo) = 0;
    }
  }
  for (auto fieldNo=0; fieldNo < J_RC_prev.num_fields(); ++fieldNo) {
    for (auto nodeNo: bulkNodes) {
      J_RC_prev.set(fieldNo, nodeNo) = 0;
    }
  }
  
  
  //                                   Advection-Diffusion 
  //===================================================================================== Advection-Diffusion 
  //                                  Density & Concentration
  //------------------------------------------------------------------------------------- Density & Concentration
  ScalarField rhoD(nDiffFields, grid.size());
  ScalarField phiD(nDiffFields, grid.size());

  // Initiate diffusive fields
  //------------------------------------------------------------------------------------- 
  
  //                                Initiate density from file
  //------------------------------------------------------------------------------------- Initiate density from file
  setScalarAttribute(rhoD, "init_rhoD_", vtklb);
  
  
    for (auto nodeNo: bulkNodes) {
      /*for (int fieldNo=0; fieldNo < rhoD.num_fields(); ++fieldNo) {
	rhoD(fieldNo, nodeNo) = 0.01*rho(0, nodeNo);
	}*/
      //rhoD(1, nodeNo) = 0.0*rho(0, nodeNo);
      /*if(grid.pos(nodeNo, 0) > 0.4*(vtklb.getGlobaDimensions(0)-2))
	rhoD(1, nodeNo) = 0.0*rho(0, nodeNo);
      */
    }
      

  
  //=====================================================================================
  //
  //                                SETUP BOUNDARY
  //
  //=====================================================================================

  HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

  //=====================================================================================
  //
  //                                  LB FIELDS
  //
  //=====================================================================================

  //                                Fluid Fields
  //===================================================================================== Fluid Fields
  //------------------------------------------------------------------------------------- 
  LbField<LT> f(nFluidFields, grid.size()); 
  LbField<LT> fTmp(nFluidFields, grid.size());

  //                           Initiate lb distributions
  //------------------------------------------------------------------------------------- Initiate lb distributions
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

  //              Initiate Total density, concentrations and lb distributions
  //------------------------------------------------------------------------------------- Initiate Total density, concentrations and lb distributions
  
  for (auto nodeNo: bulkNodes) {
    rhoTot(0, nodeNo) = 0.0;
    for (int fieldNo=0; fieldNo < rho.num_fields(); ++fieldNo) {
      rhoTot(0, nodeNo) += rho(fieldNo, nodeNo);
    }
    for (int fieldNo=0; fieldNo < rho.num_fields(); ++fieldNo) {
      phi(0, nodeNo) = rho(fieldNo, nodeNo)/rhoTot(0, nodeNo);
    }
    
    for (int fieldNo=0; fieldNo < rhoD.num_fields(); ++fieldNo) {
      phiD(fieldNo, nodeNo) = rhoD(fieldNo, nodeNo)/rhoTot(0, nodeNo);
    }
    
    for (int q = 0; q < LT::nQ; ++q) {
      fTot(0, q, nodeNo) = LT::w[q]*rhoTot(0, nodeNo);
      fTotTmp(0, q, nodeNo) = 0;
    }
  }
    

    
  //                                Diffusive Fields
  //===================================================================================== Diffusive Fields

  //------------------------------------------------------------------------------------- 
  LbField<LT> g(nDiffFields, grid.size()); 
  LbField<LT> gTmp(nDiffFields, grid.size());

  //                            Initiate lb distributions
  //------------------------------------------------------------------------------------- Initiate lb distributions
  for (int fieldNo=0; fieldNo < g.num_fields(); ++fieldNo) {
    for (auto nodeNo: bulkNodes) {
      for (int q = 0; q < LT::nQ; ++q) {	      
	g(fieldNo, q, nodeNo) = LT::w[q]*rhoD(fieldNo, nodeNo);
	gTmp(fieldNo, q, nodeNo) = 0;
      }
    }
  }
    
    
  //=====================================================================================
  //
  //                                  OUTPUT VTK
  //
  //=====================================================================================
  Output<LT> output(grid, bulkNodes, outputDir2, myRank, nProcs); 
  output.add_file("lb_run");
  output.add_scalar_variables({"rhoTot", "rho", "rhoD", "kappa_", "kappa2_", "R",     "Q", "frac_height", "grad_height"}, 
			      { rhoTot,   rho,   rhoD,   kappa,    kappa2,    Rfield,  Q,   height,   tmpGradHeight});
  output.add_vector_variables({"vel", "F", "unitNormal", "forceField", "gradHeight", "forceField2D_"}, 
			      { vel,   F,   unitNormal,   ForceField,   gradHeight, ForceField2D});


  

  //        Just ad hoc helper fields (Should be set outside the loop structure)
  //------------------------------------------------------------------------------------- ad hoc helper fields (Should be set outside the loop structure)
  std::valarray<lbBase_t> cNormInv(LT::nQ);
  std::valarray<lbBase_t> wAll(LT::nQ);
  std::valarray<lbBase_t> cNormTmp(LT::nQ);
  for (int q=0; q < LT::nQNonZero_; ++q) {
    cNormInv[q] = 1.0/LT::cNorm[q];
    wAll[q] = LT::w[q];
    cNormTmp[q] = LT::cNorm[q];
  }
  cNormInv[LT::nQNonZero_] = 0;
  cNormTmp[LT::nQNonZero_] = 0;
  wAll[LT::nQNonZero_] = LT::w[LT::nQNonZero_];
  
  //=====================================================================================
  //
  //                                  MAIN LOOP
  //
  //=====================================================================================

  for (int i = 0; i <= nIterations; i++) {
    //           Macroscopic values : rho, rhoTot and rhoRel = rho/rhoTot
    //------------------------------------------------------------------------------------- Macroscopic values : rho, rhoTot and rhoRel = rho/rhoTot
    calcDensityFields(rho, rhoRel, rhoTot, rhoD, phi, phiD, bulkNodes, f, fTot, g);
    
    
    mpiBoundary.communciateScalarField(rhoRel);
    mpiBoundary.communciateScalarField(phi);
    mpiBoundary.communciateScalarField(phiD);
    
    for (auto nodeNo: bulkNodes) {
      // Cacluate gradient
      //------------------------------------------------------------------------------------- 
      ScalarField rhoRelNode(1, nFluidFields);
      for (int fieldNo = 0; fieldNo < nFluidFields; ++fieldNo){
	rhoRelNode(0, fieldNo) = rhoRel(fieldNo, nodeNo);
      }
      const CGAttributes<LT> cgat(nFluidFields, nodeNo, cNormInv, Gamma0, rhoRelNode, rhoRel, grid);

      //-----------------------------------------------------------------------------------------------------------------------------------
      //Stored interface normals are the lower triangular part of the interface normal matrix
      //and point from phase of lower phase index toward phase of higher phase index. e.g., 0->1, 0->2, 1->2 etc.
      //-----------------------------------------------------------------------------------------------------------------------------------
      
      for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	FNorm(cnt, nodeNo)= cgat.FNorm_(0, cnt);
	F.set(cnt, nodeNo)= cgat.F_(0, cnt);
	unitNormal.set(cnt, nodeNo)= F(cnt, nodeNo)/(FNorm(cnt, nodeNo)+(FNorm(cnt, nodeNo)<lbBaseEps));
      }

      cosAng(0, nodeNo) = -1;
      cosAng(1, nodeNo) = 0.0;
      cosAng(2, nodeNo) = 1;
      
      
    }
    
    mpiBoundary.communciateVectorField_TEST(unitNormal);
    
    mpiBoundary.communciateVectorField_TEST(F);
    mpiBoundary.communciateScalarField(FNorm);
    
    /*
    //Cap Tensor Calc
    for (auto nodeNo: bulkNodes) {
    
    
    // Recoloring step
    int field_k_ind = (fieldNo*(fieldNo-1))/2;
    LbField<LT> deltaOmegaST(1,1);
    deltaOmegaST.set(0 ,0) = 0;
    for (int field_l = 0; field_l < fieldNo; ++field_l) {
    
    const int F_ind = field_k_ind + field_l;
    const int sigmaBeta_ind = fieldNo*nFluidFields + field_l;
    const auto cn = LT::cDotAll(F(F_ind, nodeNo)/(FNorm(F_ind, nodeNo)+(FNorm(F_ind, nodeNo)<lbBaseEps)));
    
    deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, 2*sigma[sigmaBeta_ind]*rhoTot(0, nodeNo), FNorm(F_ind, nodeNo), cn);
    }
    }
    */
    
    //Main calculation loop
    //------------------------------------------------------------------------------------- 
    for (auto nodeNo: bulkNodes) {
      
      Q(0, nodeNo)=0.0;


      
      
	
      
      
      // Cacluate gradient
      //------------------------------------------------------------------------------------- 
      ScalarField rhoRelNode(1, nFluidFields);
      for (int fieldNo = 0; fieldNo < nFluidFields; ++fieldNo){
	Rfield(fieldNo, nodeNo)=0.0;
	rhoRelNode(0, fieldNo) = rhoRel(fieldNo, nodeNo);
      }

      
      if(grid.pos(nodeNo, 0)>=2 && grid.pos(nodeNo, 0)<=4 && grid.pos(nodeNo, 1)>=2 && grid.pos(nodeNo, 1)<= vtklb.getGlobaDimensions(1) -5 ){
	lbBase_t fixedInletDens = 1.0; 
	Q(0, nodeNo) += 1e-3;//2*(fixedInletDens - rhoTot(0, nodeNo)); //1e-3;
	Rfield(0, nodeNo)+= 2*(rhoTot(0, nodeNo) +0.5*Q(0, nodeNo) - rho(0, nodeNo));
	Rfield(1, nodeNo)+= 2*(0 - rho(1, nodeNo));
	Rfield(2, nodeNo)+= 2*(0 - rho(2, nodeNo));
	

	for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	  FNorm(cnt, nodeNo) = 0;
	  F.set(cnt, nodeNo) = 0*F(cnt, nodeNo);
	  unitNormal.set(cnt, nodeNo)= F(cnt, nodeNo)/(FNorm(cnt, nodeNo)+(FNorm(cnt, nodeNo)<lbBaseEps));
	}
      }

      
      if ( grid.pos(nodeNo, 0)>= vtklb.getGlobaDimensions(0) -4 ){
	
	lbBase_t fixedOutletDens = 1.0; 
	Q(0, nodeNo) += 2*(fixedOutletDens - rhoTot(0, nodeNo));
	
	for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo){
	  Rfield(fieldNo, nodeNo) += Q(0, nodeNo) * rhoRel(fieldNo, nodeNo);
	}
	

	
	
	for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	  FNorm(cnt, nodeNo) = 0;
	  F.set(cnt, nodeNo) = 0*F(cnt, nodeNo);
	  unitNormal.set(cnt, nodeNo)= F(cnt, nodeNo)/(FNorm(cnt, nodeNo)+(FNorm(cnt, nodeNo)<lbBaseEps));
	}
	
	
      }

      
      
      const CGAttributes<LT> cgat(nFluidFields, nodeNo, cNormInv, Gamma0, rhoRelNode, rhoRel, grid);

      /*
      //                      Fixed Pressure Source & Fixed Phase Field 1 Source 
      //------------------------------------------------------------------------------------- Fixed Pressure Source & Fixed Phase Field 1 Source 
      if(rhoRel(1, nodeNo)> 0.99999 && rhoRel(1, nodeNo)< 0.99999999){
	Q(0, nodeNo) = 2*(1.0 - rhoTot(0, nodeNo));
	Rfield(1, nodeNo) = 2*(1.0*rhoRel(1, nodeNo) - rho(1, nodeNo));
	
	rhoTot(0, nodeNo) += 0.5*Q(0, nodeNo);
	rho(1, nodeNo) += 0.5*Rfield(1, nodeNo);
	phi(1, nodeNo) = rho(1, nodeNo)/rhoTot(0, nodeNo);
      }
      */

      
      
      
      
      VectorField<LT> IFTforceNode(1,1);
      IFTforceNode.set(0 ,0) = 0;
      ForceField.set(0, nodeNo)=0.0;
      //ForceField(0, 0, nodeNo)=1e-3;
      for (int forceNo=0; forceNo < ForceField2D.num_fields(); ++forceNo) 
	ForceField2D.set(forceNo, nodeNo)=0.0;
      
      
      //                      Lower triangular traversing of kappa and IFT
      //------------------------------------------------------------------------------------- Lower triangular traversing of kappa and IFT
      int cnt = 0;
      for (int fieldNo_k=0; fieldNo_k<nFluidFields; ++fieldNo_k) {
	for (int fieldNo_l = 0; fieldNo_l < fieldNo_k; ++fieldNo_l) {
	  const int sigmaBeta_ind = fieldNo_k*nFluidFields + fieldNo_l;
	  
	  kappa(cnt, nodeNo) = - div_test2<LT>(unitNormal, cnt, nodeNo, grid);
	  
	  lbBase_t absGradTmp = 0.5*beta[sigmaBeta_ind]*rhoRel(cnt, nodeNo)*(1-rhoRel(cnt, nodeNo));
	  
	  
	  
	  
	  //kappa2(cnt, nodeNo) = -(div_test2<LT>(F, cnt, nodeNo, grid) - LT::dot(grad<LT>(FNorm, cnt, nodeNo, grid),unitNormal(cnt, nodeNo)) );
	  
	  kappa2(cnt, nodeNo) = kappa(cnt, nodeNo);
	  
	  
	  if (FNorm(cnt, nodeNo) < 1e-3)
	    kappa2(cnt, nodeNo) = 0.0;
	  
	  //if (FNorm(cnt, nodeNo) < 1e-4)
	  //  kappa2(cnt, nodeNo) = 0.0;
	  
	  const lbBase_t IFT_threshold = std::min(1e6*rhoRel(fieldNo_k, nodeNo)*rhoRel(fieldNo_l, nodeNo),1.0);
	  
	  if (FNorm(cnt, nodeNo) > 2.5e-3 && (kappa(cnt, nodeNo)*kappa(cnt, nodeNo))< 1.8) 
	    //IFTforceNode.set(0 ,0) += 0.5*/*rhoTot(0, nodeNo)**/sigma[sigmaBeta_ind]*kappa(cnt, nodeNo)*F(cnt, nodeNo)/**IFT_threshold*/;
	    
	  //IFTforceNode.set(0 ,0) += 0.25*4/beta[sigmaBeta_ind]*sigma[sigmaBeta_ind]*kappa(cnt, nodeNo)*FNorm(cnt, nodeNo)*F(cnt, nodeNo);
	  //IFTforceNode.set(0 ,0) += 0.5*sigma[sigmaBeta_ind]*kappa2(cnt, nodeNo)*unitNormal(cnt, nodeNo);
	  //IFTforceNode.set(0 ,0) += 1.5*1/beta[sigmaBeta_ind]*sigma[sigmaBeta_ind]*kappa2(cnt, nodeNo)*F(cnt, nodeNo);
	  //IFTforceNode.set(0 ,0) += 2*1.5*4/beta[sigmaBeta_ind]*sigma[sigmaBeta_ind]*kappa2(cnt, nodeNo)*absGradTmp*absGradTmp*unitNormal(cnt, nodeNo);
	  //IFTforceNode.set(0 ,0) += 2*1.5*4/beta[sigmaBeta_ind]*sigma[sigmaBeta_ind]*kappa2(cnt, nodeNo)*absGradTmp*unitNormal(cnt, nodeNo);

	    //Quasi-2D
	    
	    if ( grid.pos(nodeNo, 0)< vtklb.getGlobaDimensions(0) -4 ){
	      ForceField2D.set(4, nodeNo) = 0.5*sigma[sigmaBeta_ind]*F(cnt, nodeNo)*2*cosAng(cnt, nodeNo)/height(0, nodeNo);

	      IFTforceNode.set(0 ,0) += 0.5*sigma[sigmaBeta_ind]*F(cnt, nodeNo)*2*cosAng(cnt, nodeNo)/height(0, nodeNo);

	      ForceField2D.set(3, nodeNo) = -0.5*sigma[sigmaBeta_ind]*FNorm(cnt, nodeNo)/height(0, nodeNo)
		*(unitNormal(cnt, nodeNo)*LT::dot(unitNormal(cnt, nodeNo),gradHeight(0,nodeNo)) - gradHeight(0,nodeNo)); //Averaged Capillary tensor effect

	      IFTforceNode.set(0 ,0) += -0.5*sigma[sigmaBeta_ind]*FNorm(cnt, nodeNo)/height(0, nodeNo)
		*(unitNormal(cnt, nodeNo)*LT::dot(unitNormal(cnt, nodeNo),gradHeight(0,nodeNo)) - gradHeight(0,nodeNo)); //Averaged Capillary tensor effect
	      
	    }
	  //Quasi-2D
	  
	    
	  cnt++;
	}
      }
	 
	  
      ForceField.set(0, nodeNo)+=IFTforceNode(0, 0)+bodyForce(0, 0);

      

      //                              Calculate local fluid viscosity
      //------------------------------------------------------------------------------------- Calculate local fluid viscosity
      lbBase_t visc_inv = 0.0;
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo){
	visc_inv += rhoRelNode(0, fieldNo)*kin_visc_inv[fieldNo];
      }
      
      lbBase_t viscNode = 1/visc_inv;

      //Quasi-2D
      
      if ( grid.pos(nodeNo, 0)< vtklb.getGlobaDimensions(0) -4 ){
	std::valarray<lbBase_t> velNodeTmp = calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo)), ForceField(0, nodeNo));

	lbBase_t velFactorNode = (1+0.5*12*viscNode/(height(0, nodeNo)*height(0, nodeNo)));

	auto rhoTotNodeTmp = (rhoTot(0, nodeNo)+0.5*Q(0, nodeNo))*1/(1+0.5*LT::dot(LT::qSumC(fTot(0, nodeNo)),gradHeight(0,nodeNo))/(velFactorNode*height(0, nodeNo)));

	velNodeTmp *= 1/(rhoTotNodeTmp*velFactorNode);
      
	ForceField.set(0, nodeNo) += -12*rhoTotNodeTmp*viscNode*velNodeTmp/(height(0, nodeNo)*height(0, nodeNo));

	ForceField2D.set(0, nodeNo) = -12*rhoTotNodeTmp*viscNode*velNodeTmp/(height(0, nodeNo)*height(0, nodeNo));

	//ForceField.set(0, nodeNo) += -(rhoTotNodeTmp-1)*LT::c2*gradHeight(0,nodeNo)/height(0, nodeNo);

	ForceField2D.set(1, nodeNo) = -(rhoTotNodeTmp-1)*LT::c2*gradHeight(0,nodeNo)/height(0, nodeNo);
	
	lbBase_t Q2DNode = -rhoTotNodeTmp*LT::dot(velNodeTmp,gradHeight(0,nodeNo))/height(0, nodeNo);
      
	Q(0, nodeNo) += Q2DNode;

      
	for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo){
	  Rfield(fieldNo, nodeNo) += Q2DNode* LT::qSum(f(fieldNo, nodeNo))/LT::qSum(fTot(0, nodeNo));
	}
      }

      if ( grid.pos(nodeNo, 0)>= vtklb.getGlobaDimensions(0) -4 ){
	auto velNodeTmp = calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo)), ForceField(0, nodeNo));
	lbBase_t velThresh = 9e-2;
	if (sqrt( LT::dot(velNodeTmp, velNodeTmp))>velThresh)
	  ForceField.set(0, nodeNo) = 2*(0 - LT::qSumC(fTot(0, nodeNo)));
      }
      
      //Quasi-2D END

      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo){
	rho(fieldNo, nodeNo) += 0.5*Rfield(fieldNo, nodeNo);
      }
      rhoTot(0, nodeNo) += 0.5*Q(0, nodeNo);

      
      /*
	std::valarray<lbBase_t> feqTotRel0Node(LT::nQ);
	for (int q=0; q < LT::nQNonZero_; ++q) {
	feqTotRel0Node[q]= wAll[q]* cgat.GammaNonZeroTotNode_;
	}
	feqTotRel0Node[LT::nQNonZero_] = wAll[LT::nQNonZero_]*cgat.Gamma0TotNode_;
      */
      
      //                              Calculate local fluid tau
      //------------------------------------------------------------------------------------- Calculate local fluid tau
      

      int fluidIntIndNode = 0;
      //lbBase fluidInt
      //-------------------------------------------------------------------------------------
      // The fluidIntIndNode states how many types of fluid interfaces that is present at present node 
      //------------------------------------------------------------------------------------- 
      
      for(int cnt=0; cnt<FNorm.num_fields(); ++cnt){
	if(FNorm(cnt, nodeNo)>1e-3){
	  fluidIntIndNode++;
	}
      }

      if(fluidIntIndNode!=0)
	visc_inv = 1/0.16667;
      
      const lbBase_t tauFlNode = LT::c2Inv/visc_inv + 0.5;

      
      auto velNodeTmp = calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo)), ForceField(0, nodeNo));
      auto strainRateNode =  1/(2*rhoTot(0,nodeNo)*LT::c2*tauFlNode)*calcStrainRateTilde<LT>(fTot(0, nodeNo), rhoTot(0, nodeNo), velNodeTmp, ForceField(0, nodeNo), Q(0, nodeNo));
      VectorField<LT> vecTmpNode(1,1);
      vecTmpNode.set(0,0)  = 2*viscNode/height(0, nodeNo)*LT::contractionLowTriVec(strainRateNode, gradHeight(0,nodeNo));

      ForceField.set(0, nodeNo) += vecTmpNode(0,0);

      ForceField2D.set(2, nodeNo) = vecTmpNode(0,0);
      
      
      
      // Calculate velocity
      //------------------------------------------------------------------------------------- 
      // Copy of local velocity distribution
      //------------------------------------------------------------------------------------- 
      auto velNode = calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo)), ForceField(0, nodeNo));
      /*
	for (int fieldNo=1; fieldNo<nFluidFields; ++fieldNo){
	velNode += calcVel<LT>(f(fieldNo, nodeNo), rhoTot(0, nodeNo));
	}
      */
      
    
      
	
      // Save density and velocity for printing
      //------------------------------------------------------------------------------------- 
      vel.set(0, nodeNo) = velNode;
      
      // BGK-collision term
      //------------------------------------------------------------------------------------- 
      const auto u2 = LT::dot(velNode, velNode);
      const auto cu = LT::cDotAll(velNode);
      // Calculate the Guo-force correction
      //------------------------------------------------------------------------------------- 
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
      const auto omegaBGKTot = newtonian.omegaBGK(tauFlNode, fTotNode, rhoTotNode, velNode, u2, cu, ForceField(0, nodeNo), 0);

      
      
      //const auto trE_Node = newtonian.trE();
      //------------------------------------------------------------------------------------- 
      
      const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tauFlNode, cu, uF, cF);
      
      const auto deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tauFlNode, cu, u2, Q(0, nodeNo));
      
      fTot.set(0, nodeNo) = fTotNode + deltaOmegaF + deltaOmegaQ0 + omegaBGKTot;
      
      
      const auto tauPhaseField = 1.0;
	
      for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	const auto gNode = g(fieldNo, nodeNo);
	const lbBase_t rhoDNode = calcRho<LT>(gNode);
	//                              Save density for printing
        //------------------------------------------------------------------------------------- Save density for printing
	rhoD(fieldNo, nodeNo) = rhoDNode;
      }


      
	
      lbBase_t Rtmp = 0.0; //2*(H*LT::c2*(rhoTot(0, nodeNo) - rhoD(1, nodeNo))*rhoTot(0, nodeNo) - rhoD(0, nodeNo))*rhoRel(0, nodeNo);//funker 
      
      
      
      lbBase_t Htemp = 0;

      
      
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {	        		
	  
	    
	//                              Fix Diff 0 conc in phase field 0   
        //------------------------------------------------------------------------------------- Fix Diff 0 conc in phase field 1  

	
	if(i>= 10000 && fieldNo==0 && (rhoRel(0, nodeNo)> (1.0-1e-4)
				      || ( rhoRel(0, nodeNo) > 0 && (kappa(0, nodeNo)*kappa(0, nodeNo))> 2.4 ) )){

	  lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rhoTot(0, nodeNo) - rhoD(1, nodeNo));

	  //lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rho(0, nodeNo) - rhoD(1, nodeNo))/(1-H*LT::c2*rhoTot(0, nodeNo));//Crashes!
	  
	  //Rtmp = 2*(rhoD0Fix - rhoD(0, nodeNo))*rhoRel(0, nodeNo);//funker

	  //Rtmp = 2*(H*LT::c2*(rhoTot(0, nodeNo) - rhoD(1, nodeNo))*rhoTot(0, nodeNo) - rhoD(0, nodeNo))*rhoRel(0, nodeNo);
	  //Rtmp = 2*(H*LT::c2*(rhoTot(0, nodeNo) - rhoD(1, nodeNo))*rhoTot(0, nodeNo) - rhoD(0, nodeNo))/rhoRel(0, nodeNo);
	  
	  //Rtmp = 2*(H*LT::c2*(rhoTot(0, nodeNo) - rhoD(1, nodeNo))*rhoTot(0, nodeNo) - rhoD(0, nodeNo))*rhoRel(0, nodeNo)*rhoRel(0, nodeNo);

	  
	  //Rtmp = 2*(rhoD0Fix - rhoD(0, nodeNo))*rhoRel(0, nodeNo);

	  /*
	  Rfield(0, nodeNo) = -Rtmp;
	  rhoD(0, nodeNo) += -0.5*Rfield(0, nodeNo);
	  phiD(0, nodeNo) += -0.5*Rfield(0, nodeNo)/rhoTot(0, nodeNo);
	  rho(0, nodeNo) += 0.5*Rfield(0, nodeNo);
	  phi(0, nodeNo) += 0.5*Rfield(0, nodeNo)/rhoTot(0, nodeNo);
	  */
	  
	}
	
	
	
	/*
	if(i>= 5000 && fieldNo==0){
	  Htemp = H;
	  Rtmp = kinConst*(H*(1-phiD(1, nodeNo)-phi(1, nodeNo))*phi(1, nodeNo)*(kappa(0, nodeNo)*sigma[1]*phi(0, nodeNo) + LT::c2*rhoTot(0, nodeNo)) 
	  - phiD(2, nodeNo)*phi(0, nodeNo));

	  Rfield(0, nodeNo) = -Rtmp;
	  //rhoD(0, nodeNo) += -0.5*Rfield(0, nodeNo);
	  //phiD(0, nodeNo) += -0.5*Rfield(0, nodeNo)/rhoTot(0, nodeNo);
	  rho(0, nodeNo) += 0.5*Rfield(0, nodeNo);
	  phi(0, nodeNo) += 0.5*Rfield(0, nodeNo)/rhoTot(0, nodeNo);
	  rhoD(2, nodeNo) += -0.5*Rfield(0, nodeNo);
	  phiD(2, nodeNo) += -0.5*Rfield(0, nodeNo)/rhoTot(0, nodeNo);
	  rho(1, nodeNo) += -0.5*Rfield(0, nodeNo);
	  phi(1, nodeNo) += -0.5*Rfield(0, nodeNo)/rhoTot(0, nodeNo);
	}
	*/

	
	const auto deltaOmegaR1   = calcDeltaOmegaR<LT>(tauPhaseField, cu, Rfield(fieldNo, nodeNo));
	
	const auto fNode = f(fieldNo, nodeNo);
	const auto rhoNode = rho(fieldNo, nodeNo);
	
	const auto feqNode = calcfeq<LT>(rhoNode, u2, cu);
	const auto omegaBGK = calcOmegaBGK_TEST<LT>(fNode, feqNode, tauPhaseField);        
	
	LbField<LT> deltaOmegaFDiff(1,1);
	deltaOmegaFDiff.set(0 ,0) = calcDeltaOmegaFDiff<LT>(tauPhaseField, phi(fieldNo, nodeNo)/*rhoRel(fieldNo, nodeNo)*/, cu, uF, cF);  // LBcollision
	
	//                                   Recoloring step
        //------------------------------------------------------------------------------------- Recoloring step
	int field_k_ind = (fieldNo*(fieldNo-1))/2;
	deltaOmegaRC.set(0, fieldNo) = 0;
	
	for (int field_l = 0; field_l < fieldNo; ++field_l) {
	  const int ind_F = field_k_ind + field_l;
	  const int ind_sigmaBeta = fieldNo*nFluidFields + field_l;
	  const lbBase_t IFT_threshold = 1.0;//std::min(1e6*rhoRel(field_l, nodeNo)*rhoRel(fieldNo, nodeNo),1.0);
	  //const lbBase_t IFT_threshold = (FNorm(ind_F, nodeNo)>1e-3);
	  
	  const auto cn = LT::cDotAll(F(ind_F, nodeNo)/(FNorm(ind_F, nodeNo)+(FNorm(ind_F, nodeNo)<lbBaseEps)));

	  
	  deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, 2*sigma[ind_sigmaBeta]/**rhoTot(0, nodeNo)*/*IFT_threshold, FNorm(ind_F, nodeNo), cn);
	  
	  //const auto delta_Dirac_x_2= 3*kx2_test*beta[sigmaBeta_ind]*LT::c2Inv*rhoRel(fieldNo, nodeNo)*rhoRel(fieldNo, nodeNo)*rhoRel(field_l, nodeNo)*rhoRel(field_l, nodeNo);
	  //const auto delta_Dirac_x_2= 2*2*kx2_test*beta[sigmaBeta_ind]*LT::c2Inv*rhoRel(fieldNo, nodeNo)*rhoRel(field_l, nodeNo);
	  
	  //deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, 2*sigma[sigmaBeta_ind]*rhoTot(0, nodeNo)/**IFT_threshold*/, delta_Dirac_x_2, cn);
	  
	  //deltaOmegaRC.set(0, fieldNo) += rhoNode*beta[ind_sigmaBeta]*rhoRel(field_l, nodeNo)*cn*cNormInv;
	  deltaOmegaRC.set(0, fieldNo) += rhoNode*beta[ind_sigmaBeta]*rhoRel(field_l, nodeNo)*cn;
	  //deltaOmegaRC.set(0, fieldNo) += rhoRelNode(0,fieldNo)*beta[ind_sigmaBeta]*rhoRel(field_l, nodeNo)*cn;
	}
	
	for (int field_l = fieldNo + 1; field_l < nFluidFields; ++field_l) {
	  const int field_k_ind = (field_l*(field_l-1))/2;
	  const int ind_F =  field_k_ind + fieldNo;
	  const int ind_sigmaBeta = fieldNo*nFluidFields + field_l;
	  
	  const auto cn = LT::cDotAll(F(ind_F, nodeNo)/(FNorm(ind_F, nodeNo)+(FNorm(ind_F, nodeNo)<lbBaseEps)));
	  
	  
	  
	  //deltaOmegaRC.set(0, fieldNo) -= rhoNode*beta[ind_sigmaBeta]*rhoRel(field_l, nodeNo)*cn*cNormInv;
	  deltaOmegaRC.set(0, fieldNo) -= rhoNode*beta[ind_sigmaBeta]*rhoRel(field_l, nodeNo)*cn;  
	  //deltaOmegaRC.set(0, fieldNo) -= rhoRelNode(0,fieldNo)*beta[ind_sigmaBeta]*rhoRel(field_l, nodeNo)*cn;  
	}
	    
	//if(fluidIntIndNode>1)
	//  deltaOmegaST.set(0 ,0) = deltaOmegaST(0,0)*1.0e-2;
	
		
	deltaOmegaRC.set(0, fieldNo) *= wAll;
	      
	/*
	//                                     TIME DERIVATIVE
        //------------------------------------------------------------------------------------- TIME DERIVATIVE
	const auto J_RC = 1.0*LT::qSumC(deltaOmegaRC(0, fieldNo)); 
	const auto J_RC_deriv = J_RC - J_RC_prev(fieldNo, nodeNo);
	
	J_RC_prev.set(fieldNo, nodeNo)= J_RC;
	
	deltaOmegaFDiff.set(0 ,0) += calcDeltaOmegaFDiff<LT>(1.0, 1.0, cu, 0.0, LT::cDotAll(J_RC_deriv));  // LBcollision
	*/
	
	//deltaOmegaRC.set(0, fieldNo) *= wAll*rhoRel(fieldNo, nodeNo);
	//deltaOmegaRC.set(0, fieldNo) *= rhoNode*feqTotRel0Node;
	
	//fTot.set(0, nodeNo) += deltaOmegaST(0, 0);
	
	/*
	  const auto TCap = tauFlNode*LT::qSumCCLowTri(deltaOmegaST(0, 0));
	  VectorField<LT> Tgradphi(1, 1);
	  Tgradphi.set(0,0)= LT::contractionLowTriVec(TCap, grad<LT>(phi, fieldNo, nodeNo, grid));
	  deltaOmegaFDiff.set(0 ,0) += calcDeltaOmegaFDiff<LT>(1, 1.0, cu, 0.0, LT::cDotAll(Tgradphi(0, 0)));  // LBcollision
	*/
	
	//                                 Calculate lb field
        //------------------------------------------------------------------------------------- Calculate lb field
	f.set(fieldNo, nodeNo) = fNode + omegaBGK + deltaOmegaRC(0, fieldNo) + deltaOmegaFDiff(0 ,0) + deltaOmegaR1;
	// Collision and propagation
	//fTmp.propagateTo(fieldNo, nodeNo, f(fieldNo, nodeNo), grid);
	
	//fTot.set(0, nodeNo) +=  deltaOmegaST(0, 0);
      }
      // Collision and propagation
      
	
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {       
	//                               Collision and propagation
	//------------------------------------------------------------------------------------- Collision and propagation
	fTmp.propagateTo(fieldNo, nodeNo, f(fieldNo, nodeNo) + deltaOmegaST(0, 0)*phi(fieldNo, nodeNo)*tauFlNode/tauPhaseField, grid);
      }
      
      
      //                               Collision and propagation
      //------------------------------------------------------------------------------------- Collision and propagation
      //fTot.set(0, nodeNo) += deltaOmegaF + omegaBGKTot;
      fTotTmp.propagateTo(0, nodeNo, fTot(0, nodeNo) + deltaOmegaST(0, 0), grid);
      
      
      
      //=====================================================================================
      //
      //                                  DIFFUSION
      //
      //=====================================================================================
      /*
      if(i>= 5000 && rhoRel(0, nodeNo)< 0.98){
	Rtmp = 2*(FixedConcFactor*H*LT::c2*rhoTot(0, nodeNo) - rhoD(0, nodeNo));
	Rfield(0, nodeNo) = -Rtmp;
	rhoD(0, nodeNo) += -0.5*Rfield(0, nodeNo);
	phiD(0, nodeNo) += -0.5*Rfield(0, nodeNo)/rhoTot(0, nodeNo);
	
      }	
      */
      
      
      
      ScalarField tauDiff_aveNode(nDiffFields,1);
      
      //lbBase_t tauDiff_aveNode = LT::c2Inv/diff_ave_inv + 0.5;
      for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	
	lbBase_t diff_ave_inv = 0.0;
	for (int fluidNo=0; fluidNo<nFluidFields; ++fluidNo){
	  diff_ave_inv += rhoRel(fluidNo, nodeNo)*diff_coef_inv[fieldNo*nFluidFields+fluidNo];
	  //diff_ave_inv += rhoRel(fluidNo, nodeNo)*diff_coef[fieldNo*nFluidFields+fluidNo];
	}
	//diff_ave_inv=1/diff_ave_inv;

	//diff_ave_inv *=(1-0.5*FNorm(0, nodeNo));
	//diff_ave_inv += 0.5*FNorm(0, nodeNo)*1/0.16667;

	
	//if(FNorm(0, nodeNo)>1e-3)
	//  diff_ave_inv = 1/0.16667;

	//for(int cnt=0; cnt<FNorm.num_fields(); ++cnt){
	//  if(FNorm(cnt, nodeNo)>1e-3){
	//    diff_ave_inv = 1/0.16667;
	//    break;
	//  }
	//}

	
	if(fluidIntIndNode!=0)
	  diff_ave_inv = 1/0.16667;
	
	
	
	tauDiff_aveNode(fieldNo, 0) = LT::c2Inv/diff_ave_inv + 0.5;
	
	const auto gNode = g(fieldNo, nodeNo);
	
	// Save density and velocity for printing
        //------------------------------------------------------------------------------------- 
	//rhoD(fieldNo, nodeNo) = rhoDNode;

	lbBase_t Rtmp2 = 0.0;
	/*
	if(i>2000 && i<4000 && fieldNo==1 && rhoRel(0, nodeNo)> 0.99){
	  Rtmp2 = 2*(InitialSaltConc*rhoTot(0, nodeNo) - rhoD(1, nodeNo));
	  
	  rhoD(1, nodeNo) += 0.5*Rtmp2;
	  phiD(1, nodeNo) += 0.5*Rtmp2/rhoTot(0, nodeNo);
	  
	}
	*/
	
	const lbBase_t rhoDNode = rhoD(fieldNo, nodeNo); //calcRho<LT>(gNode);
	
	const auto geqNode = calcfeq<LT>(rhoDNode, u2, cu);
	const auto omegaBGK_Diff = calcOmegaBGK_TEST<LT>(gNode, geqNode, tauDiff_aveNode(fieldNo, 0));
	lbBase_t RNode = 0.0;
	if(fieldNo==0)
	  RNode = -Rfield(0, nodeNo);
	if(fieldNo==1)
	  RNode = Rtmp2;
	const auto deltaOmegaRDiff   = calcDeltaOmegaR<LT>(tauDiff_aveNode(fieldNo, 0), cu, RNode);
	
	
	
	g.set(fieldNo, nodeNo) = gNode + omegaBGK_Diff + deltaOmegaRDiff;
      }
	
      //------------------------------------------------------
      //Interaction between diffusive fields and fluid fields
      //------------------------------------------------------
      
      LbField<LT> deltaOmegaDI(1, nDiffFields);
      for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	deltaOmegaDI.set(0, fieldNo) = 0;
      }
	
      //const CGAttributes<LT> cgat(nFluidFields, nodeNo, cNormInv, Gamma0, rhoRelNode, rhoRel, grid);
      
      lbBase_t W;//, W_1, W_2, W1, W2;
      int diffPhaseInd;
      int solventPhaseInd;
      LbField<LT> cosPhiTmp(1, 1);
      LbField<LT> cDotFRCTmp(1, 1);
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
	*/
      /*
	diffPhaseInd = 0;
	solventPhaseInd = 0;
	//soluble in phase 0
	W = rhoRel(0, nodeNo)-1;   
	cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1)+((cgat.FNorm_(0,0)+cgat.FNorm_(0,1))<lbBaseEps));
	
	deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0); 
	//-----------------------------------------------------------------------------------------------------------------------------------
	*/
      
      diffPhaseInd = 1;
	      
      //soluble in phase 0
      solventPhaseInd = 0;
      W = rhoRel(solventPhaseInd, nodeNo) - 1;
      //W = rho(0,nodeNo)/(rho(0,nodeNo)+rho(1,nodeNo)) - 1;
      //cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1));
      cosPhiTmp.set(0, 0) = cgat.cosPhi_(0, 0);
      cDotFRCTmp.set(0, 0) = cgat.cDotFRC_(0, 0);
      //cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1)+((cgat.FNorm_(0,0)+cgat.FNorm_(0,1))<lbBaseEps));
      deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0)*cNormTmp;
      //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 1]*W*cDotFRCTmp(0, 0);

      //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 1);
      
      
      /*
	diffPhaseInd = 0;
	//not soluble in phase 0
	solventPhaseInd = 2;
	//W = rhoRel(solventPhaseInd, nodeNo);
	W = - rhoRel(2, nodeNo);
	//W = rho(0,nodeNo)/(rho(0,nodeNo)+rho(1,nodeNo)) - 1;
	//cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1));
	cosPhiTmp.set(0, 0) = - cgat.cosPhi_(0, 0);
	deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 0]*W*cosPhiTmp(0, 0);   
      */

      diffPhaseInd = 2;

      //soluble in phase 1
      solventPhaseInd = 1;
      //W = rhoRel(solventPhaseInd, nodeNo) - 1;
      //Htemp=0;
      W = (4*kinConst*Htemp*(1-phiD(1, nodeNo)-phi(1, nodeNo))*phi(1, nodeNo)*(kappa(0, nodeNo)*sigma[1]*phi(0, nodeNo) + LT::c2*rhoTot(0, nodeNo))/(phiD(2, nodeNo)+(phiD(2, nodeNo)<lbBaseEps)) - phi(0, nodeNo));
      //W = rho(0,nodeNo)/(rho(0,nodeNo)+rho(1,nodeNo)) - 1;
      //cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1));
      cosPhiTmp.set(0, 0) = -cgat.cosPhi_(0, 0);
      //cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1)+((cgat.FNorm_(0,0)+cgat.FNorm_(0,1))<lbBaseEps));
      deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 0]*W*cosPhiTmp(0, 0);
      
      //deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solventPhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 1);
      
      /*
	  
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
	//diffPhaseInd = 5;
	*/
      
      /*
	diffPhaseInd = 0;
	solventPhaseInd = 1;
	//soluble in phase 0 & 1 (or, i.e., not soluble phase 2) 
	//W = rhoRel(0, nodeNo) - 1;
	
	W_1 = -rhoRel(2, nodeNo); // At interface 0-2, not soluble in phase 2
	W_2 = -rhoRel(2, nodeNo); // At interface 1-2, not soluble in phase 2
	//deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 1]*W_1*cgat.cosPhi_(0, 0);
	cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,2)*cgat.cosPhi_(0, 2) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,2)+cgat.FNorm_(0,1)+((cgat.FNorm_(0,2)+cgat.FNorm_(0,1))<lbBaseEps));
	deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 2]*W_1*cosPhiTmp(0, 0);
	//deltaOmegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 1*nFluidFields + 2]*W_2*cgat.cosPhi_(0, 2);
	*/
      //-----------------------------------------------------------------------------------------------------------------------------------
      /*
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
	  
	  
	LbField<LT> deltaOmegaFDiff(1,1);
	deltaOmegaFDiff.set(0 ,0) = calcDeltaOmegaFDiff<LT>(tauDiff_aveNode(fieldNo, 0), rhoD(fieldNo, nodeNo)/rhoTot(0, nodeNo), cu, uF, cF);
	
	/*
	//TIME DERIVATIVE
	const auto J_DI = tauDiff_aveNode(fieldNo, 0)*LT::qSumC(deltaOmegaDI(0, fieldNo));
	
	const auto J_DI_deriv = J_DI - J_DI_prev(fieldNo, nodeNo);
	
	J_DI_prev.set(fieldNo, nodeNo)= J_DI;
	
	deltaOmegaFDiff.set(0 ,0) += calcDeltaOmegaFDiff<LT>(tauDiff_aveNode(fieldNo, 0), 1.0, cu, 0.0, LT::cDotAll(J_DI_deriv));  // LBcollision
	*/
	
	//if (fieldNo==0){
	/*
	  const auto TCap = tauFlNode*LT::qSumCCLowTri(deltaOmegaST(0, 0));
	  VectorField<LT> Tgradphi(1, 1);
	  Tgradphi.set(0,0)= LT::contractionLowTriVec(TCap, grad<LT>(phiD, fieldNo, nodeNo, grid));
	  deltaOmegaFDiff.set(0 ,0) += calcDeltaOmegaFDiff<LT>(tauDiff_aveNode(fieldNo, 0), 1.0, cu, 0.0, LT::cDotAll(Tgradphi(0, 0)));  // LBcollision
	*/	
	//}
	
	deltaOmegaDI.set(0, fieldNo) *= wAll*rhoD(fieldNo, nodeNo);
	  
	gTmp.propagateTo(fieldNo, nodeNo, g(fieldNo, nodeNo) + deltaOmegaDI(0, fieldNo) + deltaOmegaFDiff.set(0 ,0)
			 + deltaOmegaST(0, 0)*phiD(fieldNo, nodeNo)*tauFlNode/tauDiff_aveNode(fieldNo, 0), grid);
      }
      //------------------------------------------------------------------------------------- END DIFFUSION

	
	
    }//------------------------------------------------------------------------------------- End nodes
    // Swap data_ from fTmp to f;
    //------------------------------------------------------------------------------------- 
    fTot.swapData(fTotTmp);  // LBfield
    f.swapData(fTmp);  // LBfield
    g.swapData(gTmp);  // LBfield
    
   
    //=====================================================================================
    //
    //                               BOUNDARY CONDITIONS
    //
    //=====================================================================================
    // Mpi
    //------------------------------------------------------------------------------------- 
    mpiBoundary.communicateLbField(fTot, grid);

    // Half way bounce back
    //------------------------------------------------------------------------------------- 
    bounceBackBnd.apply(fTot, grid);
    
    
    mpiBoundary.communicateLbField(f, grid);
    // Half way bounce back
    bounceBackBnd.apply(f, grid);
    
    mpiBoundary.communicateLbField(g, grid);
    // Half way bounce back
    bounceBackBnd.apply(g, grid);
    
   
    //=====================================================================================
    //
    //                                 WRITE TO FILE
    //
    //=====================================================================================
    if ( ((i % nItrWrite) == 0)  ) {
      
      output.write(i);
      if (myRank==0) {
	std::cout << "PLOT AT ITERATION : " << i << std::endl;
      }
    }
    
  } //-------------------------------------------------------------------------------------  End iterations
  
  MPI_Finalize();
  
  return 0;
}
