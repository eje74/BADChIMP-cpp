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

#include <omp.h>

//                                SET THE LATTICE TYPE
//------------------------------------------------------------------------------------- SET THE LATTICE TYPE
#define LT D2Q9

#define VTK_CELL VTK::pixel
//#define LT D3Q19
//#define VTK_CELL VTK::voxel

int main(int argc, char *argv[])
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
  std::string outputDir = chimpDir + "output";
  
  std::string exePath = argv[0];
  
  
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
  std::valarray<lbBase_t> Gamma0 = input["fluid"]["alpha"];
  std::valarray<lbBase_t> GammaNonZero = (1-LT::w0*Gamma0)/(1-LT::w0);
  
  //                                 Phase separation: fluids
  //------------------------------------------------------------------------------------- Phase separation: fluids
  //std::valarray<lbBase_t> beta = input["fluid"]["phaseseparation"]["beta"];

  std::valarray<lbBase_t> WInterface = input["fluid"]["phaseseparation"]["W"];

  std::valarray<lbBase_t> WInv(nFluidFields*nFluidFields);
  for (int i=0; i < nFluidFields; ++i) {
    for (int j=0; j < nFluidFields; ++j) {	
      WInv[i*nFluidFields + j] = 1/WInterface[i*nFluidFields + j];
    }
  }
  
  //                                 Phase separation: Diffusive mixing fluids
  //------------------------------------------------------------------------------------- Phase separation: Diffusive mixing fluids

  std::valarray<lbBase_t> diffCoefFluids = input["fluid"]["phaseseparation"]["diffCoef"];
  
  //                                        Surface tension
  //------------------------------------------------------------------------------------- Surface tension
  std::valarray<lbBase_t> sigma = input["fluid"]["phaseseparation"]["sigma"];

  std::valarray<lbBase_t> sigmaSigmaTmp(nFluidFields*nFluidFields);
  for (int i=0; i < nFluidFields; ++i) {
    for (int j=0; j < nFluidFields; ++j) {
      sigmaSigmaTmp[i*nFluidFields + j] = 0;
      for (int k=0; k < nFluidFields; ++k) {
	sigmaSigmaTmp[i*nFluidFields + j] += sigma[i*nFluidFields + k]* sigma[k*nFluidFields + j];
      }
    }
  }


  std::valarray<lbBase_t> Xkl(nFluidFields*nFluidFields);
  for (int i=0; i < nFluidFields; ++i) {
    for (int j=0; j < nFluidFields; ++j) {	
      Xkl[i*nFluidFields + j] = 0.0;
    }
  }
  
  if (nFluidFields==3){
    //std::cout<<"test"<<std::endl;
    lbBase_t denom;
    denom = 2*sigma[2*nFluidFields + 0]*sigma[2*nFluidFields + 1];
    Xkl[0*nFluidFields + 1] = (sigma[2*nFluidFields + 0]*sigma[2*nFluidFields + 0] + sigma[2*nFluidFields + 1]*sigma[2*nFluidFields + 1] - sigma[0*nFluidFields + 1]*sigma[0*nFluidFields + 1])
      /(denom + (denom<lbBaseEps));
    Xkl[1*nFluidFields + 0] = Xkl[0*nFluidFields + 1];

    denom = 2*sigma[1*nFluidFields + 0]*sigma[1*nFluidFields + 2];
    Xkl[0*nFluidFields + 2] = (sigma[1*nFluidFields + 0]*sigma[1*nFluidFields + 0] + sigma[1*nFluidFields + 2]*sigma[1*nFluidFields + 2] - sigma[0*nFluidFields + 2]*sigma[0*nFluidFields + 2])
      /(denom + (denom<lbBaseEps));
    Xkl[2*nFluidFields + 0] = Xkl[0*nFluidFields + 2];

    denom = 2*sigma[0*nFluidFields + 1]*sigma[0*nFluidFields + 2];
    Xkl[1*nFluidFields + 2] = (sigma[0*nFluidFields + 1]*sigma[0*nFluidFields + 1] + sigma[0*nFluidFields + 2]*sigma[0*nFluidFields + 2] - sigma[1*nFluidFields + 2]*sigma[1*nFluidFields + 2])
      /(denom + (denom<lbBaseEps));
    Xkl[2*nFluidFields + 1] = Xkl[1*nFluidFields + 2];
  }
    
  /*
  //                                           Diffusion
  //===================================================================================== Diffusion

  // Number of diffusive fields
  //------------------------------------------------------------------------------------- Number of diffusive fields
  const int nDiffFields = input["diffusion"]["numfields"];
  
  // Diffusion coefficients
  //------------------------------------------------------------------------------------- Diffusion coefficients
  std::valarray<lbBase_t> diff_coef = input["diffusion"]["fluidinteraction"]["diffcoef"];
  std::valarray<lbBase_t> diff_coef_inv(nDiffFields*nFluidFields);
  for (int i=0; i < nDiffFields; ++i) {
    for (int j=0; j < nFluidFields; ++j) {
      diff_coef_inv[i*nFluidFields + j] = 1./diff_coef[i*nFluidFields + j];
    }
  }

  //                                 Phase separation: Diffusive fields
  //------------------------------------------------------------------------------------- Phase separation: Diffusive fields
  std::valarray<lbBase_t> betaDiff = input["diffusion"]["fluidinteraction"]["beta"];
  */
  
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

  //                                   Start time for partial mixing 
  //------------------------------------------------------------------------------------- Kinetic Reaction Constant
  const lbBase_t mixStartTime = input["Partial-Misc"]["mixStartTime"];
  
  //                                    Output directory number
  //------------------------------------------------------------------------------------- Output directory number
  std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));  
  std::string outputDir2 = outputDir + "/out" + dirNum;
  std::string progMainFilePath = input["out"]["progMainFile"];


    
  //                                Recoloration Constant 2k
  //------------------------------------------------------------------------------------- Recoloration Constant 2k
  lbBase_t kx2 = 0;
  lbBase_t cTest = 0;
  lbBase_t cTest2 = 0;
  lbBase_t cTest3 = 0;
  for (int q = 0; q < LT::nQNonZero_; ++q) {
    kx2 += LT::w[q] * LT::c(q,0)*LT::c(q,0) /  LT::cNorm[q];
    cTest += LT::w[q] * LT::c(q,0)*LT::c(q,0);
    cTest2 += LT::w[q] * LT::c(q,0)*LT::c(q,0)/  (LT::cNorm[q]*LT::cNorm[q]);
    cTest3 += LT::w[q] * LT::c(q,0)*LT::c(q,0)*LT::c(q,0)*LT::c(q,0)/  (LT::cNorm[q]*LT::cNorm[q]);
  }
  //------------------------------------------------------------------------------------- END READ FROM INPUT



  
  //=====================================================================================
  //
  //                               WRITE PARAMETERS TO SCREEN
  //
  //=====================================================================================

  if (myRank==0) {
    //                              Path to executable
    //------------------------------------------------------------------------------------- Path to executable
    std::cout<< "Path to executable: "<<exePath<<std::endl;
    

    std::cout<< "Path to main.cpp: "<<progMainFilePath<<std::endl;
    std::cout << std::endl;
    
    
    //                              Number of fluid fields
    //------------------------------------------------------------------------------------- Number of fluid fields
    std::cout << "nFluidFields = " << nFluidFields << std::endl;
    std::cout << std::endl;

    // kinematic viscosity
    //------------------------------------------------------------------------------------- kinematic viscosity
    std::cout << "kin. visc. = " << std::endl;
    for (int n = 0; n < nFluidFields; ++n){
      std::cout << " " <<  kin_visc[n];
    }
    std::cout << std::endl;

    std::cout << "kin. visc. inv. = " << std::endl;
    for (int n = 0; n < nFluidFields; ++n){
      std::cout << " " <<  kin_visc_inv[n];
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
    /*
    std::cout << "beta = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
      for (int j=0; j < nFluidFields; ++j) {
	std::cout << " " << beta[i*nFluidFields + j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    */

    std::cout << "WInterface = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
      for (int j=0; j < nFluidFields; ++j) {
	std::cout << " " << WInterface[i*nFluidFields + j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "WInv = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
      for (int j=0; j < nFluidFields; ++j) {
	std::cout << " " << WInv[i*nFluidFields + j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "diffCoefFluids = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
      for (int j=0; j < nFluidFields; ++j) {
	std::cout << " " << diffCoefFluids[i*nFluidFields + j];
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

    //                                     Xkl
    //------------------------------------------------------------------------------------- Xkl
    std::cout << "X_kl = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
      for (int j=0; j < nFluidFields; ++j) {
	std::cout << " " << Xkl[i*nFluidFields + j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    /*
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
    */
    //                                Recoloration Constant 2k
    //------------------------------------------------------------------------------------- Recoloration Constant 2k
    std::cout<<"Recoloration Constant 2k = "<< kx2 <<std::endl;
    std::cout << std::endl;

    std::cout<<"lattice Constant c2 = "<< cTest <<std::endl;
    std::cout << std::endl;


    std::cout<<"Sum over lattice unit vectors squared = "<< cTest2 <<std::endl;
    std::cout << std::endl;

    std::cout<<"Sum over lattice unit vectors ^4 = "<< cTest3 <<std::endl;
    std::cout << std::endl;
    
    //                                 Partial Miscibility: Henry's Constant
    //------------------------------------------------------------------------------------- Partial Miscibility: Henry's Constant
    std::cout << "Henry's Constant H = " << H << std::endl;
    std::cout << std::endl;
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
  //ScalarField rhoRel(nFluidFields, grid.size());
  ScalarField phi(nFluidFields, grid.size());
  ScalarField rhoTot(2, grid.size());
  
  ScalarField posField(LT::nD, grid.size());
  
  //                                      Sources
  //------------------------------------------------------------------------------------- Sources
  ScalarField Rfield(nFluidFields, grid.size());
  ScalarField Qfield(1, grid.size());

  //                             Effective kinematic viscosity
  //------------------------------------------------------------------------------------- Effective kinematic viscosity
  ScalarField kinViscEffField(1, grid.size());
  
  //                              Initiate density from file
  //------------------------------------------------------------------------------------- Initiate density from file
  setScalarAttribute(rho, "init_rho_", vtklb);
  
  /*
    for (auto nodeNo: bulkNodes) {
      rho(1, nodeNo) += rho(2, nodeNo);
      rho(2, nodeNo) = 0.0;

      if(grid.pos(nodeNo, 0) < 0.5*(vtklb.getGlobaDimensions(0)-2)){
	rho(1, nodeNo) += rho(0, nodeNo);
	rho(0, nodeNo) = 0.0;
      }
    }
  */

  
  
  
    
  //                                    Wall wettability
  //------------------------------------------------------------------------------------- Wall wettability
  setScalarAttributeWall(phi, "init_rho_wall_", vtklb, nodes);
  
  //                                   Scale wettability
  //------------------------------------------------------------------------------------- scale wettability
  normelizeScalarField(phi, solidBoundaryNodes);
  
    
    
    
  //                                       Velocity
  //------------------------------------------------------------------------------------- Velocity
  VectorField<LT> vel(1, grid.size());
  //                                   Initiate velocity
  //------------------------------------------------------------------------------------- Initiate velocity
  for (auto nodeNo: bulkNodes) {
    vel.set(0, nodeNo) = 0;
    //vel(0, 0, nodeNo) = 1e-4;
  }
    

  //                                    Phase separation
  //                             Color gradients for fluid pair
  //                     Absolute values of color gradients for fluid pair
  //------------------------------------------------------------------------------------- Phase separation
  ScalarField FNorm(nFluidFields*(nFluidFields-1)/2, grid.size());
  ScalarField divF(nFluidFields*(nFluidFields-1)/2, grid.size());
  
  VectorField<LT> F(nFluidFields*(nFluidFields-1)/2, grid.size());
  VectorField<LT> ForceField(1, grid.size());
  
  //VectorField<LT> J_RC_prev(nFluidFields, grid.size());
  //VectorField<LT> J_DI_prev(nDiffFields, grid.size());
  
  VectorField<LT> unitNormal(nFluidFields*(nFluidFields-1)/2, grid.size());
  

  //VectorField<LT>TgradphiField(nFluidFields, grid.size());
  VectorField<LT>gradphiField(nFluidFields, grid.size());

  //VectorField<LT>gradphi2Field(nFluidFields, grid.size());

  //VectorField<LT>phiGradNeq(nFluidFields, grid.size());
  
  //                         Curvature of interface for fluid pair
  //                                 (lower triangular matrix)
  //------------------------------------------------------------------------------------- Curvature of interface for fluid pair
  ScalarField kappa(nFluidFields*(nFluidFields-1)/2, grid.size());
  ScalarField kappa2(nFluidFields*(nFluidFields-1)/2, grid.size());
  //VectorField<LT> gradKappa(nFluidFields*(nFluidFields-1)/2, grid.size());
  //ScalarField nGradKappa(nFluidFields*(nFluidFields-1)/2, grid.size());
  //ScalarField cosAng(nFluidFields*(nFluidFields-1)/2, grid.size());
  
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
  /*
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
  */
  
  //                                   Advection-Diffusion 
  //===================================================================================== Advection-Diffusion 
  //                                  Density & Concentration
  //------------------------------------------------------------------------------------- Density & Concentration
  //ScalarField rhoD(nDiffFields, grid.size());
  //ScalarField phiD(nDiffFields, grid.size());

  ScalarField diffCoefEffField(1, grid.size());

  // Initiate diffusive fields
  //------------------------------------------------------------------------------------- 
  
  //                                Initiate density from file
  //------------------------------------------------------------------------------------- Initiate density from file
  //setScalarAttribute(rhoD, "init_rhoD_", vtklb);
  
  
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
      const auto velNode = vel(0, nodeNo);
      const auto rhoNode = rho(fieldNo, nodeNo);
      const auto u2 = LT::dot(velNode, velNode);
      const auto cu = LT::cDotAll(velNode);
      const auto feqNode = calcfeq<LT>(rhoNode, u2, cu);
      f.set(fieldNo, nodeNo) = feqNode;
      fTmp.set(fieldNo, nodeNo) = 0.0;
      /*
      for (int q = 0; q < LT::nQ; ++q) {
	f(fieldNo, q, nodeNo) = LT::w[q]*rho(fieldNo, nodeNo);
	fTmp(fieldNo, q, nodeNo) = 0;
      }
      */
    }
  }
  
  LbField<LT> fTot(1, grid.size()); 
  LbField<LT> fTotTmp(1, grid.size());

  //              Initiate Total density, concentrations and lb distributions
  //------------------------------------------------------------------------------------- Initiate Total density, concentrations and lb distributions
  //#pragma omp parallel for num_threads(2)
  for (auto nodeNo: bulkNodes) {
    rhoTot(0, nodeNo) = 0.0;
    for (int d = 0; d < LT::nD; ++d) {
      posField(d, nodeNo) = grid.pos(nodeNo, d);
    }
    
    for (int fieldNo=0; fieldNo < rho.num_fields(); ++fieldNo) {
      rhoTot(0, nodeNo) += rho(fieldNo, nodeNo);
    }
    for (int fieldNo=0; fieldNo < rho.num_fields(); ++fieldNo) {
      phi(0, nodeNo) = rho(fieldNo, nodeNo)/rhoTot(0, nodeNo);
    }

    /*
    for (int fieldNo=0; fieldNo < rhoD.num_fields(); ++fieldNo) {
      phiD(fieldNo, nodeNo) = rhoD(fieldNo, nodeNo)/rhoTot(0, nodeNo);
    }
    */
    /*
    for (int q = 0; q < LT::nQ; ++q) {
      fTot(0, q, nodeNo) = LT::w[q]*rhoTot(0, nodeNo);
      fTotTmp(0, q, nodeNo) = 0;
    }
    */
  }
    

  /*  
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
  */
    
  //=====================================================================================
  //
  //                                  OUTPUT VTK
  //
  //=====================================================================================
  Output<LT> output(grid, bulkNodes, outputDir2, myRank, nProcs); 
  output.add_file("lb_run");
  /*output.add_scalar_variables({"rhoTot", "rho", "rhoD", "diffCoef",       "kinVisc",       "divF", "nGradKappa", "kappa_", "kappa2_", "R",     "Q",    "phi", "TCap", "pos"}, 
			      { rhoTot,   rho,   rhoD,   diffCoefEffField, kinViscEffField, divF,   nGradKappa,   kappa,    kappa2,    Rfield,  Qfield, phi,   TCap,   posField});
  output.add_vector_variables({"vel", "F", "unitNormal", "gradKappa", "forceField", "Tgradphi",     "gradphi",   "gradphi2",    "phiGradNeq"}, 
			      { vel,   F,   unitNormal,   gradKappa,   ForceField,   TgradphiField, gradphiField, gradphi2Field, phiGradNeq });
  */
  output.add_scalar_variables({"rhoTot", "rho",  "diffCoef",       "kinVisc",       "divF",  "kappa_", "kappa2_", "R",     "Q",    "phi", "TCap", "pos"}, 
			      { rhoTot,   rho,    diffCoefEffField, kinViscEffField, divF,    kappa,    kappa2,    Rfield,  Qfield, phi,   TCap,   posField});
  output.add_vector_variables({"vel", "F", "unitNormal", "forceField", "gradphi"}, 
			      { vel,   F,   unitNormal,   ForceField,   gradphiField});

  if (myRank==0) {
    system(("mkdir "+outputDir2+"/input").c_str());
    system(("cp ./input/input.dat "+outputDir2+"/input").c_str());
    system(("cp -a ./input/mpi/ "+outputDir2+"/input").c_str());
    system(("cp "+exePath+" "+outputDir2+"/input").c_str());
    system(("cp "+progMainFilePath+" "+outputDir2+"/input").c_str());
  }

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


    //#pragma omp parallel for num_threads(2)
    
    for (auto nodeNo: bulkNodes) {
      fTot.set(0, nodeNo) = 0.0;
      for (int fieldNo=0; fieldNo < rho.num_fields(); ++fieldNo) {
	fTot.set(0, nodeNo) += f(fieldNo, nodeNo);
      }
      
    }
    
    
    //           Macroscopic values : rho, rhoTot and rhoRel = rho/rhoTot
    //------------------------------------------------------------------------------------- Macroscopic values : rho, rhoTot and rhoRel = rho/rhoTot
    //calcDensityFields(rho, rhoRel, rhoTot, rhoD, phi, phiD, bulkNodes, f, fTot, g);

    
    calcDensityFields(rho, rhoTot, phi, bulkNodes, f);
    
    //mpiBoundary.communciateScalarField(rhoRel);
    mpiBoundary.communciateScalarField(phi);
    //mpiBoundary.communciateScalarField(phiD);

    int rampTimesteps = 5000;
    
    const lbBase_t ramp{ 0.5 * (1-std::cos(3.14159*std::min(i, rampTimesteps)/rampTimesteps)) };
    
    //#pragma omp parallel for num_threads(2)
    for (auto nodeNo: bulkNodes) {
      // Cacluate gradient
      //------------------------------------------------------------------------------------- 
      ScalarField phiNode(1, nFluidFields);
      for (int fieldNo = 0; fieldNo < nFluidFields; ++fieldNo){
	phiNode(0, fieldNo) = phi(fieldNo, nodeNo);
      }
      const CGAttributes<LT> cgAttr(i, nFluidFields, nodeNo, cNormInv, Gamma0, phiNode, phi, grid);

     
      
      //-----------------------------------------------------------------------------------------------------------------------------------
      //Stored interface normals are the lower triangular part of the interface normal matrix
      //and point from phase of lower phase index toward phase of higher phase index. e.g., 0->1, 0->2, 1->2 etc.
      //-----------------------------------------------------------------------------------------------------------------------------------
      
      for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	lbBase_t tmpVal = H*LT::c2*rhoTot(0,nodeNo);
	
	FNorm(cnt, nodeNo)= cgAttr.FNorm_(0, cnt);
	
	F.set(cnt, nodeNo)= cgAttr.F_(0, cnt);
	divF(cnt, nodeNo)= cgAttr.divF_(0, cnt);
	gradphiField.set(cnt, nodeNo)= cgAttr.gradNode_(0, cnt);
	
	
	
	/*
	if(rhoRelNode(0, 0)< tmpVal){
	  FNorm(cnt, nodeNo) = 0.0;
	  F.set(cnt, nodeNo) = 0.0;
	}
	*/
	
	unitNormal.set(cnt, nodeNo)= F(cnt, nodeNo)/(FNorm(cnt, nodeNo)+(FNorm(cnt, nodeNo)<lbBaseEps));
      }

      
      
      /*
      cosAng(0, nodeNo) = -1;
      cosAng(1, nodeNo) = 0;
      cosAng(2, nodeNo) = 1;
      */
      
    }

    mpiBoundary.communciateVectorField_TEST(unitNormal);

    
    mpiBoundary.communciateVectorField_TEST(F);
    mpiBoundary.communciateScalarField(FNorm);

    
    
    for (auto nodeNo: bulkNodes) {
      //                      Lower triangular traversing of kappa 
      //------------------------------------------------------------------------------------- Lower triangular traversing of kappa and IFT
      int cnt = 0;
      
      
      for (int fieldNo_k=0; fieldNo_k<nFluidFields; ++fieldNo_k) {
	for (int fieldNo_l = 0; fieldNo_l < fieldNo_k; ++fieldNo_l) {
	  kappa(cnt, nodeNo) = - div_test2<LT>(unitNormal, cnt, nodeNo, grid);

	 

	  cnt++;
	}
      }
    }
    
    
    
    

    

    
    
    //Main calculation loop
    //-------------------------------------------------------------------------------------
    //#pragma omp parallel for num_threads(4)
    for (auto nodeNo: bulkNodes) {
      
      Qfield(0, nodeNo)=0.0;


      
      
	
      
      
      // Cacluate gradient
      //------------------------------------------------------------------------------------- 
      ScalarField rhoRelNode(1, nFluidFields);
      for (int fieldNo = 0; fieldNo < nFluidFields; ++fieldNo){
	Rfield(fieldNo, nodeNo)=0.0;
	//rhoRelNode(0, fieldNo) = rhoRel(fieldNo, nodeNo);
      }

      
      
      
    
      
      
      //const CGAttributes<LT> cgAttr(nFluidFields, nodeNo, cNormInv, Gamma0, rhoRelNode, rhoRel, grid);

      /*
      //                      Fixed Pressure Source & Fixed Phase Field 1 Source 
      //------------------------------------------------------------------------------------- Fixed Pressure Source & Fixed Phase Field 1 Source 
      if(rhoRel(1, nodeNo)> 0.99999 && rhoRel(1, nodeNo)< 0.99999999){
	Qfield(0, nodeNo) = 2*(1.0 - rhoTot(0, nodeNo));
	Rfield(1, nodeNo) = 2*(1.0*rhoRel(1, nodeNo) - rho(1, nodeNo));
	
	rhoTot(0, nodeNo) += 0.5*Qfield(0, nodeNo);
	rho(1, nodeNo) += 0.5*Rfield(1, nodeNo);
	phi(1, nodeNo) = rho(1, nodeNo)/rhoTot(0, nodeNo);
      }
      */

      
      
      
      
      VectorField<LT> IFTforceNode(1,1);
      IFTforceNode.set(0 ,0) = 0;
      ForceField.set(0, nodeNo)=0.0;
      //ForceField(0, 0, nodeNo)=1e-3;

      //VectorField<LT> unitNormalOrder1Node(1,1);
      

      
      
      
      //                      Lower triangular traversing of kappa, IFT, and DEff
      //------------------------------------------------------------------------------------- Lower triangular traversing of kappa, IFT, and DEff
      int cnt = 0;
      diffCoefEffField(0, nodeNo) = 0.0;
      kinViscEffField(0, nodeNo) = 0.0;
      
      lbBase_t TmpDenominator = 0.0;
      for (int fieldNo_k=0; fieldNo_k<nFluidFields; ++fieldNo_k) {
	for (int fieldNo_l = 0; fieldNo_l < fieldNo_k; ++fieldNo_l) {
	  const int sigmaBeta_ind = fieldNo_k*nFluidFields + fieldNo_l;

	  //lbBase_t norm1 = vecNorm<LT>(gradphiField(fieldNo_k, nodeNo));
	  //lbBase_t norm2 = vecNorm<LT>(gradphiField(fieldNo_l, nodeNo));

	  //lbBase_t norm1 = phi(fieldNo_k, nodeNo);
	  //lbBase_t norm2 = phi(fieldNo_l, nodeNo);

	  lbBase_t norm1 = phi(fieldNo_k, nodeNo)*vecNorm<LT>(gradphiField(fieldNo_k, nodeNo));
	  lbBase_t norm2 = phi(fieldNo_l, nodeNo)*vecNorm<LT>(gradphiField(fieldNo_l, nodeNo));
	  
	  diffCoefEffField(0, nodeNo) += norm1*norm2*diffCoefFluids[sigmaBeta_ind];
	  TmpDenominator += norm1*norm2;

	  
	  
	  //kappa(cnt, nodeNo) = - div_test2<LT>(unitNormal, cnt, nodeNo, grid);

       
	  
	  //kappa(cnt, nodeNo) = - div_test<LT>(unitNormal, cnt, nodeNo, grid);

	 
	  //kappa2(cnt, nodeNo) = -(div_test2<LT>(F, cnt, nodeNo, grid) - LT::dot(grad<LT>(FNorm, cnt, nodeNo, grid),unitNormal(cnt, nodeNo)) );
	  
	  kappa2(cnt, nodeNo) = kappa(cnt, nodeNo);
	  
	  
	  if (FNorm(cnt, nodeNo) < 5e-4)
	    kappa2(cnt, nodeNo) = 0.0;
	  
	  
	  //if (FNorm(cnt, nodeNo) < 1e-4)
	  //  kappa2(cnt, nodeNo) = 0.0;
	  
	  const lbBase_t IFT_threshold = 1.0; //std::min(1e6*rhoRel(fieldNo_k, nodeNo)*rhoRel(fieldNo_l, nodeNo),1.0);
	  
	  //if (FNorm(cnt, nodeNo) > 2.5e-3 && (kappa(cnt, nodeNo)*kappa(cnt, nodeNo))< 1.8) 

	  
	  //IFTforceNode.set(0 ,0) += 0.5*rhoTot(0, nodeNo)*sigma[sigmaBeta_ind]*kappa(cnt, nodeNo)*F(cnt, nodeNo)/**IFT_threshold*/;

	  //IFTforceNode.set(0 ,0) += 0.5*sigma[sigmaBeta_ind]*kappa(cnt, nodeNo)*F(cnt, nodeNo)/**IFT_threshold*/;

	  /*
	  IFTforceNode.set(0 ,0) += sigma[sigmaBeta_ind]*kappa(cnt, nodeNo)*4*WInv[sigmaBeta_ind]*rhoRel(fieldNo_k, nodeNo)*rhoRel(fieldNo_l, nodeNo)*unitNormal(cnt, nodeNo);
	  */
	  
	  
	  
	  
	  //IFTforceNode.set(0 ,0) += 0.5*sigma[sigmaBeta_ind]*kappa(cnt, nodeNo)*F(cnt, nodeNo)/**IFT_threshold*/;

	  

	  //lbBase_t rhoInit=1.0;
	  //Rfield(fieldNo_k, nodeNo)+=-0.25*beta[sigmaBeta_ind]*rhoInit*LT::c2*tauPhaseField*kappa2(cnt, nodeNo);
	    
	  cnt++;
	}
      }
      diffCoefEffField(0, nodeNo) /= (TmpDenominator + (TmpDenominator<lbBaseEps));
      diffCoefEffField(0, nodeNo) += 0.5*LT::c2*(TmpDenominator<lbBaseEps);
      
      ForceField.set(0, nodeNo) +=  bodyForce(0, 0);
      ForceField.set(0, nodeNo) += IFTforceNode(0, 0);

      
      
      //---------------------------------------------
      const auto tauPhaseField =LT::c2Inv*diffCoefEffField(0, nodeNo) + 0.5 ; //0.7;//0.875;//1.0;//1.5;//2.0/3.0;//1.0; //10.0;
      //---------------------------------------------
      

      //                              Calculate local fluid viscosity
      //------------------------------------------------------------------------------------- Calculate local fluid viscosity
      /*
      lbBase_t visc_inv = 0.0;
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo){
	visc_inv += phi(fieldNo, nodeNo)*kin_visc_inv[fieldNo];
      }
      
      lbBase_t viscNode = 1/visc_inv;
      */
      lbBase_t viscNode = 1.0;
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo){
	viscNode *= pow(kin_visc[fieldNo], phi(fieldNo, nodeNo));
      }
      lbBase_t visc_inv = 1/viscNode;
	
      kinViscEffField(0, nodeNo) = viscNode;

      
      //Outlet velocity control
      //if ( grid.pos(nodeNo, 0)>= vtklb.getGlobaDimensions(0) -4 ){
      /*
      if ( grid.pos(nodeNo, 0) >= vtklb.getGlobaDimensions(0)-7
	   && grid.pos(nodeNo, 0) <= vtklb.getGlobaDimensions(0)-5
	   && grid.pos(nodeNo, 1)>=2 && grid.pos(nodeNo, 1)<= height(0, nodeNo)-1 ){	
	auto velNodeTmp = calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo)), ForceField(0, nodeNo));
	lbBase_t velThresh = 5e-2;
	if (sqrt( LT::dot(velNodeTmp, velNodeTmp))>velThresh)
	  ForceField.set(0, nodeNo) = 2*(velThresh - LT::qSumC(fTot(0, nodeNo)));
      }
      */
     
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo){ 
	rho(fieldNo, nodeNo) += 0.5*Rfield(fieldNo, nodeNo);
	
      }

      

      
   
      
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

      
      //if(fluidIntIndNode!=0)
      //visc_inv = 1/0.16667;
      
      const lbBase_t tauFlNode = LT::c2Inv/visc_inv + 0.5;

      
      
      
      
      
      // Calculate velocity
      //------------------------------------------------------------------------------------- 
      // Copy of local velocity distribution
      //-------------------------------------------------------------------------------------
      
      
      auto velNode = calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo)), ForceField(0, nodeNo));
      
     
	
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
      
    
      
      
      
      LbField<LT> deltaOmegaRC(1, nFluidFields);
	
      LbField<LT> deltaOmegaST(1,1);
      deltaOmegaST.set(0 ,0) = 0;

      LbField<LT> deltaOmegaIFPressurePerb(1,1);
      deltaOmegaIFPressurePerb.set(0 ,0) = 0;
      
      const auto rhoTotNode = rhoTot(0, nodeNo);

  

      
      //const auto trE_Node = newtonian.trE();
      //------------------------------------------------------------------------------------- 

      //const auto deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tauFlNode, cu, u2, Qfield(0, nodeNo));
      const std::valarray<lbBase_t> deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tauFlNode, cu, u2, Qfield(0, nodeNo));
      //const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tauFlNode, cu, uF, cF);

      
      LbField<LT> deltaOmegaF(1,1);
      deltaOmegaF.set(0,0) = 0.0;

      //deltaOmegaF.set(0,0) = calcDeltaOmegaFTRT<LT>(tauFlNode, tauPhaseField, 1.0, cu, uF, cF);  // LBcollision
      
      
      // Calculate the TCap correction
      //------------------------------------------------------------------------------------- 
      const auto uF_IFT = LT::dot(velNode, IFTforceNode(0, 0));
      const auto cF_IFT = LT::cDotAll(IFTforceNode(0, 0));

      /*
      std::valarray<lbBase_t> deltaOmegaGradTCap(LT::nQ);

      lbBase_t tau_factor = (1 - 0.5 / tauFlNode);

      for (int q = 0; q < LT::nQ; ++q){
	deltaOmegaGradTCap[q] = LT::w[q]*tau_factor * LT::c4Inv * ( cF_IFT[q] * cu[q] - LT::c2 * uF_IFT);
      }
      */
      //const auto deltaOmegaGradTCap = wAll * (1 - 0.5 / tauFlNode) * LT::c4Inv * ( cF_IFT * cu - LT::c2 * uF_IFT);
      
      //fTot.set(0, nodeNo) = fTotNode + omegaBGKTot + deltaOmegaQ0 + deltaOmegaF /*+ deltaOmegaGradTCap*/;
      
      



      /*
      for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	const auto gNode = g(fieldNo, nodeNo);
	const lbBase_t rhoDNode = calcRho<LT>(gNode);
	//                              Save density for printing
        //------------------------------------------------------------------------------------- Save density for printing
	rhoD(fieldNo, nodeNo) = rhoDNode;
      }
      */

      
	
      lbBase_t Rtmp = 0.0; //2*(H*LT::c2*(rhoTot(0, nodeNo) - rhoD(1, nodeNo))*rhoTot(0, nodeNo) - rhoD(0, nodeNo))*rhoRel(0, nodeNo);//funker 
      
      
      
      

      /*
      lbBase_t windowFunc = 1.0;

      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {
	windowFunc *= phi(fieldNo, nodeNo);
      }
      windowFunc*=100;
      windowFunc=std::min(windowFunc, 1.0);
      windowFunc=0.0;
      */
      
      //Rtmp = (H*LT::c2*rho(0, nodeNo)-rho(3, nodeNo))*4*WInv[1]*rho(0, nodeNo)*rho(1, nodeNo)/(rho(0, nodeNo)+rho(1, nodeNo));

      

      /*
      if(i>= 2000){
	lbBase_t phiA = phi(0, nodeNo) + phi(2, nodeNo) + phi(3, nodeNo);
	if((phiA > (1.0-1e-6))
	    || ( phiA > phi(3, nodeNo) && (kappa(0, nodeNo)*kappa(0, nodeNo))> 2.4 ) ){

	  lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rhoTot(0, nodeNo) - rho(2, nodeNo));
	  //lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rho(0, nodeNo)+rho(3, nodeNo));
	  Rtmp = 2*(rhoD0Fix - rho(3, nodeNo))*(phiA-phi(3, nodeNo));
	
	}
      }
      */
      if(i>= mixStartTime){
	lbBase_t phi023 = phi(0, nodeNo) + phi(2, nodeNo) + phi(3, nodeNo);
	lbBase_t phi02 = phi(0, nodeNo) + phi(2, nodeNo);
	lbBase_t phi03 = phi(0, nodeNo) + phi(3, nodeNo);
	
	lbBase_t rho0134 = rhoTot(0, nodeNo) - rho(2, nodeNo);
	lbBase_t phi0134 = rho0134/rhoTot(0, nodeNo);
	

	
	
	lbBase_t pressureDiffRemote = 0.0;


	lbBase_t deltaTmp = 0;
	
	for (int field_l = 1; field_l < nFluidFields; ++field_l) {
	  if(field_l!=3){
	    const int sigmaBeta_ind = 0*nFluidFields + field_l;
	    const int field_k_ind = (field_l*(field_l-1))/2;
	    const int ind_kappa =  field_k_ind + 0;
	    pressureDiffRemote -= sigma[sigmaBeta_ind]*kappa(ind_kappa, nodeNo)*phi(field_l, nodeNo);
	    
	    deltaTmp += 0.5*FNorm(ind_kappa, nodeNo);
	  }
	}
	lbBase_t pressureRemote = LT::c2*rhoTot(0, nodeNo) + pressureDiffRemote;
	lbBase_t rhoTotDiffRemote = LT::c2Inv*pressureDiffRemote;
	
	lbBase_t rhoTotRemote = rhoTot(0, nodeNo) + rhoTotDiffRemote;

	
	rhoTot(1, nodeNo) = rhoTotRemote;

	if(phi02<1e-3)
	  rhoTot(1, nodeNo) = rhoTot(0, nodeNo);
	
	lbBase_t srcFilter = 0.0;

	/*
	if(phiA<(1-phiB) && phiB>(1-phiA)) srcFilter = 0.5;
	else if(phiB<(1-phiA)) srcFilter = phiB;
	*/
	
	srcFilter = phi023;

	
	
	if(phi(1, nodeNo)>0.5){	
	  srcFilter = phi02;	  
	}

	/*
	lbBase_t phi023Threshold = 0.99;//(1-1e-3);
	
	if(phi023>phi023Threshold){
	  srcFilter = 1.0;
	}
	else
	  srcFilter = 0.0;
	*/
	

	lbBase_t srcFilter2 = pow(srcFilter, 2);
	
	
	

	//lbBase_t rhoD0Fix=(H*LT::c2*rhoTot(0, nodeNo)*rhoWater*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2));

	//lbBase_t rhoD0Fix=(H*rhoTot(0, nodeNo)*(LT::c2*rhoWater + phi(0, nodeNo)*pressureRemote)*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2));
	//lbBase_t rhoD0Fix=(H*rhoTot(0, nodeNo)*pressureRemote*phiC*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2));
	//lbBase_t rhoD0Fix=(H*rhoTot(0, nodeNo)*(LT::c2*rhoWater + pressureRemote)*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2));
	//lbBase_t rhoD0Fix=H*rhoTot(0, nodeNo)*(LT::c2*rhoWater + pressureDiffRemote)*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2);
	//lbBase_t rhoD0Fix=H*rhoTot(0, nodeNo)*(LT::c2*rhoWater + pressureDiffRemote);

	//lbBase_t rhoD0Fix=H*rhoTot(0, nodeNo)*(LT::c2*rhoWater + pressureDiffRemote)*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2);
	//lbBase_t rhoD0Fix=H*LT::c2*rhoTotRemote*rhoTotRemote*phiWater*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2);
	//lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*rhoTotRemote*phiWater*(srcFilter2 + /*2*4*WInv[1]**/phiB*(1-phiA)) + rho(3, nodeNo)*(1-srcFilter2);
	//lbBase_t rhoD0Fix=H*LT::c2*rhoTotRemote*rhoWater*(srcFilter2 + kinConst*2*4*WInv[1]*phiB*(1-phiA)) + rho(3, nodeNo)*(1-srcFilter2);

	//rhoD0Fix += ( 2*H*LT::c2*rhoTotRemote*rhoTotRemote*(1.0-phi(2, nodeNo)) - phi(3, nodeNo) )* 4e-2*4*WInv[1]*phiB*(1-phiA);
	
	//lbBase_t rhoD0Fix= H*pressureRemote*rhoWater*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2));
	
	//lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rhoTot(0, nodeNo) - rho(2, nodeNo))*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2);

	lbBase_t rhoTotInHenry = rhoTot(0, nodeNo);//rhoTotRemote;
	lbBase_t phiWater = phi03;
	//lbBase_t rhoD0Fix= H*LT::c2*rhoTot(0, nodeNo)*(rhoWater1)*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2);
	lbBase_t rhoD0Fix= H*LT::c2*rhoTotInHenry*rhoTotInHenry*phiWater*srcFilter2 + rho(3, nodeNo)*(1-srcFilter2);
	//lbBase_t rhoD0Fix= H*rhoTotInHenry*phiWater1*(LT::c2*rhoTotInHenry*srcFilter2 + pressureDiffRemote) + rho(3, nodeNo)*(1-srcFilter2);
	//lbBase_t rhoD0Fix= H*LT::c2*rhoTotInHenry*rhoTotInHenry*phi(0, nodeNo);
	
	
	//if(phiA > 0.5)
	//  rhoD0Fix*=1.0/phiA;
	
	//if(phi(0, nodeNo) > lbBaseEps && phi(0, nodeNo) < phi(3, nodeNo))
	//  rhoD0Fix = rho(0, nodeNo);

	
	//else
	//  rhoD0Fix*=1.0/phiA;
	
	Rtmp = 2*(rhoD0Fix - rho(3, nodeNo))/**FNorm(0, nodeNo)*FNorm(2, nodeNo)*/  /**srcFilter2*/;

	//Rtmp = 2*(H*LT::c2*rhoTotInHenry*rhoTotInHenry - rho(3, nodeNo))*phi03*kinConst;
	
	
	//if (Rtmp<0 && phi(1, nodeNo))
	
	//Rtmp = (H * pressureRemote * rhoWater - rho(3, nodeNo)) * phiB*(1-phiA)*4;
	 

	/*
	if( phiA > phiAThreshold ){
	
	  lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rhoTot(0, nodeNo) - rho(2, nodeNo))/phiA;
	  //Rtmp = 2*(rhoD0Fix - rho(3, nodeNo))*phiB;
	  //lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rho(0, nodeNo) + rho(3, nodeNo))/phiA;
	  Rtmp = 2*(rhoD0Fix - rho(3, nodeNo));
	  
	  //Rtmp = 2*(rhoD0Fix - rho(3, nodeNo))*(1-phi(1, nodeNo));
        
	}
	*/
	
	/*
	else if(phi(0, nodeNo)>lbBaseEps){
	  lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rhoTot(0, nodeNo) - rho(2, nodeNo))/phiA;
	  //lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rho(0, nodeNo) + rho(3, nodeNo));
	  Rtmp = 0.666667*(rhoD0Fix - rho(3, nodeNo))*phiB*phiB;
	}
	*/
	
	/*
	else if(//kappa2(0, nodeNo)< -0.5){  
	//&&  phi(0, nodeNo) > lbBaseEps){
	//else if(kappa2(0, nodeNo)< -0.5  &&  phi(0, nodeNo) > lbBaseEps){   
	  //lbBase_t kappa0 = - div_test2<LT>(unitNormal, 0, nodeNo, grid);
	  //if(kappa0*kappa0>2.4){
	  lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rhoTot(0, nodeNo) - rho(2, nodeNo));
	  //lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rho(0, nodeNo) + rho(3, nodeNo));
	  //Rtmp = 2*(rhoD0Fix - rho(3, nodeNo))*phiA;
	  Rtmp = 2*(rhoD0Fix - rho(3, nodeNo))*phi(0, nodeNo)*(1-phi(0, nodeNo))*2;
	    //Rtmp = 2*(rhoD0Fix - rho(3, nodeNo))*(1-phi(1, nodeNo));
	    //}
	}
      */
      }

      /*
      if(i>= 2000){
	lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rhoTot(0, nodeNo) - rho(2, nodeNo));
	Rtmp = kinConst*(rhoD0Fix - rho(3, nodeNo))*FNorm(0, nodeNo)*0.5;
      }
      */
      
      
      
      //Fluid Field Loop
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {	        		

	
	
	    
	//                              Fix Diff 0 conc in phase field 0   
        //------------------------------------------------------------------------------------- Fix Diff 0 conc in phase field 1  

	
	//if(i>= 10000 && fieldNo==0 && (rhoRel(0, nodeNo)> (1.0-1e-4)
	//			      || ( rhoRel(0, nodeNo) > 0 && (kappa(0, nodeNo)*kappa(0, nodeNo))> 2.4 ) )){

	  //lbBase_t rhoD0Fix=H*LT::c2*rhoTot(0, nodeNo)*(rhoTot(0, nodeNo) - rhoD(1, nodeNo));

	  
	  
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
	  
	//}

	if(i>= mixStartTime){
	  if(fieldNo==0)
	    Rfield(fieldNo, nodeNo) = -Rtmp;
	  if(fieldNo==3)
	    Rfield(fieldNo, nodeNo) = Rtmp;
	}

	rho(fieldNo, nodeNo) += 0.5*Rfield(fieldNo, nodeNo);
	phi(fieldNo, nodeNo) += 0.5*Rfield(fieldNo, nodeNo)/rhoTot(0, nodeNo);
	
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

	
	
	
	//const auto deltaOmegaR1  = calcDeltaOmegaR<LT>(tauPhaseField, cu, Rfield(fieldNo, nodeNo));
	//const auto deltaOmegaR1  = calcDeltaOmegaRTRT<LT>(tauFlNode, tauPhaseField, cu, Rfield(fieldNo, nodeNo));
	const auto deltaOmegaR1  = calcDeltaOmegaQTRT<LT>(tauFlNode, tauPhaseField, cu, u2, Rfield(fieldNo, nodeNo));
	
	
	auto fNode = f(fieldNo, nodeNo);
	const auto rhoNode = rho(fieldNo, nodeNo);
	
	const auto feqNode = calcfeq<LT>(rhoNode, u2, cu);

	const auto cJphi = LT::cDotAll(LT::qSumC(f(fieldNo, nodeNo)));

	
	//phiGradNeq.set(0, nodeNo) = 4*beta[1]*rhoRel(0, nodeNo)*(1-rhoRel(0, nodeNo))*unitNormal(0, nodeNo);

	//phiGradNeq.set(0, nodeNo) = 2*2*LT::c2Inv*(LT::qSumC(f(0, nodeNo)/*-feqNode*/) - phi(0, nodeNo)*LT::qSumC(fTot(fieldNo, nodeNo)))
	  /*+ LT::c2Inv*beta[1]*rho(fieldNo, nodeNo)*(1-phi(fieldNo, nodeNo))*unitNormal(0, nodeNo)*/;

	/*
	if(i>0){
	  F.set(0, nodeNo) = phiGradNeq(0, nodeNo);
	  FNorm(0, nodeNo) = LT::dot(phiGradNeq(0, nodeNo),phiGradNeq(0, nodeNo)); 
	}
	*/

	
	//LB Regularization
	const auto fNeqNode =  fNode - feqNode;
	const auto PiNeqLowTri = LT::qSumCCLowTri(fNeqNode);
	const auto M_iNeq = LT::qSumC(fNeqNode);
	const auto MNeq = LT::qSum(fNeqNode);
	fNode = calcRegDist<LT>(feqNode, MNeq, M_iNeq, PiNeqLowTri);
	//LB Regularization END
	

	//const auto omegaBGK = calcOmegaBGK_TEST<LT>(fNode, feqNode, tauPhaseField);        

	const auto omegaBGK = calcOmegaBGKTRT<LT>(fNode, tauFlNode, tauPhaseField, rhoNode, u2, cu);

	//lbBase_t rho0Tmp = 1.0;
	//const auto deltaOmegaPureAdv = calcDeltaOmegaPureAdvectionTRT<LT>(tauFlNode, tauPhaseField, phi(fieldNo, nodeNo), rhoTotNode, rho0Tmp, cNormTmp);
	

	
	LbField<LT> deltaOmegaFDiff(1,1);
	deltaOmegaFDiff.set(0 ,0) = calcDeltaOmegaFTRT<LT>(tauFlNode, tauPhaseField, phi(fieldNo, nodeNo), cu, uF, cF);  // LBcollision

	
	deltaOmegaF.set(0,0) += deltaOmegaFDiff.set(0 ,0);

	//                                   Recoloring step
        //------------------------------------------------------------------------------------- Recoloring step
	int field_k_ind = (fieldNo*(fieldNo-1))/2;
	deltaOmegaRC.set(0, fieldNo) = 0;

	
	
	for (int field_l = 0; field_l < fieldNo; ++field_l) {
	  const int ind_F = field_k_ind + field_l;
	  const int ind_sigmaBeta = fieldNo*nFluidFields + field_l; 
	  
	  const lbBase_t IFT_threshold = std::min(1e6*phi(fieldNo, nodeNo)*phi(field_l, nodeNo),1.0);
	  
	  const auto cn = LT::cDotAll(F(ind_F, nodeNo)/(FNorm(ind_F, nodeNo)+(FNorm(ind_F, nodeNo)<lbBaseEps)));
	  const auto un = LT::dot(velNode, F(ind_F, nodeNo)/(FNorm(ind_F, nodeNo)+(FNorm(ind_F, nodeNo)<lbBaseEps)));

	  lbBase_t potential = phi(field_l, nodeNo);//rho(field_l, nodeNo);//

	  /*
	  if(fieldNo==4 && field_l==3){
	    potential = phi(field_l, nodeNo)*(1-phi(field_l, nodeNo)-phi(fieldNo, nodeNo));	    
	  }
	  */
	  
	  /*
	  if(field_l==3){
	    potential = phi(field_l, nodeNo)*(1-phi(field_l, nodeNo)-phi(fieldNo, nodeNo));	    
	  }
	  */

	  deltaOmegaRC.set(0, fieldNo) += calcDeltaOmegaRC5<LT>(tauFlNode, tauPhaseField, WInv[ind_sigmaBeta],  rhoTotNode, phi(fieldNo, nodeNo), potential, cn, cu, un, kx2);

        
	  if(sigma[ind_sigmaBeta]!=0.0 && i>10){
	    lbBase_t phiTotDenomInv = 1.0;
	    
	    /*
	    if(i>10){
	      lbBase_t phiTotDenom = phi(field_l,nodeNo) + phi(fieldNo, nodeNo);
	      phiTotDenom = phiTotDenom*phiTotDenom;
	    //phiTotDenom = phiTotDenom + (phiTotDenom<lbBaseEps);
	    //phiTotDenom = 1.0;

	    //lbBase_t phiTotDenomInv = (phiTotDenom>1e-8)/(phiTotDenom*phiTotDenom);
	    
	    //if(sigma[ind_sigmaBeta]!=0.0)
	    //  phiTotDenomInv = (phiTotDenom>1e-6)/phiTotDenom;
	    //phiTotDenomInv = (phiTotDenom>1e-3)/phiTotDenom;
	      phiTotDenomInv = (phi(fieldNo, nodeNo)*phi(field_l, nodeNo)>1e-6)/phiTotDenom;
	    }
	    */
	    /*
	    if(phi(field_l,nodeNo) * phi(fieldNo, nodeNo)>5e-3){
	      lbBase_t phiTotDenom = phi(field_l,nodeNo) + phi(fieldNo, nodeNo);
	      phiTotDenomInv = 1.0/(phiTotDenom*phiTotDenom);
	    }
	    */
	    /*  
	    if(fieldNo==1 && field_l==0)
	      phiTotDenomInv =1/(1-phi(3, nodeNo));
	    */
	    deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, IFT_threshold*sigma[ind_sigmaBeta], FNorm(ind_F, nodeNo)*phiTotDenomInv, cn);

	    deltaOmegaIFPressurePerb.set(0 ,0) += calcDeltaOmegaIFPressurePerturb<LT>(tauFlNode, IFT_threshold*sigma[ind_sigmaBeta], FNorm(ind_F, nodeNo)*phiTotDenomInv);
	  }
	  
	}
	
	for (int field_l = fieldNo + 1; field_l < nFluidFields; ++field_l) {
	  const int field_k_ind = (field_l*(field_l-1))/2;
	  const int ind_F =  field_k_ind + fieldNo;
	  const int ind_sigmaBeta = fieldNo*nFluidFields + field_l;
	  
	  const lbBase_t IFT_threshold = std::min(1e6*phi(fieldNo, nodeNo)*phi(field_l, nodeNo),1.0);
	  
	  
	  const auto cn = LT::cDotAll(F(ind_F, nodeNo)/(FNorm(ind_F, nodeNo)+(FNorm(ind_F, nodeNo)<lbBaseEps)));
	  const auto un = LT::dot(velNode, F(ind_F, nodeNo)/(FNorm(ind_F, nodeNo)+(FNorm(ind_F, nodeNo)<lbBaseEps)));

	  lbBase_t potential = phi(field_l, nodeNo);// rho(field_l, nodeNo); //
	  /*
	  if(fieldNo==3 && field_l==4){
	    potential = phi(field_l, nodeNo)*(1-phi(fieldNo, nodeNo)-phi(field_l, nodeNo));
	  }
	  */
	  /*
	  if(fieldNo==3){
	    potential = phi(field_l, nodeNo)*(1-phi(fieldNo, nodeNo)-phi(field_l, nodeNo));
	  }
	  */
	  
	  deltaOmegaRC.set(0, fieldNo) -= calcDeltaOmegaRC5<LT>(tauFlNode, tauPhaseField, WInv[ind_sigmaBeta], rhoTotNode, phi(fieldNo, nodeNo), potential, cn, cu, un, kx2);
	  
	
	  if(sigma[ind_sigmaBeta]!=0.0 && i>10){
	    lbBase_t phiTotDenomInv = 1.0;
	    /*
	    if(i>10){
	      lbBase_t phiTotDenom = phi(field_l,nodeNo) + phi(fieldNo, nodeNo);
	      phiTotDenom = phiTotDenom*phiTotDenom;
	    //phiTotDenom = phiTotDenom + (phiTotDenom<lbBaseEps);
	    //phiTotDenom = 1.0;

	    
	    //if(sigma[ind_sigmaBeta]!=0.0)
	      phiTotDenomInv = (phi(fieldNo, nodeNo)*phi(field_l, nodeNo)>1e-6)/phiTotDenom;
	    }
	    */
	    /*
	    if(phi(field_l,nodeNo) * phi(fieldNo, nodeNo)>5e-3){
	      lbBase_t phiTotDenom = phi(field_l,nodeNo) + phi(fieldNo, nodeNo);
	      phiTotDenomInv = 1.0/(phiTotDenom*phiTotDenom);
	    }
	    */
	    /*
	    if(fieldNo==0 && field_l==1)
	      phiTotDenomInv =1/(1-phi(3, nodeNo));
	    */
	    deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, IFT_threshold*sigma[ind_sigmaBeta], FNorm(ind_F, nodeNo)*phiTotDenomInv, cn);
	    deltaOmegaIFPressurePerb.set(0 ,0) += calcDeltaOmegaIFPressurePerturb<LT>(tauFlNode, IFT_threshold*sigma[ind_sigmaBeta], FNorm(ind_F, nodeNo)*phiTotDenomInv);
	  }
	}
	//END Recoloring step ---------------------------------------------------------------------------------


	//                                 Calculate lb field
        //------------------------------------------------------------------------------------- Calculate lb field
	f.set(fieldNo, nodeNo) = /*fTot(0, nodeNo)*phi(fieldNo, nodeNo)*/  fNode + omegaBGK /*+ deltaOmegaPureAdv */+ deltaOmegaRC(0, fieldNo)  /*+ deltaOmegaFDiff(0 ,0)*/ + deltaOmegaR1;
	//f.set(fieldNo, nodeNo) = fTotNode*phi(fieldNo, nodeNo) + omegaBGK + deltaOmegaRC(0, fieldNo) + deltaOmegaFDiff(0 ,0) + deltaOmegaR1;
	// Collision and propagation
      }
      // Collision and propagation

	
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {       
	//                               Collision and propagation
	//------------------------------------------------------------------------------------- Collision and propagation
	//NB
	//Capillary Tensor, through deltaOmegaST, must be added here, after all fluid phases
	//have been traversed
	//NB
	//-------------------------------------------------------------------------------------
	fTmp.propagateTo(fieldNo, nodeNo, f(fieldNo, nodeNo) + phi(fieldNo, nodeNo)*deltaOmegaST(0, 0) + phi(fieldNo, nodeNo)*deltaOmegaF(0 ,0), grid);
      }
      
	
	
    }//------------------------------------------------------------------------------------- End nodes
    // Swap data_ from fTmp to f;
    //------------------------------------------------------------------------------------- 
    f.swapData(fTmp);  // LBfield
    //g.swapData(gTmp);  // LBfield
    
   
    //=====================================================================================
    //
    //                               BOUNDARY CONDITIONS
    //
    //=====================================================================================

    // Mpi
    //------------------------------------------------------------------------------------- 
    //mpiBoundary.communicateLbField(fTot, grid);
    // Half way bounce back
    //------------------------------------------------------------------------------------- 
    //bounceBackBnd.apply(fTot, grid);
    
    // Mpi
    //------------------------------------------------------------------------------------- 
    mpiBoundary.communicateLbField(f, grid);
    // Half way bounce back
    //------------------------------------------------------------------------------------- 
    bounceBackBnd.apply(f, grid);

    // Mpi
    //------------------------------------------------------------------------------------- 
    //mpiBoundary.communicateLbField(g, grid);
    // Half way bounce back
    //------------------------------------------------------------------------------------- 
    //bounceBackBnd.apply(g, grid);
    
   
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
