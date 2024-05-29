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


void zouHePressureBoundary(const std::vector<int> &bndNodes, LbField<LT> &f, const lbBase_t rhoPrime, const VectorField<LT> &force,  const Grid<LT> &grid)
{
  #pragma omp parallel for num_threads(2)
  for (const auto &nodeNo: bndNodes) {
    const std::valarray<lbBase_t> fNode = f(0, nodeNo);
    const std::valarray<lbBase_t> forceNode = force(0, nodeNo);
    const lbBase_t rho_ux = rhoPrime - 0.5*forceNode[0] - (fNode[2] + fNode[6] + fNode[8] + 2*(fNode[3] + fNode[4] + fNode[5]));
    f(0, 0, nodeNo) = fNode[4] + (2./3.)*rho_ux - (1./3.)*forceNode[0];
    f(0, 1, nodeNo) = fNode[5] + 0.5*(fNode[6] - fNode[2]) + (1./6.)*rho_ux + (5./12.)*forceNode[0] + (1./4.)*forceNode[1]; 
    f(0, 7, nodeNo) = fNode[3] + 0.5*(fNode[2] - fNode[6]) + (1./6.)*rho_ux + (5./12.)*forceNode[0] - (1./4.)*forceNode[1];
    
    // LT::qSum(f(0,nodeNo))
  }
}

void zouHePressureBoundaryRight(const std::vector<int> &bndNodes, LbField<LT> &f, const lbBase_t rhoPrime, const VectorField<LT> &force, const Grid<LT> &grid)
{
  #pragma omp parallel for num_threads(2)
  for (const auto &nodeNo : bndNodes)
  {
    const std::valarray<lbBase_t> fNode = f(0, nodeNo);
    const std::valarray<lbBase_t> forceNode = force(0, nodeNo);
    const lbBase_t rho_ux = rhoPrime - 0.5 * forceNode[0] - (fNode[2] + fNode[6] + fNode[8] + 2 * (fNode[0] + fNode[1] + fNode[7]));
    f(0, 4, nodeNo) = fNode[0] + (2. / 3.) * rho_ux - (1. / 3.) * forceNode[0];
    f(0, 3, nodeNo) = fNode[7] + 0.5 * (fNode[6] - fNode[2]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] + (1. / 4.) * forceNode[1];
    f(0, 5, nodeNo) = fNode[1] + 0.5 * (fNode[2] - fNode[6]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] - (1. / 4.) * forceNode[1];

    // LT::qSum(f(0,nodeNo))
  }
}

void zouHeConcBoundaryRight(const std::vector<int> &bndNodes, LbField<LT> &f, const ScalarField &phi,  const lbBase_t rhoPrime, const VectorField<LT> &F, const VectorField<LT> &force, const Grid<LT> &grid)
{
  #pragma omp parallel for num_threads(2)
  for (const auto &nodeNo : bndNodes)
  {
    const std::valarray<lbBase_t> forceNode = force(0, nodeNo);  
    
    std::valarray<lbBase_t> fNode = f(0, nodeNo);
    lbBase_t phi0Node = phi(0, nodeNo) + F(0, 0, nodeNo);
    if (phi0Node > 1.0) phi0Node = 1.0;
    else if (phi0Node < 0.0) phi0Node = 0.0;
    lbBase_t rho_ux = rhoPrime*phi0Node - 0.5 * forceNode[0]*phi0Node*0 - (fNode[2] + fNode[6] + fNode[8] + 2 * (fNode[0] + fNode[1] + fNode[7]));
    f(0, 4, nodeNo) = fNode[0] + (2. / 3.) * rho_ux - (1. / 3.) * forceNode[0];
    f(0, 3, nodeNo) = fNode[7] + 0.5 * (fNode[6] - fNode[2]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0]*phi0Node + (1. / 4.) * forceNode[1]*phi0Node;
    f(0, 5, nodeNo) = fNode[1] + 0.5 * (fNode[2] - fNode[6]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0]*phi0Node - (1. / 4.) * forceNode[1]*phi0Node;


    fNode = f(1, nodeNo);  
    //phi0Node = phi(1, nodeNo) - F(0, 0, nodeNo);
    phi0Node = (1-phi0Node);
    if (phi0Node > 1.0) phi0Node = 1.0;
    else if (phi0Node < 0.0) phi0Node = 0.0;
    rho_ux = rhoPrime*phi0Node - 0.5 * forceNode[0]*phi0Node*0 - (fNode[2] + fNode[6] + fNode[8] + 2 * (fNode[0] + fNode[1] + fNode[7]));
    f(1, 4, nodeNo) = fNode[0] + (2. / 3.) * rho_ux - (1. / 3.) * forceNode[0];
    f(1, 3, nodeNo) = fNode[7] + 0.5 * (fNode[6] - fNode[2]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0]*phi0Node + (1. / 4.) * forceNode[1]*phi0Node;
    f(1, 5, nodeNo) = fNode[1] + 0.5 * (fNode[2] - fNode[6]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0]*phi0Node - (1. / 4.) * forceNode[1]*phi0Node;

    
    // LT::qSum(f(0,nodeNo))
  }
}

//const lbBase_t PI = 3.14159265359;
const lbBase_t PI = 3.14159265358979323846;

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
  

  std::vector<int> OutletBoundaryNodes;
  std::vector<int> InletBoundaryNodes;

  for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {

    if (grid.pos(nodeNo, 0) == vtklb.getGlobaDimensions(0)-3)
      OutletBoundaryNodes.push_back(nodeNo);

    if (grid.pos(nodeNo, 0) == 1)
      InletBoundaryNodes.push_back(nodeNo);
    
  }

  /*
  for (auto nodeNo: OutletBoundaryNodes) {
    std::cout <<"Outlet; node no: "<<nodeNo<< std::endl;
  }
  for (auto nodeNo: InletBoundaryNodes) {
    std::cout <<"Inlet; node no: "<<nodeNo<< std::endl;
  }
  */
  
  //std::vector<int> OutletBoundaryNodes = findFluidBndNodes(nodes, const ScalarField &markerField, const lbBase_t markerVal)

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

  //                                     Fixed mean inlet flux
  //------------------------------------------------------------------------------------- Fixed inlet flux
  lbBase_t meanInletFlux = input["fluid"]["meanInletFlux"];
  lbBase_t inletXEnd = 1 + input["fluid"]["inletWidth"];
  
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
    std::cout << std::endl;

    std::cout<<"System dimensions ("<<vtklb.getGlobaDimensions(0)<<","<<vtklb.getGlobaDimensions(1)<<")"<< std::endl;
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
  ScalarField Qfield(1, grid.size());
  
  //                              Initiate density from file
  //------------------------------------------------------------------------------------- Initiate density from file
  setScalarAttribute(rho, "init_rho_", vtklb);
  
  
  
  for (auto nodeNo: bulkNodes) {
    if(/*grid.pos(nodeNo, 0) <= 2 //&& grid.pos(nodeNo, 0) >= 0.5*vtklb.getGlobaDimensions(0) //||*/ /*grid.pos(nodeNo, 0) < 0.97*vtklb.getGlobaDimensions(0)*/ /*grid.pos(nodeNo, 0) < -10*/ grid.pos(nodeNo, 0) < inletXEnd+5 
       //&& grid.pos(nodeNo, 1) >= 1//2
       //	 && grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 4//5
	 ){
      rho(0, nodeNo) += rho(1, nodeNo);
      rho(1, nodeNo) = 0.0;
    }
  }
   

  
  
  //                              Initiate frac geo from file
  //------------------------------------------------------------------------------------- Initiate frac geo from file
  setScalarAttribute(height, "frac_height_", vtklb);
  setScalarAttribute(tmpGradHeight, "gradH_", vtklb);
  for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {	
    gradHeight(0, 0, nodeNo) = tmpGradHeight(0, nodeNo); 
    gradHeight(0, 1, nodeNo) = tmpGradHeight(1, nodeNo);
  }
  /*
  for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {	
    if(grid.pos(nodeNo, 0)<=3 ){//|| grid.pos(nodeNo, 0) >= vtklb.getGlobaDimensions(0)-4){
      height(0, nodeNo) = 1;
      gradHeight.set(0, nodeNo) = 0;
    }
  }
  */
  
  VectorField<LT> ForceField2D(5, grid.size());
  
    
  //                                    Wall wettability
  //------------------------------------------------------------------------------------- Wall wettability
  //setScalarAttributeWall(rhoRel, "init_rho_wall_", vtklb, nodes);
  setScalarAttributeWall(phi, "init_rho_wall_", vtklb, nodes);
  
  //                                   Scale wettability
  //------------------------------------------------------------------------------------- scale wettability
  //normelizeScalarField(rhoRel, solidBoundaryNodes);
  normelizeScalarField(phi, solidBoundaryNodes);
    
    
    
  //                                       Velocity
  //------------------------------------------------------------------------------------- Velocity
  VectorField<LT> vel(1, grid.size());
  //                                   Initiate velocity
  //------------------------------------------------------------------------------------- Initiate velocity
  //std::srand(std::time(nullptr));
  std::srand(0);
  for (auto nodeNo: bulkNodes) {
    
    for (int dim=0; dim<LT::nD; ++dim){
      lbBase_t velNoise = 1e-6*((2.0*std::rand())/lbBase_t(RAND_MAX) - 1);
      //std::cout << nodeNo <<",  "<< velNoise << std::endl;
      vel(0, dim, nodeNo) = 0 /*+ velNoise*/;
    }
    
    //vel.set(0, nodeNo) = 0;
  }
    

  //                                    Phase separation
  //                             Color gradients for fluid pair
  //                     Absolute values of color gradients for fluid pair
  //------------------------------------------------------------------------------------- Phase separation
  ScalarField FNorm(nFluidFields*(nFluidFields-1)/2, grid.size());
  VectorField<LT> F(nFluidFields*(nFluidFields-1)/2, grid.size());
  ScalarField divF(nFluidFields*(nFluidFields-1)/2, grid.size());
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
  ScalarField EffRadiusCoefInv(nFluidFields*(nFluidFields-1)/2, grid.size());
  ScalarField normalPlaneAngleTop1(nFluidFields*(nFluidFields-1)/2, grid.size());
  
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
  
  
  //  for (auto nodeNo: bulkNodes) {
      /*for (int fieldNo=0; fieldNo < rhoD.num_fields(); ++fieldNo) {
	rhoD(fieldNo, nodeNo) = 0.01*rho(0, nodeNo);
	}*/
      //rhoD(1, nodeNo) = 0.0*rho(0, nodeNo);
      /*if(grid.pos(nodeNo, 0) > 0.4*(vtklb.getGlobaDimensions(0)-2))
	rhoD(1, nodeNo) = 0.0*rho(0, nodeNo);
      */
  //  }
      

  
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
  //omp_set_dynamic(0);     // Explicitly disable dynamic teams
  //omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
  //std::cout<<"threads="<<omp_get_num_threads()<<std::endl;

  omp_set_num_threads(2);

  std::cout << "before parallel section: " << std::endl;
  std::cout << "Num threads = " << omp_get_num_threads() << std::endl;
  std::cout << "Max threads = " << omp_get_max_threads() << std::endl;
  
  
  for (int fieldNo=0; fieldNo < f.num_fields(); ++fieldNo) {
    
    #pragma omp parallel for num_threads(2)  
    for (auto nodeNo: bulkNodes) {
      
      auto u2 = LT::dot(vel(0, nodeNo), vel(0, nodeNo));
      auto cu = LT::cDotAll(vel(0, nodeNo));
      f.set(fieldNo, nodeNo) = calcfeq<LT>(rho(fieldNo, nodeNo), u2, cu);
      fTmp.set(fieldNo, nodeNo) = 0;
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
  #pragma omp parallel for num_threads(2)
  for (auto nodeNo: bulkNodes) {
    rhoTot(0, nodeNo) = 0.0;
    for (int fieldNo=0; fieldNo < rho.num_fields(); ++fieldNo) {
      rhoTot(0, nodeNo) += rho(fieldNo, nodeNo);
    }
    
    
    for (int fieldNo=0; fieldNo < rho.num_fields(); ++fieldNo) {
      phi(fieldNo, nodeNo) = rho(fieldNo, nodeNo)/rhoTot(0, nodeNo);
    }
    
    for (int fieldNo=0; fieldNo < rhoD.num_fields(); ++fieldNo) {
      phiD(fieldNo, nodeNo) = rhoD(fieldNo, nodeNo)/rhoTot(0, nodeNo);
    }
    
    auto u2 = LT::dot(vel(0, nodeNo), vel(0, nodeNo));
    auto cu = LT::cDotAll(vel(0, nodeNo));
    fTot.set(0, nodeNo) = 0.0;//calcfeq<LT>(rhoTot(0, nodeNo), u2, cu);
    fTotTmp.set(0, nodeNo) = 0.0;
    
    /*
    for (int q = 0; q < LT::nQ; ++q) {
      fTot(0, q, nodeNo) = LT::w[q]*rhoTot(0, nodeNo);
      fTotTmp(0, q, nodeNo) = 0;
    }
    */
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



 //                            Set unitNormal for boundaries
 //------------------------------------------------------------------------------------- Set unitNormal for boundaries
  #pragma omp parallel for num_threads(2)
  for (auto nodeNo:  solidBoundaryNodes){
    
    if(grid.pos(nodeNo, 0) == 0){
      //std::cout<<" Left solid "<<nodes.getType(nodeNo)<<" ";
      for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	unitNormal.set(cnt, nodeNo)= 0;
	unitNormal(cnt, 0, nodeNo)= 1;	
      }
    }
    
    /*
    if(grid.pos(nodeNo, 1) == 0){
      //std::cout<<" Bottom solid "<<nodes.getType(nodeNo)<<" ";
      for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	unitNormal.set(cnt, nodeNo)= 0;
	unitNormal(cnt, 1, nodeNo)= -1;
      }
      
      if(grid.pos(nodeNo, 0) <= inletXEnd ){
	//std::cout<<" Bottom solid inlet "<<nodes.getType(nodeNo)<<" ";
	phi(0, nodeNo) = 1.0;
	phi(1, nodeNo) = 0.0;
	for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	  unitNormal.set(cnt, nodeNo) = 0;
	  unitNormal(cnt, 0, nodeNo) = 1;
	}
      }
      
    }
    */

    /*
    if(grid.pos(nodeNo, 1) == vtklb.getGlobaDimensions(1) - 3){
      //std::cout<<" Top solid "<<nodes.getType(nodeNo)<<" ";
      for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	unitNormal.set(cnt, nodeNo)= 0;
	unitNormal(cnt, 1, nodeNo)= 1;
      }
      
      if(grid.pos(nodeNo, 0) <= inletXEnd ){
	//std::cout<<" Top solid inlet "<<nodes.getType(nodeNo)<<" ";
	phi(0, nodeNo) = 1.0;
	phi(1, nodeNo) = 0.0;
	for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	  unitNormal.set(cnt, nodeNo) = 0;
	  unitNormal(cnt, 0, nodeNo) = 1;
	}
      }
    }
    */
    
  }
  
  
#pragma omp parallel for num_threads(2)
  for (auto nodeNo:  OutletBoundaryNodes){
    //std::cout<<"Type "<< nodes.getType(nodeNo) <<", x= "<< grid.pos(nodeNo, 0) << " ";
    for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
      for(int neighDir=0; neighDir<2; ++neighDir){
	  unitNormal.set(cnt, grid.neighbor(neighDir, nodeNo))= 0;
	  unitNormal(cnt, 0, grid.neighbor(neighDir, nodeNo))= 1;	  
	}
	unitNormal.set(cnt, grid.neighbor(7, nodeNo))= 0;
	unitNormal(cnt, 0, grid.neighbor(7, nodeNo))= 1;
    }
  }
  
  
  /*
  for (auto nodeNo: bulkNodes) {
  
    for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
      
      if(grid.pos(nodeNo, 0)== 1){
	for(int neighDir=3; neighDir<6; ++neighDir){
	  unitNormal.set(cnt, grid.neighbor(neighDir, nodeNo))= 0;
	  unitNormal(cnt, 0, grid.neighbor(neighDir, nodeNo))= 1;
	}
      }
      
      //if(grid.pos(nodeNo, 0)== vtklb.getGlobaDimensions(0) - 3){
      //for(int neighDir=0; neighDir<2; ++neighDir){
      //  unitNormal.set(cnt, grid.neighbor(neighDir, nodeNo))= 0;
      //  unitNormal(cnt, 0, grid.neighbor(neighDir, nodeNo))= 1;
      //  
      //}
      //unitNormal.set(cnt, grid.neighbor(7, nodeNo))= 0;
      //unitNormal(cnt, 0, grid.neighbor(7, nodeNo))= 1;
      //}
      
      
      if(grid.pos(nodeNo, 1) == 1){
	std::cout<<"Type "<<nodes.getType(nodeNo)<<" ";
	for(int neighDir=5; neighDir<8; ++neighDir){
	  unitNormal.set(cnt, grid.neighbor(neighDir, nodeNo))= 0;
	  unitNormal(cnt, 1, grid.neighbor(neighDir, nodeNo))= -1;
	}
	if(grid.pos(nodeNo, 0) >= 1 && grid.pos(nodeNo, 0) <= inletXEnd -1 ){
	  for(int neighDir=5; neighDir<8; ++neighDir){
	    unitNormal.set(cnt, grid.neighbor(neighDir, nodeNo)) = 0;
	    unitNormal(cnt, 0, grid.neighbor(neighDir, nodeNo)) = 1;
	  }
	}
      }
      if(grid.pos(nodeNo, 1) == vtklb.getGlobaDimensions(1) - 4){
	for(int neighDir=1; neighDir<4; ++neighDir){
	  unitNormal.set(cnt, grid.neighbor(neighDir, nodeNo))= 0;
	  unitNormal(cnt, 1, grid.neighbor(neighDir, nodeNo))= 1;
	}
	if(grid.pos(nodeNo, 0) >= 1 && grid.pos(nodeNo, 0) <= inletXEnd -1 ){
	  for(int neighDir=1; neighDir<4; ++neighDir){
	    unitNormal.set(cnt, grid.neighbor(neighDir, nodeNo)) = 0;
	    unitNormal(cnt, 0, grid.neighbor(neighDir, nodeNo)) = 1;
	  }
	}
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
  output.add_scalar_variables({"rhoTot", "rho", "rhoD", "phi", "divF", "kappa", "kappa2", "R",     "Q",     "frac_height", "grad_height", "cosTheta", "EffRCurveCoefInv", "normalPlaneAngleTop"}, 
			      { rhoTot,   rho,   rhoD,   phi,   divF,   kappa,    kappa2,    Rfield,  Qfield,  height,        tmpGradHeight, cosAng,     EffRadiusCoefInv,   normalPlaneAngleTop1});
  output.add_vector_variables({"vel", "F", "unitNormal", "forceField", "gradHeight", "forceField2D_"}, 
			      { vel,   F,   unitNormal,   ForceField,   gradHeight,   ForceField2D});

  if (myRank==0) {
    system(("cp ./input/input.dat "+outputDir2).c_str());
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

    #pragma omp parallel for num_threads(2)
    for (auto nodeNo: bulkNodes) {
      fTot.set(0, nodeNo) = f(0, nodeNo) + f(1, nodeNo);
      /*
      if( grid.pos(nodeNo, 0) <= inletXEnd
	 && grid.pos(nodeNo, 1) >= 4
	 && grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 7
	 ){
	f.set(0, nodeNo) = fTot(0, nodeNo);
	f.set(1, nodeNo) = 0.0;
      }
      */
    }
    
    //           Macroscopic values : rho, rhoTot and rhoRel = rho/rhoTot
    //------------------------------------------------------------------------------------- Macroscopic values : rho, rhoTot and rhoRel = rho/rhoTot
    calcDensityFields(rho, rhoRel, rhoTot, rhoD, phi, phiD, bulkNodes, f, fTot, g);

    
    mpiBoundary.communciateScalarField(rhoRel);
    mpiBoundary.communciateScalarField(phi);
    //mpiBoundary.communciateScalarField(phiD);

    int rampTimesteps = 5000;
    
    const lbBase_t ramp{ 0.5 * (1-std::cos(PI*std::min(i, rampTimesteps)/rampTimesteps)) };
    
    lbBase_t inletNodeNo = 0;
    lbBase_t inletDensityAve = 0.0;
    lbBase_t inletDensityAveGlobal;
    LbField<LT> sumFC_aveInlet(1,1);
    sumFC_aveInlet.set(0 ,0) = 0;
    LbField<LT> sumFC_aveInletGlobal(1,1);
    sumFC_aveInletGlobal.set(0 ,0) = 0;


    #pragma omp parallel for num_threads(2)
    for (auto nodeNo: bulkNodes) {
      // Cacluate gradient
      //------------------------------------------------------------------------------------- 
      
      //phi(1, nodeNo) = 1-phi(0, nodeNo);

      
            




      ScalarField phiNode(1, nFluidFields);
      for (int fieldNo = 0; fieldNo < nFluidFields; ++fieldNo){
	
	phiNode(0, fieldNo) = phi(fieldNo, nodeNo);
      }

      
      
      
      const CGAttributes<LT> cgat(nFluidFields, nodeNo, cNormInv, Gamma0, phiNode, phi, grid);
     
      cosAng(0, nodeNo) = 0.0;
      EffRadiusCoefInv(0, nodeNo) = 0.0;
      normalPlaneAngleTop1(0, nodeNo) = 0.0;
      //-----------------------------------------------------------------------------------------------------------------------------------
      //Stored interface normals are the lower triangular part of the interface normal matrix
      //and point from phase of lower phase index toward phase of higher phase index. e.g., 0->1, 0->2, 1->2 etc.
      //-----------------------------------------------------------------------------------------------------------------------------------
      
      for (int cnt=0; cnt<(nFluidFields*(nFluidFields-1)/2); ++cnt){
	FNorm(cnt, nodeNo)= cgat.FNorm_(0, cnt);
	F.set(cnt, nodeNo)= cgat.F_(0, cnt);
	//if (FNorm(cnt, nodeNo)>lbBaseEps){
	//  unitNormal.set(cnt, nodeNo)= F(cnt, nodeNo)/FNorm(cnt, nodeNo);
	//}
	unitNormal.set(cnt, nodeNo)= F(cnt, nodeNo)/(FNorm(cnt, nodeNo)+(FNorm(cnt, nodeNo)<lbBaseEps));
	divF(cnt, nodeNo)= cgat.divF_(0, cnt);
	
	/*
	if(phiNode(0, 0)*phiNode(0, 1)<1e-4){
	  FNorm(cnt, nodeNo) = 0;
	  F.set(cnt, nodeNo) = 0;
	  unitNormal.set(cnt, nodeNo)= 0;
	}
	*/


	/*
	if(grid.pos(nodeNo, 0)>= vtklb.getGlobaDimensions(0) - 2){
	  std::cout << "on boundary" << " ";
	    
	  for (int dim=0; dim<LT::nD; ++dim){
	    F(cnt, dim, nodeNo) = 2*F(cnt, dim, grid.neighbor(4, nodeNo)) - F(cnt, dim, grid.neighbor(4, grid.neighbor(4, nodeNo)));
	  }
	  if(grid.pos(nodeNo, 1) < 2 || grid.pos(nodeNo, 1) > vtklb.getGlobaDimensions(1) - 5){
	    F.set(cnt, nodeNo) = 0;
	  }
	  //F(cnt, 1, nodeNo) = F(cnt, 1, grid.neighbor(4, nodeNo)) + (F(cnt, 1, grid.neighbor(3, nodeNo)) -  F(cnt, 1, grid.neighbor(5, nodeNo)))*0.5;
	  FNorm(cnt, nodeNo) = sqrt(LT::dot(F(cnt, nodeNo),F(cnt, nodeNo)));
	  unitNormal.set(cnt, nodeNo) = F(cnt, nodeNo)/(FNorm(cnt, nodeNo)+(FNorm(cnt, nodeNo)<lbBaseEps));
	  
	}
	*/
	
	/*
	if(grid.pos(nodeNo, 0) >= 1 && grid.pos(nodeNo, 0) <= inletXEnd -1 ){ 
	  unitNormal.set(cnt, nodeNo) = 0;
	  unitNormal(cnt, 0, nodeNo) = 1;
	  
	}
	*/
	

	
      }

      
      
      

      if (grid.pos(nodeNo, 0) > inletXEnd){

	lbBase_t contactAngleTop = /*0.0;*/PI;/*0.5*PI;*/
	lbBase_t contactAngleBttm = /*0.0;*/PI;/*0.5*PI;*/
      
	const lbBase_t gradSquared = LT::dot(gradHeight(0,nodeNo),gradHeight(0,nodeNo));
	lbBase_t tmp = 1/sqrt(gradSquared+1);	
	lbBase_t tmp2 = LT::dot(F(0, nodeNo),gradHeight(0,nodeNo));
	lbBase_t tmp3 = sqrt(tmp2*tmp2);
	lbBase_t normalPlaneAngleTop = 0;
	if(tmp3>lbBaseEps) 
	  normalPlaneAngleTop = std::acos(tmp)*tmp2/(tmp3+ (tmp3<lbBaseEps));
	lbBase_t normalPlaneAngleBttm = 0;
	
	lbBase_t effContactAngleTop = contactAngleTop + normalPlaneAngleTop;
	/*
	lbBase_t effRadiusCoefTop;
	
	if(effContactAngleTop > PI || effContactAngleTop < 0){
	  effRadiusCoefTop = tmp;
	}
	else
	  effRadiusCoefTop = 1/cos(effContactAngleTop);
	*/  
	lbBase_t effContactAngleBttm = contactAngleBttm + normalPlaneAngleBttm;
	/*
	lbBase_t effRadiusCoefBttm;
	
	if(effContactAngleBttm > PI || effContactAngleBttm < 0){
	  effRadiusCoefBttm = tmp;
	}
	else
	  effRadiusCoefBttm = 1/cos(effContactAngleBttm);
	*/
	
	lbBase_t effContactAngle = 0.5*(effContactAngleTop + effContactAngleBttm);

	/*
	if(effContactAngle > PI) effContactAngle = PI;
	if(effContactAngle < 0) effContactAngle = 0;
	*/
	
	cosAng(0, nodeNo) = std::cos(effContactAngle);
	normalPlaneAngleTop1(0, nodeNo) = normalPlaneAngleTop;
	EffRadiusCoefInv(0, nodeNo) = std::cos(effContactAngle);
	//EffRadiusCoefInv(0, nodeNo) = 2/(effRadiusCoefTop + effRadiusCoefBttm);
      }
      
      //cosAng(1, nodeNo) = 0;
      //cosAng(2, nodeNo) = 1;

      /*
      //Calculate sum f and sum fc average at inlet
      if((grid.pos(nodeNo, 0) >= 2 && grid.pos(nodeNo, 0) <= 4)
	 && grid.pos(nodeNo, 1) >= 2
	 && grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 5
	 ){
     
	inletNodeNo++;
	inletDensityAve += LT::qSum(fTot(0, nodeNo));
	sumFC_aveInlet.set(0 ,0) += LT::qSumC(fTot(0, nodeNo));
      
      }
      */
    }
    /*
    if(inletNodeNo != 0){
      inletDensityAve /= inletNodeNo;
      sumFC_aveInlet.set(0 ,0) = sumFC_aveInlet(0 ,0)/ inletNodeNo;
    }
    //Communicate average inelet sum f and sum fc to all mpi nodes
    MPI_Allreduce(&inletDensityAve, &inletDensityAveGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int dim = 0; dim<LT::nD; dim++){
      lbBase_t tmp = sumFC_aveInlet(0, dim, 0);
      lbBase_t tmpGlobal;
      MPI_Allreduce(&tmp, &tmpGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      sumFC_aveInletGlobal(0, dim, 0) = tmpGlobal;
    }
    */

    
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

    /*
    int centerSrcPointInlet1 = 5;
    int centerSrcPointOutlet1 = vtklb.getGlobaDimensions(0)-8;
    lbBase_t inletdens = 1.0;
    
    lbBase_t outletdens = 1.0;
    lbBase_t epsilonDelta = 2.5;
    ScalarField deltaSrc(1, grid.size());
    ScalarField deltaSrc2(1, grid.size());
    
    for (auto nodeNo: bulkNodes) {
      if(grid.pos(nodeNo, 0) == centerSrcPointInlet1){
	inletdens = rhoTot(0, nodeNo);
	break;
      }
      if(grid.pos(nodeNo, 0) == centerSrcPointOutlet1){
	outletdens = rhoTot(0, nodeNo);
	break;
      }
      
    }
    
    for (auto nodeNo: bulkNodes) {
      lbBase_t xdiff = grid.pos(nodeNo, 0)-centerSrcPointInlet1;
      lbBase_t xdiff2 = grid.pos(nodeNo, 0)-centerSrcPointOutlet1;
      if(xdiff*xdiff<=epsilonDelta*epsilonDelta){
	deltaSrc(0, nodeNo) = 1/(2*epsilonDelta)*(1+std::cos(3.14159*xdiff/epsilonDelta));
      }
      else
	deltaSrc(0, nodeNo) = 0.0;

      if(xdiff2*xdiff2<=epsilonDelta*epsilonDelta){
	deltaSrc2(0, nodeNo) = 1/(2*epsilonDelta)*(1+std::cos(3.14159*xdiff2/epsilonDelta));
      }
      else
	deltaSrc2(0, nodeNo) = 0.0;
      
    }
    */
    
    //Main calculation loop
    //------------------------------------------------------------------------------------- 

    #pragma omp parallel for num_threads(2)
    for (auto nodeNo: bulkNodes) {
      
      Qfield(0, nodeNo)=0.0;


      
      
	
      
      
      // Cacluate gradient
      //------------------------------------------------------------------------------------- 
      
      ScalarField phiNode(1, nFluidFields);
      for (int fieldNo = 0; fieldNo < nFluidFields; ++fieldNo){
	Rfield(fieldNo, nodeNo)=0.0;
	phiNode(0, fieldNo) = phi(fieldNo, nodeNo);
      }
      //phi(1,nodeNo) = 1- phi(0, nodeNo);
 
      

      
      

      int centerSrcPointInlet = 5;
      lbBase_t ux_mean = 1e-4;
      lbBase_t QInletMean = meanInletFlux;//2.0e-3;/*FUNKER MED velThreshold=9e-2*/  //5e-4;//1e-4; 
      lbBase_t channelHeight = 13;
      lbBase_t uProfilePreFactor=12*ux_mean/(channelHeight*channelHeight);
      LbField<LT> uInletMean(1,1);
      uInletMean.set(0 ,0) = 0;
      uInletMean(0, 0 ,0) = QInletMean;
            
      ForceField.set(0, nodeNo)=0.0;

    
      
      if( (grid.pos(nodeNo, 0) >= 2 && grid.pos(nodeNo, 0) <= inletXEnd )
	  /*&& grid.pos(nodeNo, 1) >= 5
	    && grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 8*/

	  //&& grid.pos(nodeNo, 1)<= height(0, nodeNo)-1
	 ){


	
	Qfield(0, nodeNo) = QInletMean*rhoTot(0, nodeNo)*ramp;
	//Qfield(0, nodeNo) = QInletMean*ramp;//*12/(channelHeight*channelHeight)*(channelHeight*yGlobal-yGlobal*yGlobal);
	//uInletMean(0, 0 ,0) = uProfilePreFactor*(channelHeight*yGlobal-yGlobal*yGlobal)*ramp;
	//uInletMean(0, 0 ,0) = Qfield(0, nodeNo);

	//Rfield(0, nodeNo) = 2*(rhoTot(0, nodeNo) +0.5*Qfield(0, nodeNo) - rho(0, nodeNo));
	//Rfield(0, nodeNo) = Qfield(0, nodeNo);
	//Rfield(1, nodeNo) = 2*(0 - rho(1, nodeNo));//Qfield(0, nodeNo);//2*(rhoTot(0, nodeNo)+0.5*Q(0, nodeNo) - rho(1, nodeNo));
	//Rfield(1, nodeNo) = -Qfield(0, nodeNo);
	//Rfield(2, nodeNo) = 2*(0 - rho(2, nodeNo));

	
	//ForceField.set(0, nodeNo) = 2*(uInletMean(0, 0)*LT::qSum(fTot(0, nodeNo)) - calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo))));

	//ForceField.set(0, nodeNo) = 2*(uInletMean(0, 0) - calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo))));

	//ForceField.set(0, nodeNo) = 2*(uInletMean(0, 0) * inletDensityAveGlobal - sumFC_aveInletGlobal(0 ,0));

	//ForceField(0, 1, nodeNo) = 0.0;
	
	
	

      }

      /*
      if(grid.pos(nodeNo, 0) == 1
	 || (grid.pos(nodeNo, 1) < 2 && grid.pos(nodeNo, 0) <= inletXEnd -1)
	 || (grid.pos(nodeNo, 1) > vtklb.getGlobaDimensions(1) - 5 && grid.pos(nodeNo, 0) <= inletXEnd -1)){
	for (int fluidPairNo=0; fluidPairNo<(nFluidFields*(nFluidFields-1)/2); ++fluidPairNo){
	  FNorm(fluidPairNo, nodeNo) = 0;
	  F.set(fluidPairNo, nodeNo) = 0;
	  unitNormal.set(fluidPairNo, nodeNo)= 0;
	  
	}	
      }
      */
      
      
      
      
      if( grid.pos(nodeNo, 0) <= inletXEnd 

	 /*
	 && grid.pos(nodeNo, 1) >= 2
	 && grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 5
	 */
	 
	 ){
	Rfield(0, nodeNo) = /*0.0;*/Qfield(0, nodeNo);
	Rfield(1, nodeNo) = /*Qfield(0, nodeNo);*/0.0;//2*(0.0 - rho(1, nodeNo));//0.0;//Qfield(0, nodeNo)*(1-phi(0, nodeNo));
      }
      
   
    
      
      
      
      const CGAttributes<LT> cgat(nFluidFields, nodeNo, cNormInv, Gamma0, phiNode, phi, grid);

      
      VectorField<LT> IFTforceNode(1,1);
      IFTforceNode.set(0 ,0) = 0;
      
      //ForceField(0, 0, nodeNo)=1e-3;
      for (int forceNo=0; forceNo < ForceField2D.num_fields(); ++forceNo) 
	ForceField2D.set(forceNo, nodeNo)=0.0;
      
      
      //                      Lower triangular traversing of kappa and IFT
      //------------------------------------------------------------------------------------- Lower triangular traversing of kappa and IFT
      int cnt = 0;
      
      for (int fieldNo_k=0; fieldNo_k<nFluidFields; ++fieldNo_k) {
	for (int fieldNo_l = 0; fieldNo_l < fieldNo_k; ++fieldNo_l) {
	  const int sigmaBeta_ind = fieldNo_k*nFluidFields + fieldNo_l;
	  
	  //if(divF(cnt, nodeNo)>lbBaseEps){
	  
	  
	  kappa(cnt, nodeNo) = - div_test2<LT>(unitNormal, cnt, nodeNo, grid) ;
	  //kappa2(cnt, nodeNo) = - (div_test2<LT>(F, cnt, nodeNo, grid)*FNorm(cnt, nodeNo) - LT::dot(F(cnt, nodeNo),grad<LT>(FNorm, cnt, nodeNo, grid)))/(FNorm(cnt, nodeNo)*FNorm(cnt, nodeNo));
	  kappa2(cnt, nodeNo) = - (divF(cnt, nodeNo) - LT::dot(unitNormal(cnt, nodeNo),grad<LT>(FNorm, cnt, nodeNo, grid)))/(FNorm(cnt, nodeNo)+(FNorm(cnt, nodeNo)<lbBaseEps));
	  
	  //}
	  /*
	  else{
	    kappa(cnt, nodeNo) = 0.0;
	    kappa2(cnt, nodeNo) = 0.0;
	  }
	  */

	  /*
	  if( grid.pos(nodeNo, 0) <= inletXEnd){
	    kappa(cnt, nodeNo) = 0.0;
	    kappa2(cnt, nodeNo) = 0.0;
	  }
	  */
	  
	  if(abs(kappa(cnt, nodeNo)) >1.8)
	    kappa(cnt, nodeNo) = 0;
	  
	  if(abs(kappa2(cnt, nodeNo)) >2.0)
	    kappa2(cnt, nodeNo) *= FNorm(cnt, nodeNo);
	  
	  if(grid.pos(nodeNo, 0)>= vtklb.getGlobaDimensions(0) - 5){
	    kappa(cnt, nodeNo) = 0.0;
	  }

	  if ((phi(0, nodeNo)<1e-6 || phi(0, nodeNo)> (1-1e-6) || abs(kappa(cnt, nodeNo)) > 1.0)
	      && FNorm(cnt, nodeNo)<1e-2
	      /*&& grid.pos(nodeNo, 1) >= 2
		&& grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 5*/ ){
	    //if(phi(0, nodeNo) >0.9 && phi(0, nodeNo)<1.0){
	    //FNorm(cnt, nodeNo) = 0.0;
	    //F.set(cnt, nodeNo) =FNorm(0, nodeNo)*F(cnt, nodeNo) ;
	    //unitNormal.set(cnt, nodeNo) = 0.0;
	      //Rfield(0, nodeNo) = 2*(rhoTot(0, nodeNo) - rho(0, nodeNo));
	    //kappa(cnt, nodeNo) = 0.0;
	      
	    /*}
	    else if(phi(0, nodeNo) <0.1 && phi(0, nodeNo) > 0.0){
	      //FNorm(cnt, nodeNo) = 0.0;
	      //F.set(cnt, nodeNo) = 0;
	      unitNormal.set(cnt, nodeNo)= 0;
	      //Rfield(0, nodeNo) = 2*(0.0 - rho(0, nodeNo));
	      kappa(cnt, nodeNo) =0.0;
	    }  
	   */ 
	  }
	  
	  
	  //if (sqrt(kappa(cnt, nodeNo)*kappa(cnt, nodeNo))>1.0){
	  /*
	  if (FNorm(cnt, nodeNo)<1e-2){
	    kappa(cnt, nodeNo) *= FNorm(cnt, nodeNo);
	    //FNorm(cnt, nodeNo) = 0.0;
	    //F.set(cnt, nodeNo) = 0.0;
	    //unitNormal.set(cnt, nodeNo) = 0.0;
	  }
	  */

	  
	  /*
	  if(grid.pos(nodeNo, 0)>= vtklb.getGlobaDimensions(0) - 3){
 	    kappa(cnt, nodeNo) = 2*kappa(cnt, grid.neighbor(4, nodeNo)) - kappa(cnt, grid.neighbor(4, grid.neighbor(4, nodeNo)));
	  }
	  */
	  
	  /*
	  if(grid.pos(nodeNo, 1) > vtklb.getGlobaDimensions(1) - 5){
	    kappa(cnt, nodeNo) = kappa(cnt, grid.neighbor(6, nodeNo));
	  }
	  else if(grid.pos(nodeNo, 1) < 2){
	    kappa(cnt, nodeNo) = kappa(cnt, grid.neighbor(2, nodeNo));
	  }
	  */
	  
	  /*
	  if(FNorm(cnt, nodeNo)<1e-5){
	    kappa2(cnt, nodeNo) = 0.0;
	    kappa(cnt, nodeNo) = 0.0;
	  }
	  */
	  /*
	  if(phi(0, nodeNo)*(1-phi(0, nodeNo))<1e-3){ 
	    kappa(cnt, nodeNo) = 0.0;
	    kappa2(cnt, nodeNo) = 0.0;
	  }
	  */
	  
	  //kappa(cnt, nodeNo) = kappa2(cnt, nodeNo);

	  lbBase_t IFT_threshold =   1.0; 
	  IFTforceNode.set(0 ,0) += 0.5*rhoTot(0, nodeNo)*sigma[sigmaBeta_ind]*kappa(cnt, nodeNo)*F(cnt, nodeNo)*IFT_threshold;
	  ForceField2D.set(3, nodeNo) = IFTforceNode(0 ,0);  
	    
	  if (/*grid.pos(nodeNo, 0) > inletXEnd &&*/ grid.pos(nodeNo, 0)< vtklb.getGlobaDimensions(0) -5//){
	      /*&& grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 5 && grid.pos(nodeNo, 1) >= 2*/){
	         //std::min(1e4*phi(0, nodeNo)*(1-phi(0, nodeNo)),1.0);
	    //if(phi(0, nodeNo)*(1-phi(0, nodeNo))>1e-4){
	      
	      //}
	    //Quasi-2D
	    
	      //if ( grid.pos(nodeNo, 0)< vtklb.getGlobaDimensions(0)  ){  
	    //if( grid.pos(nodeNo, 1) >= 2 && grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 5){
	      //ForceField2D.set(4, nodeNo) += 0.5*sigma[sigmaBeta_ind]*rhoTot(0, nodeNo)*F(cnt, nodeNo)*2*cosAng(cnt, nodeNo)/height(0, nodeNo); 
	      ForceField2D.set(4, nodeNo) += 0.5*sigma[sigmaBeta_ind]*rhoTot(0, nodeNo)*F(cnt, nodeNo)*EffRadiusCoefInv(0, nodeNo)*2/height(0, nodeNo);					  
	      IFTforceNode.set(0 ,0) += ForceField2D(4, nodeNo);
	      //}
	    
	     
	      //ForceField2D.set(3, nodeNo) = -0.0*0.5*sigma[sigmaBeta_ind]*FNorm(cnt, nodeNo)/height(0, nodeNo)
	      //	*(unitNormal(cnt, nodeNo)*LT::dot(unitNormal(cnt, nodeNo),gradHeight(0,nodeNo)) - gradHeight(0,nodeNo)); //Averaged Capillary tensor effect

	      //ForceField2D.set(3, nodeNo) = 0.0;

	    //IFTforceNode.set(0 ,0) += ForceField2D(3, nodeNo); //Averaged Capillary tensor effect
	      
	  }
	  
	  //Quasi-2D
	  
	    
	  cnt++;
	}
      }
      cnt = 0;
      
      ForceField.set(0, nodeNo) += IFTforceNode(0, 0)*ramp + bodyForce(0, 0);


      //                              Calculate local fluid viscosity
      //------------------------------------------------------------------------------------- Calculate local fluid viscosity
      lbBase_t visc_inv = 0.0;
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo){
	visc_inv += phiNode(0, fieldNo)*kin_visc_inv[fieldNo];
      }
      
      lbBase_t viscNode = 1/visc_inv;
      //lbBase_t viscNode = 1/kin_visc_inv[0];

      //Quasi-2D
      
      if ( grid.pos(nodeNo, 0) > inletXEnd && grid.pos(nodeNo, 0)< vtklb.getGlobaDimensions(0) -5
	   /*&& grid.pos(nodeNo, 1) >= 2
	     && grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 5*/
	  ){
      //if ( grid.pos(nodeNo, 0)< vtklb.getGlobaDimensions(0)  ){  
	std::valarray<lbBase_t> velNodeTmp = calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo)), ForceField(0, nodeNo));

	lbBase_t velFactorNode = (1+0.5*12*viscNode/(height(0, nodeNo)*height(0, nodeNo)));

	auto rhoTotNodeTmp = (rhoTot(0, nodeNo)+0.5*Qfield(0, nodeNo))*1/(1+0.5*LT::dot(LT::qSumC(fTot(0, nodeNo)),gradHeight(0,nodeNo))/(velFactorNode*height(0, nodeNo)));

	velNodeTmp *= 1/(rhoTotNodeTmp*velFactorNode);

	std::valarray<lbBase_t>  momTmp = 1/velFactorNode*(LT::qSumC(fTot(0, nodeNo))+0.5*ForceField(0, nodeNo));
	
	//ForceField2D.set(0, nodeNo) = -12*rhoTotNodeTmp*viscNode*velNodeTmp/(height(0, nodeNo)*height(0, nodeNo));
	//ForceField2D.set(0, nodeNo) = -12*viscNode/velFactorNode*(LT::qSumC(fTot(0, nodeNo))+0.5*ForceField(0, nodeNo))/(height(0, nodeNo)*height(0, nodeNo))*(1-0.5*FNorm(0, nodeNo));
	ForceField2D.set(0, nodeNo) = -12*viscNode/(height(0, nodeNo)*height(0, nodeNo))*momTmp;
	
	
	ForceField.set(0, nodeNo) += ForceField2D(0, nodeNo);


	ForceField2D.set(1, nodeNo) = 0.0;//-(rhoTotNodeTmp-1)*LT::c2*gradHeight(0,nodeNo)/height(0, nodeNo)*0;

	ForceField.set(0, nodeNo) += ForceField2D(1, nodeNo);


	
	//lbBase_t Q2DNode = -rhoTotNodeTmp*LT::dot(velNodeTmp,gradHeight(0,nodeNo))/height(0, nodeNo);
	lbBase_t Q2DNode;
	Q2DNode = - LT::dot(momTmp,gradHeight(0,nodeNo))/height(0, nodeNo);
	//if(FNorm(0, nodeNo)>1e-2) Q2DNode = 0.0;
	
	
	Qfield(0, nodeNo) += Q2DNode;
	//if(rho(0, nodeNo)>0.9*rhoTot(0, nodeNo))
	//  Rfield(0, nodeNo) = Qfield(0, nodeNo);
	Rfield(0, nodeNo) = Qfield(0, nodeNo)*phi(0, nodeNo);
	Rfield(1, nodeNo) = Qfield(0, nodeNo)*phi(1, nodeNo);

	/*
	Rfield(0, nodeNo) = 2*((rhoTot(0, nodeNo) +0.5*Qfield(0, nodeNo))*phi(0, nodeNo) - rho(0, nodeNo));
	Rfield(1, nodeNo) = 2*((rhoTot(0, nodeNo) +0.5*Qfield(0, nodeNo))*phi(1, nodeNo) - rho(1, nodeNo));
	*/
	
      }
      
      
      
      //Quasi-2D END

      
      //Outlet velocity control
      //if ( grid.pos(nodeNo, 0)>= vtklb.getGlobaDimensions(0) -4 ){
   

      lbBase_t phi0Tmp = rho(0, nodeNo)/rhoTot(0, nodeNo);

      
      
      rhoTot(0, nodeNo) += 0.5*Qfield(0, nodeNo);

      
      //if((rho(0, nodeNo) + 0.5*Rfield(0, nodeNo))/rhoTot(0, nodeNo) != phi0Tmp){

      /*
      if(abs((rho(0, nodeNo) + 0.5*Rfield(0, nodeNo)) - rhoTot(0, nodeNo)*phi(0, nodeNo)) > lbBaseEps){
	//Rfield(0, nodeNo) = 2*(rhoTot(0, nodeNo)*phi0Tmp - rho(0, nodeNo));
	std::cout<<"Not concentration conservation"<<std::endl;
      }
      */
      
      /*
      if(rho(0, nodeNo)+0.5*Rfield(0, nodeNo) >rhoTot(0, nodeNo)){
	Rfield(0, nodeNo) = 2*(rhoTot(0, nodeNo) - rho(0, nodeNo));
      }
      else if(rho(0, nodeNo)+0.5*Rfield(0, nodeNo) < 0.0){
	Rfield(0, nodeNo) = 2*(0.0 - rho(0, nodeNo));
      }
      */
      
      /*
      if (FNorm(0, nodeNo)<1e-2 //sqrt(kappa(cnt, nodeNo)*kappa(cnt, nodeNo)) > 0.75
	  && grid.pos(nodeNo, 1) >= 2
	  && grid.pos(nodeNo, 1) <= vtklb.getGlobaDimensions(1) - 5 ){
	if(phi0Tmp >0.9 && phi0Tmp<1.0){
	  //FNorm(cnt, nodeNo) = 0.0;
	  //F.set(cnt, nodeNo) = 0;
	  unitNormal.set(cnt, nodeNo)= 0;
	  //Rfield(0, nodeNo) = 2*(rhoTot(0, nodeNo) - rho(0, nodeNo));
	  kappa(cnt, nodeNo) =0.0;
	}
	else if(phi0Tmp <0.1 && phi0Tmp > 0.0){
	  //FNorm(cnt, nodeNo) = 0.0;
	  //F.set(cnt, nodeNo) = 0;
	  unitNormal.set(cnt, nodeNo)= 0;
	  //Rfield(0, nodeNo) = 2*(0.0 - rho(0, nodeNo));
	  kappa(cnt, nodeNo) =0.0;
	}  
      }
      */
     
      
      
      for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo){
	rho(fieldNo, nodeNo) += 0.5*Rfield(fieldNo, nodeNo);
      }
      
      rhoTot(0, nodeNo) = rho(0, nodeNo) + rho(1, nodeNo);
      
      
      
      
      //rho(1, nodeNo) = rhoTot(0, nodeNo) - rho(0, nodeNo);

      
      phi(0, nodeNo) = rho(0, nodeNo)/rhoTot(0, nodeNo);
      phi(1, nodeNo) = rho(1, nodeNo)/rhoTot(0, nodeNo);
      
      
      //---------------------Velocity Instability reduction-------------------------------------------
      
      auto velNodeTmp1 = calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo)), ForceField(0, nodeNo));
      lbBase_t velThresh = 9e-2;
      if(grid.pos(nodeNo, 0)>= vtklb.getGlobaDimensions(0) -5) velThresh = 5e-2;
      lbBase_t speedNodeTmp = sqrt( LT::dot(velNodeTmp1, velNodeTmp1)); 
      if (speedNodeTmp > velThresh){
	auto velUnitVectorTmp = velNodeTmp1/speedNodeTmp;
	ForceField.set(0, nodeNo) = 2*(velThresh*velUnitVectorTmp*rhoTot(0, nodeNo) - LT::qSumC(fTot(0, nodeNo)));
      }
      
      //----------------------------------------------------------------------------------------------

      //---------------------Velocity perpendicular to walls in y-direction forced to zero------------------------------------------- 

      /*
      if((grid.pos(nodeNo, 1) < 2 || grid.pos(nodeNo, 1) > vtklb.getGlobaDimensions(1) - 5) && grid.pos(nodeNo, 0) <= inletXEnd ){
	std::valarray<lbBase_t>  ForceTmp = 2*(0 - LT::qSumC(fTot(0, nodeNo)));
	ForceField(0, 1, nodeNo) = ForceTmp[1];
      }
      */
      
      //---------------------Velocity perpendicular to walls in y-direction forced to zero------------------------------------------- 
      
      if( grid.pos(nodeNo, 0) < 2 ){
	std::valarray<lbBase_t>  ForceTmp = 2*(0 - LT::qSumC(fTot(0, nodeNo)));
	//ForceField(0, 0, nodeNo) = ForceTmp[0];
	ForceField.set(0, nodeNo) = ForceTmp;
      }
      
      
      //---------------------Velocity perpendicular to outlet direction forced to zero------------------------------------------- 

      /*
      if(grid.pos(nodeNo, 0) > vtklb.getGlobaDimensions(0) - 5){
	std::valarray<lbBase_t>  ForceTmp = 2*(0 - LT::qSumC(fTot(0, nodeNo)));
	ForceField(0, 1, nodeNo) = ForceTmp[1];
      }
      */
      //----------------------------------------------------------------------------------------------
      
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

      
      //if(fluidIntIndNode!=0)
      //visc_inv = 1/0.16667;
      
      //const lbBase_t tauFlNode = LT::c2Inv/visc_inv + 0.5;
      const lbBase_t tauFlNode = LT::c2Inv*viscNode + 0.5;
      
      //auto velNodeTmp = calcVel<LT>(fTot(0, nodeNo), LT::qSum(fTot(0, nodeNo)), ForceField(0, nodeNo));

      //auto strainRateNode =  1/(2*rhoTot(0,nodeNo)*LT::c2*tauFlNode)*calcStrainRateTildeLowTri<LT>(fTot(0, nodeNo), rhoTot(0, nodeNo), velNodeTmp, ForceField(0, nodeNo), Qfield(0, nodeNo));
      VectorField<LT> vecTmpNode(1,1);
      vecTmpNode.set(0,0)  = 0.0;//2*rhoTot(0,nodeNo)*viscNode/height(0, nodeNo)*LT::contractionLowTriVec(strainRateNode, gradHeight(0,nodeNo))*0.0;

      //ForceField2D.set(2, nodeNo) = vecTmpNode(0,0);
      
      //ForceField.set(0, nodeNo) += ForceField2D(2, nodeNo);

      
      
      
      
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

      /*
      if(LT::dot(velNode, velNode)<1e-15)
	velNode = 0*velNode;
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
      
      auto fTotNode = fTot(0, nodeNo);
      const auto rhoTotNode = rhoTot(0, nodeNo);
      
      const auto fToteqNode = calcfeq<LT>(rhoTotNode, u2, cu);
      
      //LB Regularization
      //if(grid.pos(nodeNo, 0) > vtklb.getGlobaDimensions(0) - 5 || grid.pos(nodeNo, 1) < 2 || grid.pos(nodeNo, 1) > vtklb.getGlobaDimensions(1) - 5){//LB Regularization close to outlet
        const auto fTotNeqNode =  fTotNode - fToteqNode;
	const auto PiTotNeqLowTri = LT::qSumCCLowTri(fTotNeqNode);
	const auto M_iNeq = -0.5*LT::qSumC(fTotNeqNode);
	const auto MNeq = -0.5*LT::qSum(fTotNeqNode);
	fTotNode = calcRegDist<LT>(fToteqNode, MNeq, M_iNeq, PiTotNeqLowTri);
	//}
      
      //LB Regularization END
      
      //const auto omegaBGKTot = newtonian.omegaBGK(tauFlNode, fTotNode, rhoTotNode, velNode, u2, cu, ForceField(0, nodeNo), Qfield(0, nodeNo));

      const std::valarray<lbBase_t> omegaBGKTot = calcOmegaBGK<LT>(fTotNode, tauFlNode, rhoTotNode, u2, cu);
      
      
      //const auto trE_Node = newtonian.trE();
      //------------------------------------------------------------------------------------- 

      //const auto deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tauFlNode, cu, u2, Qfield(0, nodeNo));
      const std::valarray<lbBase_t> deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tauFlNode, cu, u2, Qfield(0, nodeNo));
      const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tauFlNode, cu, uF, cF);
      
      
      fTot.set(0, nodeNo) = fTotNode + deltaOmegaQ0 + deltaOmegaF + omegaBGKTot;
      
      
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
	  

	const auto rhoNode = rho(fieldNo, nodeNo);

	/*
	const auto deltaOmegaR1   = calcDeltaOmegaR<LT>(tauPhaseField, cu, Rfield(fieldNo, nodeNo));
	
	auto fNode = f(fieldNo, nodeNo);
	
	
	const auto feqNode = calcfeq<LT>(rhoNode, u2, cu);
	
	//LB Regularization
	//if(grid.pos(nodeNo, 0) > vtklb.getGlobaDimensions(0) - 5 || grid.pos(nodeNo, 1) < 2 || grid.pos(nodeNo, 1) > vtklb.getGlobaDimensions(1) - 5){ //LB Regularization close to outlet
	  const auto fNeqNode =  fNode - feqNode;
	  const auto PiNeqLowTri = LT::qSumCCLowTri(fNeqNode);
	  const auto M_iNeq = -0.5*LT::qSumC(fNeqNode);
	  const auto MNeq = -0.5*LT::qSum(fNeqNode);
	  fNode = calcRegDist<LT>(feqNode, MNeq, M_iNeq, PiNeqLowTri);  
	  //}
	
	//LB Regularization END
	
	const auto omegaBGK = calcOmegaBGK_TEST<LT>(fNode, feqNode, tauPhaseField);        
	*/
	
	//LbField<LT> deltaOmegaFDiff(1,1);
	//deltaOmegaFDiff.set(0 ,0) = calcDeltaOmegaFDiff<LT>(tauPhaseField, phi(fieldNo, nodeNo)/*rhoRel(fieldNo, nodeNo)*/, cu, uF, cF);  // LBcollision
	
	//                                   Recoloring step
        //------------------------------------------------------------------------------------- Recoloring step
	int field_k_ind = (fieldNo*(fieldNo-1))/2;
	deltaOmegaRC.set(0, fieldNo) = 0;
	
	for (int field_l = 0; field_l < fieldNo; ++field_l) {
	  const int ind_F = field_k_ind + field_l;
	  const int ind_sigmaBeta = fieldNo*nFluidFields + field_l;
	  const lbBase_t IFT_threshold = 1.0;//std::min(1e6*rhoRel(field_l, nodeNo)*rhoRel(fieldNo, nodeNo),1.0);
	  //const lbBase_t IFT_threshold = (FNorm(ind_F, nodeNo)>1e-3);
	  
	  //const auto cn = LT::cDotAll(F(ind_F, nodeNo)/(FNorm(ind_F, nodeNo)+(FNorm(ind_F, nodeNo)<lbBaseEps)));
	  const auto cn = LT::cDotAll(unitNormal(ind_F, nodeNo));

	  //if (1 || (grid.pos(nodeNo, 0) > inletXEnd && grid.pos(nodeNo, 0)< vtklb.getGlobaDimensions(0) -4 && FNorm(ind_F, nodeNo) > 1.0e-3))
	  //if(grid.pos(nodeNo, 1) > vtklb.getGlobaDimensions(1) - 5 || grid.pos(nodeNo, 1) < 2)
	  //deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, 2*sigma[ind_sigmaBeta]*rhoTot(0, nodeNo)*IFT_threshold, FNorm(ind_F, nodeNo), cn);
	  
	
	  
	  //deltaOmegaRC.set(0, fieldNo) += rhoNode*beta[ind_sigmaBeta]*phi(field_l, nodeNo)*cn*cNormInv;
	  
	  deltaOmegaRC.set(0, fieldNo) += rhoNode*beta[ind_sigmaBeta]*phi(field_l, nodeNo)*cn;
	  //deltaOmegaRC.set(0, fieldNo) += rhoTotNode*phi(fieldNo, nodeNo)*beta[ind_sigmaBeta]*phi(field_l, nodeNo)*cn;
	  //deltaOmegaRC.set(0, fieldNo) += rhoRelNode(0,fieldNo)*beta[ind_sigmaBeta]*rhoRel(field_l, nodeNo)*cn;
	  
	}
	
	for (int field_l = fieldNo + 1; field_l < nFluidFields; ++field_l) {
	  const int field_k_ind = (field_l*(field_l-1))/2;
	  const int ind_F =  field_k_ind + fieldNo;
	  const int ind_sigmaBeta = fieldNo*nFluidFields + field_l;
	  
	  //const auto cn = LT::cDotAll(F(ind_F, nodeNo)/(FNorm(ind_F, nodeNo)+(FNorm(ind_F, nodeNo)<lbBaseEps)));
	  const auto cn = LT::cDotAll(unitNormal(ind_F, nodeNo));
	  
	  
	  //deltaOmegaRC.set(0, fieldNo) -= rhoNode*beta[ind_sigmaBeta]*phi(field_l, nodeNo)*cn*cNormInv;
	  
	  deltaOmegaRC.set(0, fieldNo) -= rhoNode*beta[ind_sigmaBeta]*phi(field_l, nodeNo)*cn;
        
	  //deltaOmegaRC.set(0, fieldNo) -= rhoRelNode(0,fieldNo)*beta[ind_sigmaBeta]*rhoRel(field_l, nodeNo)*cn;  
	}
	    

		
	deltaOmegaRC.set(0, fieldNo) *= wAll;

	
	//if ( grid.pos(nodeNo, 0) >= vtklb.getGlobaDimensions(0) - 4 /*|| grid.pos(nodeNo, 0) < inletXEnd */ ){
	//  deltaOmegaRC.set(0, fieldNo) *= 0*wAll;
	//}
	

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
	//f.set(fieldNo, nodeNo) = fNode + omegaBGK + deltaOmegaRC(0, fieldNo) + deltaOmegaFDiff(0 ,0) + deltaOmegaR1;
	f.set(fieldNo, nodeNo) = fTot(0, nodeNo) *phi(fieldNo, nodeNo) + deltaOmegaRC(0, fieldNo) /*+ deltaOmegaR1*/;
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


      //fTotTmp.propagateTo(0, nodeNo, fTot(0, nodeNo) + deltaOmegaST(0, 0), grid);

      //fTotTmp.propagateTo(0, nodeNo, fTotNode + deltaOmegaQ0 + deltaOmegaF + omegaBGKTot + deltaOmegaST(0, 0), grid);
      
      
      

	
	
    }//------------------------------------------------------------------------------------- End nodes



    
    // Swap data_ from fTmp to f;
    //------------------------------------------------------------------------------------- 
    //fTot.swapData(fTotTmp);  // LBfield
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

    //zouHePressureBoundaryRight(OutletBoundaryNodes, fTot, 1.0, ForceField, grid);
    
    //zouHePressureBoundary(InletBoundaryNodes, fTot, 1.0+0.001*ramp, ForceField, grid);
    
    // Mpi
    //------------------------------------------------------------------------------------- 
    mpiBoundary.communicateLbField(f, grid);
    // Half way bounce back
    //------------------------------------------------------------------------------------- 
    bounceBackBnd.apply(f, grid);
    
    zouHeConcBoundaryRight(OutletBoundaryNodes, f, phi, 1.0, F, ForceField, grid);
    
    
    
   
    //=====================================================================================
    //
    //                                 WRITE TO FILE
    //
    //=====================================================================================
    if ( int(0.2*nItrWrite) > 0){
      if (((i % int(0.2*nItrWrite)) == 0) && myRank==0 ) {
      
	std::cout << "Iteration: " << i << std::endl;
      
      }
    }
    else if(myRank==0 ) {
      std::cout << "Iteration: " << i << std::endl;
    }
      
    if ( ((i % nItrWrite) == 0)  ) {
      
      output.write(i);
      if (myRank==0) {
	std::cout << "PLOT AT ITERATION: " << i << std::endl;
      }
    }
    
  } //-------------------------------------------------------------------------------------  End iterations
  
  MPI_Finalize();
  
  return 0;
}
