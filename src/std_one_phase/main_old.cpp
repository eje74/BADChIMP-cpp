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


//                                SET THE LATTICE TYPE
//------------------------------------------------------------------------------------- SET THE LATTICE TYPE
#define LT D2Q9
template <typename T>
using Output = LBOutputUnstructured<LT, T, VTK::BINARY, VTK::voxel>;
#define VTK_CELL VTK::pixel
//#define LT D3Q19
//#define VTK_CELL VTK::voxel
 

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

  lbBase_t tau = kin_visc_inv[0]*LT::c2Inv + 0.5;

  //                                     Body force
  //------------------------------------------------------------------------------------- Body force
  VectorField<LT> bodyForce(1, 1, input["fluid"]["bodyforce"]);
  
  //                          Parameters for compressibility
  //------------------------------------------------------------------------------------- Parameters for compressibility
  std::valarray<lbBase_t> alpha = LT::w0*inputAsValarray<lbBase_t>(input["fluid"]["alpha"]);
  std::valarray<lbBase_t> Gamma0 = inputAsValarray<lbBase_t>(input["fluid"]["alpha"]);
  std::valarray<lbBase_t> GammaNonZero = (1-LT::w0*Gamma0)/(1-LT::w0);

  //                                    Output directory number
  //------------------------------------------------------------------------------------- Output directory number
  std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));  
  std::string outputDir2 = outputDir + "/out" + dirNum;




  
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
  
  //                              Initiate density from file
  //------------------------------------------------------------------------------------- Initiate density from file
  vtklb.toAttribute("init_rho");
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
    rho(0, n) = vtklb.getScalarAttribute<lbBase_t>();
  }
    
 

  //                                      Sources
  //------------------------------------------------------------------------------------- Sources
  ScalarField Rfield(nFluidFields, grid.size());
  ScalarField Qfield(1, grid.size());
  
 
    
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
      vel(0, dim, nodeNo) = 0 /*+ velNoise*/;
    }
    
    //vel.set(0, nodeNo) = 0;
  }
    

  //                                       Force Field
  //------------------------------------------------------------------------------------- Force Field
  
  VectorField<LT> ForceField(1, grid.size());

  //                                   Initiate Force Field
  //------------------------------------------------------------------------------------- Initiate Force Field
  
  for (auto fieldNo=0; fieldNo < ForceField.num_fields(); ++fieldNo) {
    for (auto nodeNo: bulkNodes) {
      ForceField.set(fieldNo, nodeNo) = 0;
    }
  }
  
  //=====================================================================================
  //
  //                                SETUP BOUNDARY
  //
  //=====================================================================================

  HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

  //=====================================================================================
  //
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

 
    
  //=====================================================================================
  //
  //                                  OUTPUT VTK
  //
  //=====================================================================================
  Output<double> output(grid, bulkNodes, outputDir2, myRank, nProcs); 
  output.add_file("lb_run");
  output.add_scalar_variables({"rho", "Q"    }, 
			      { rho,   Qfield});
  output.add_vector_variables({"vel", "forceField"}, 
			      { vel,   ForceField });

  if (myRank==0) {
    system(("cp ./input/input.dat "+outputDir2).c_str());
  }

  
  //=====================================================================================
  //
  //                                  MAIN LOOP
  //
  //=====================================================================================

  for (int i = 0; i <= nIterations; i++) {

    
    int rampTimesteps = 5000; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////NB!
    
    const lbBase_t ramp{ 0.5 * (1-std::cos(PI*std::min(i, rampTimesteps)/rampTimesteps)) };
    
  
    
    //                                   Main calculation loop
    //------------------------------------------------------------------------------------- 


    for (auto nodeNo: bulkNodes) {


      Qfield(0, nodeNo)=0.0;
      ForceField.set(0, nodeNo)=0.0;   
      ForceField.set(0, nodeNo) += bodyForce(0, 0);

	 
      //                            Copy of local velocity diestirubtion
      //------------------------------------------------------------------------------------- Copy of local velocity diestirubtion
      const std::valarray<lbBase_t> fNode = f(0, nodeNo);
      
      //                                      Macroscopic values
      //------------------------------------------------------------------------------------- Macroscopic values
      const lbBase_t rhoNode = calcRho<LT>(fNode);
      const auto velNode = calcVel<LT>(fNode, rhoNode, bodyForce(0, 0));
      
      //                            Save density and velocity for printing
      //------------------------------------------------------------------------------------- Save density and velocity for printing
      rho(0, nodeNo) = rhoNode;
      vel.set(0, nodeNo) = velNode;
      
      //                                    BGK-collision term
      //------------------------------------------------------------------------------------- BGK-collision term
      
      
      const lbBase_t u2 = LT::dot(velNode, velNode);
      const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
      
      //const auto feqNode = calcfeq<LT>(rhoNode, u2, cu);
      //const auto omegaBGKTot = newtonian.omegaBGK(tau, fNode, rhoNode, velNode, u2, cu, ForceField(0, nodeNo), Qfield(0, nodeNo));
      
      const std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);
            
      //                           Calculate the Guo-force correction
      //------------------------------------------------------------------------------------- Calculate the Guo-force correction
      const lbBase_t uF = LT::dot(velNode, bodyForce(0, 0));
      const std::valarray<lbBase_t> cF = LT::cDotAll(bodyForce(0, 0));
      const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);
      
      //                           Calculate the mass-source correction
      //------------------------------------------------------------------------------------- Calculate the mass-source correction
      const std::valarray<lbBase_t> deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tau, cu, u2, Qfield(0, nodeNo));
      
      
      //                               Collision and propagation
      //------------------------------------------------------------------------------------- Collision and propagation
      fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF + deltaOmegaQ0, grid);
      
      
      
      
    }//------------------------------------------------------------------------------------- End nodes
    // Swap data_ from fTmp to f;
    //------------------------------------------------------------------------------------- 
    
    f.swapData(fTmp);  // LBfield
    
    
   
    //=====================================================================================
    //
    //                               BOUNDARY CONDITIONS
    //
    //=====================================================================================

  
    
    // Mpi
    //------------------------------------------------------------------------------------- 
    mpiBoundary.communicateLbField(f, grid);
    // Half way bounce back
    //------------------------------------------------------------------------------------- 
    bounceBackBnd.apply(f, grid);
    

    
    
    
   
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
