
// //////////////////////////////////////////////
//
// RELATIVE PERMEABILITY FROM LS SIMULATIONS
//   -  No co-flowing boundaries
//
// For documentation see:
//    doc/documentation.pdf
//
// //////////////////////////////////////////////

#include "../LBSOLVER.h"
#include "../IO.h"

//------------------------------------------------------------------------------------- SET THE LATTICE TYPE
#define LT D3Q19

//##################################################################################### Boundary conditions
//===================================================================================== Boundary setup
//------------------------------------------------------------------------------------- solid-fluid
//------------------------------------------------------------------------------------- pressure-fluid

//===================================================================================== Boundary conditions implementation
//------------------------------------------------------------------------------------- solid-fluid
//------------------------------------------------------------------------------------- pressure-fluid
/* Here we are using anti bounce back.
 * eq. 5.53 (p. 200) in Krugers book, but without velocity interpolation
 */

int main()
{
  //===================================================================================== Setup mpi
  MPI_Init(NULL, NULL);
  int nProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  //===================================================================================== Setup paths
  std::string chimpDir = "./../";
  std::string mpiDir = chimpDir + "input/mpi/";
  std::string inputDir = chimpDir + "input/";
  Input input(inputDir + "input.dat");
  std::string outputDir = chimpDir + "output/";

  //===================================================================================== Grid and Geometry setup
  LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
  Grid<LT> grid(vtklb);
  Nodes<LT> nodes(vtklb, grid);
  BndMpi<LT> mpiBoundary(vtklb, nodes, grid);
  // Set bulk nodes
  std::vector<int> bulkNodes = findBulkNodes(nodes);


  //===================================================================================== Read input-file
  //------------------------------------------------------------------------------------- Number of iterations
  int nIterations = input["iterations"]["max"];
  //------------------------------------------------------------------------------------- Write interval
  int nItrWrite = input["iterations"]["write"];
  //------------------------------------------------------------------------------------- Relaxation time
  lbBase_t tau = input["fluid"]["tau"];
  //------------------------------------------------------------------------------------- Body force
  VectorField<LT> bodyForce(1, 1, input["fluid"]["bodyforce"]);
  //------------------------------------------------------------------------------------- Write to screen
  if (myRank == 0)
  {
    std::cout << std::endl;
    std::cout << "INPUT PARAMETERS:" << std::endl;
    std::cout << "  nIterations = " << nIterations << std::endl;
    std::cout << "  nItrWrite = " << nItrWrite << std::endl;
    std::cout << "  tau = " << tau << std::endl;
    std::cout << "  bodyforce = " << "[";
    std::cout << bodyForce(0, 0, 0) << ", " << bodyForce(0, 1, 0) << ", " << bodyForce(0, 2, 0);
    std::cout << "]" << std::endl;
    std::cout << std::endl;
  }
  //===================================================================================== Macroscopic fields
  //------------------------------------------------------------------------------------- Density
  ScalarField rho(1, grid.size());
  //                              Initiate density
  //------------------------------------------------------------------------------------- Initiate density
  for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    rho(0, n) = 1.0;
  }
  //------------------------------------------------------------------------------------- Velocity
  VectorField<LT> vel(1, grid.size());
  //                                   Initiate velocity
  //------------------------------------------------------------------------------------- Initiate velocity
  for (auto nodeNo : bulkNodes)
  {
    vel.set(0, nodeNo) = 0;
  }

  //===================================================================================== Boundary setup
  //------------------------------------------------------------------------------------- solid-fluid
  HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);
  //------------------------------------------------------------------------------------- pressure-fluid

  //===================================================================================== Lb fields
  //------------------------------------------------------------------------------------- declaration
  LbField<LT> f(1, grid.size());
  LbField<LT> fTmp(1, grid.size());
  //                           Initiate lb distributions
  //------------------------------------------------------------------------------------- Initiate lb distributions
  for (auto nodeNo : bulkNodes)
  {
    auto u2 = LT::dot(vel(0, nodeNo), vel(0, nodeNo));
    auto cu = LT::cDotAll(vel(0, nodeNo));
    f.set(0, nodeNo) = calcfeq<LT>(rho(0, nodeNo), u2, cu);
    fTmp.set(0, nodeNo) = 0;
  }
  //===================================================================================== Output vtk
  Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs);
  output.add_file("lb_run");
  output.add_scalar_variables(
      {"rho"},
      {rho});
  output.add_vector_variables(
      {"vel"},
      {vel});

  // ###################################################################################### MAIN LOOP
  std::time_t start, end;

  if (myRank == 0)
  {
    std::time(&start);
  }

  for (int i = 0; i <= nIterations; i++)
  {
    //------------------------------------------------------------------------------------- Rampup
    //   int rampTimesteps = 5000; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////NB!
    //   const lbBase_t ramp{ 0.5 * (1-std::cos(PI*std::min(i, rampTimesteps)/rampTimesteps)) };
    //                                   Main calculation loop
    //-------------------------------------------------------------------------------------
    //===================================================================================== For each node
    for (auto nodeNo : bulkNodes)
    {
      //------------------------------------------------------------------------------------- Copy of local velocity distribution
      const std::valarray<lbBase_t> fNode = f(0, nodeNo);
      //                                      Macroscopic values
      //------------------------------------------------------------------------------------- Macroscopic values
      lbBase_t rhoNode = calcRho<LT>(fNode);
      const std::valarray<lbBase_t> forceNode = bodyForce(0, 0);
      const auto velNode = calcVel<LT>(fNode, rhoNode, forceNode);
      //                            Save density and velocity for printing
      //------------------------------------------------------------------------------------- Save density and velocity for printing
      rho(0, nodeNo) = rhoNode;
      vel.set(0, nodeNo) = velNode;

      //                                    BGK-collision term
      //------------------------------------------------------------------------------------- BGK-collision term
      const lbBase_t u2 = LT::dot(velNode, velNode);
      const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
      const std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);
      //                           Calculate the Guo-force correction
      //------------------------------------------------------------------------------------- Calculate the Guo-force correction
      const lbBase_t uF = LT::dot(velNode, forceNode);
      const std::valarray<lbBase_t> cF = LT::cDotAll(forceNode);
      const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);
      //------------------------------------------------------------------------------------- Collision and propagation
      fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);
      // fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF + deltaOmegaQ0, grid);
    } //------------------------------------------------------------------------------------- End for bulkNodes
    //------------------------------------------------------------------------------------- Swap lb fields
    f.swapData(fTmp);

    //=====================================================================================
    //
    //                               BOUNDARY CONDITIONS
    //
    //=====================================================================================
    //                                     MPI
    //------------------------------------------------------------------------------------- MPI
    mpiBoundary.communicateLbField(f, grid);
    // mpiBoundary.communciateVectorField_TEST(vel);
    //                             Solid-fluid boundary
    //------------------------------------------------------------------------------------- solid-fluid
    bounceBackBnd.apply(f, grid);
    //------------------------------------------------------------------------------------- pressure-fluid
    //------------------------------------------------------------------------------------- fluid-fluid

    //===================================================================================== Write to file
    if (((i % nItrWrite) == 0))
    {
      if (myRank == 0)
      {
        std::cout << "INTERATION:  " << i << std::endl;
        output.write(i);
      }
      // if (((i % (50 * nItrWrite)) == 0))
      // {
      //   output.write(i);
      // }
      // std::vector<lbBase_t> massFluxLocal(2, 0.0);
      // for (auto n : pressureFluidNodes)
      // {
      //   int fluidPhase = (nodes.getTag(n) & 3) - 1;
      //   if ((fluidPhase != 0) && (fluidPhase != 1))
      //   {
      //     std::cout << "Fluid phase = " << fluidPhase << " in write mass flux" << std::endl;
      //     MPI_Finalize();
      //     exit(1);
      //   }
      //   massFluxLocal[fluidPhase] += vel(0, 2, n) * rho(0, n);
      // }
      // std::vector<lbBase_t> massFluxGlobal(2, 0.0);
      // MPI_Allreduce(massFluxLocal.data(), massFluxGlobal.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      /*      if (myRank==0) {
              lbBase_t q1 = 0.5*massFluxGlobal[0];
              lbBase_t q2 = 0.5*massFluxGlobal[1];
              lbBase_t q1_change = (q1 - oldMassFlux[0])/(q1 + 1e-15);
              lbBase_t q2_change = (q2 - oldMassFlux[1])/(q2 + 1e-15);
        std::cout << "PLOT AT ITERATION: " << i << std::endl;
              std::cout << "q1 = " << q1 << " (" << q1_change << ")" << std::endl;
              std::cout << "q2 = " << q2 << " (" << q2_change << ")" << std::endl;
        std::ofstream myfile;
        myfile.open(nameOutputFolder + ".flux", std::ios::out | std::ios::app);
        myfile  << "PLOT AT ITERATION: " << i << "\n";
              myfile  << "q1 = " << q1 << " (" << q1_change << ")" << "\n";
              myfile  << "q2 = " << q2 << " (" << q2_change << ")" << "\n";
        myfile.close();

              oldMassFlux[0] = q1;
              oldMassFlux[1] = q2;
            } */
      //                           Update mass source
      //------------------------------------------------------------------------------------- Update mass source
      // MPI_Allreduce(massChangeLocal.data(), massChange.data(), massChange.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      // Sett local to zeros
      // std::fill(massChangeLocal.begin(), massChangeLocal.end(), 0.0);
    }
  } //----------------------------------------------------------------------------------------  End for nIterations
  if (myRank == 0)
  {
    std::time(&end);
    double runTime = double(end - start);
    std::cout << "Run time = " << runTime << "\n";
  }

  MPI_Finalize();

  return 0;
}
