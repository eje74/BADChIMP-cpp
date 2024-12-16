// ///////////////////////////////////////////////////////////
//
//  RANS solver
//
//  TO DO
//        - new viscosity model (new RANS class)
//            - at least TRT
//            - Other possibilites : Cumulant, Regularized LB..
//        - 2 scalar advective diffusive - fields k and \epsilon
//                -- with source terms
//
//        -Focus on bulk:
//           - First standard bounce back
//             Then:
//             - Boundary cond. logarithmic law of the wall
//
//        -Lattice Boltzmann model:
//            - See eg. Filippova ... and succi 2001
//            - Advection diffusion use LB algorithm
//
//
// ///////////////////////////////////////////////////////////

#include "../LBSOLVER.h"
#include "../IO.h"
#include "LBrans.h"
#include <chrono>
#include <numeric>
#include "LBBoundaryInterpolation.h"
#include "LBBoundaryRans.h"
#include "ransmain.h"

// SET THE LATTICE TYPE
//------------------------------------------------------------------------------------- SET THE LATTICE TYPE
#define LT D2Q9

//=======================================================================================
//
//    RUN Iteration
//
//=======================================================================================
//=====================================================================================
//
//                  R E G U L A R I Z E D   D I S T R I B U T I O N
//
//=====================================================================================
template<int DIM>
inline lbBase_t lowerDiagCcCont(const std::valarray<lbBase_t> &m, const std::valarray<lbBase_t> & c)
{
    std::cout << "Error in template specialization" << std::endl;
    exit(1);
    return 0;
}

template<>
inline lbBase_t lowerDiagCcCont<2>(const std::valarray<lbBase_t> &m, const std::valarray<lbBase_t> & c)
{
    return c[0]*c[0]*m[0] + 2*c[0]*c[1]*m[1] + c[1]*c[1]*m[2];
}

template<>
inline lbBase_t lowerDiagCcCont<3>(const std::valarray<lbBase_t> &m, const std::valarray<lbBase_t> & c)
{
    return c[0]*c[0]*m[0] + 2*c[1]*c[0]*m[1] + c[1]*c[1]*m[2] + 2*c[2]*c[0]*m[3] + 2*c[2]*c[1]*m[4] + c[2]*c[2]*m[5];
}

template<int DIM>
inline lbBase_t lowerDiagTrace(const std::valarray<lbBase_t> &m)
{
    std::cout << "Error in template specialization" << std::endl;
    exit(1);
    return 0;
}

template<>
inline lbBase_t lowerDiagTrace<2>(const std::valarray<lbBase_t> &m)
{
    return m[0] + m[2];   
}

template<>
inline lbBase_t lowerDiagTrace<3>(const std::valarray<lbBase_t> &m)
{
    return m[0] + m[2] + m[5];   
}

template<typename DXQY>
std::valarray<lbBase_t> fRegularized(const std::valarray<lbBase_t> &f, const lbBase_t dRhoNode)
{
    const lbBase_t M = DXQY::qSum(f) + dRhoNode;
    const auto Mi = DXQY::qSumC(f);
    const auto Mij = DXQY::qSumCCLowTri(f);

    std::valarray<lbBase_t> fReg(DXQY::nQ);
    const lbBase_t c2traceM = DXQY::c2*lowerDiagTrace<DXQY::nD>(Mij);
    for (int q=0; q<DXQY::nQ; ++q) {
        const lbBase_t Qdelta = DXQY::cNorm[q]*DXQY::cNorm[q] - DXQY::nD*DXQY::c2;
        const lbBase_t cM = DXQY::cDotRef(q, Mi)*DXQY::c2Inv;        
        const lbBase_t QM = DXQY::c4Inv0_5*(lowerDiagCcCont<DXQY::nD>(Mij, DXQY::cValarray(q)) - c2traceM - DXQY::c2*Qdelta*M);
        fReg[q] = DXQY::w[q]*(M + cM + QM);
    }

    return fReg;
}

//=======================================================================================
//
//                                      M A I N
//
//=======================================================================================
int main()
{
  //=====================================================================================
  //
  //                                   SETUP MPI
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

  std::string chimpDir = "../";
  // std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
  std::string mpiDir = chimpDir + "input/mpi/";
  std::string inputDir = chimpDir + "input/";
  std::string outputDir = chimpDir + "output/";

  //=====================================================================================
  //
  //                             SETUP GRID AND GEOMETRY
  //
  //=====================================================================================

  Input input(inputDir + "rans_input.dat");
  LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
  Grid<LT> grid(vtklb);
  Nodes<LT> nodes(vtklb, grid);
  BndMpi<LT> mpiBoundary(vtklb, nodes, grid);

  // Set bulk nodes
  //------------------------------------------------------------------------------------- Set bulk nodes
  std::vector<int> bulkNodes = findBulkNodes(nodes);

  std::vector<int> outletBoundaryNodes;
  std::vector<int> inletBoundaryNodes;

  for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo)
  {

    if (grid.pos(nodeNo, 0) == vtklb.getGlobaDimensions(0) - 3)
      outletBoundaryNodes.push_back(nodeNo);

    if (grid.pos(nodeNo, 0) == 0)
      inletBoundaryNodes.push_back(nodeNo);
  }

  for (auto nodeNo : outletBoundaryNodes)
  {
    std::cout << "Outlet; node no: " << nodeNo << std::endl;
  }
  for (auto nodeNo : inletBoundaryNodes)
  {
    std::cout << "Inlet; node no: " << nodeNo << std::endl;
  }

  //                                 Set solid boundary
  //------------------------------------------------------------------------------------- Set solid boundary
  // std::vector<int> solidBoundaryNodes = findSolidBndNodes(nodes);
  //------------------------------------------------------------------------------------- Fluid-Solid
  auto fluidSolidBoundaryNodes = findFluidBndNodes(nodes);
  Boundary<LT> boundary(fluidSolidBoundaryNodes, nodes, grid);
  std::vector<InterpolationElement> bndInterp(boundary.size());
  readBoundaryNodeFile(mpiDir + "boundary" + std::to_string(myRank) + ".txt", bndInterp, boundary, nodes, grid, vtklb);


  //=====================================================================================
  //
  //                                 READ FROM INPUT
  //
  //=====================================================================================

  //                               Number of iterations
  //------------------------------------------------------------------------------------- Number of iterations
  int nIterations = input["iterations"]["max"];

  // //                                  Write interval
  // //------------------------------------------------------------------------------------- Write interval
  int nItrWrite = input["iterations"]["write"];

  //int TRT = input["TRT"]["TRT_true"];

  //                                  Relaxation time
  //------------------------------------------------------------------------------------- Relaxation time                          // FJERNE?
  // lbBase_t tau = LT::c2Inv * input["fluid"]["viscosity"] + 0.5;

  //                                    Body force
  //------------------------------------------------------------------------------------- Body force
  VectorField<LT> bodyForce(1, 1);
  bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);
  // bodyForce.set(0, 0) = input["fluid"]["bodyforce"];

  // VectorField<LT> bodyForce(1, 1, input["fluid"]["bodyforce"]);

  lbBase_t u_ref = input["RANS"]["inlet"]["u_ref"];
  lbBase_t l_turb = input["RANS"]["inlet"]["l_turb"];
  lbBase_t I_turb = input["RANS"]["inlet"]["I_turb"];
  int rampTimesteps = input["RANS"]["inlet"]["rampTimeSteps"];
  lbBase_t kInlet = 1.0 * std::pow(I_turb * u_ref, 2);
  lbBase_t epsilonInlet = std::pow(input["RANS"]["k-epsilonCoef"]["C_mu"], 1.0) * std::pow(kInlet, 1.5) / l_turb;
  lbBase_t y_pluss_cut_off = input["RANS"]["wall"]["y_pluss_cutoff"];

  if (myRank == 0)
  {
    std::cout << "-------------RANS Inlet-----------------" << std::endl;
    std::cout << "kInlet = " << kInlet << std::endl;
    std::cout << "epsilonInlet = " << epsilonInlet << std::endl;
    std::cout << "-------------RANS wall -----------------" << std::endl;
    std::cout << "y_pluss_cutoff = " << y_pluss_cut_off << std::endl;
  }

  //                                    Output directory number
  //------------------------------------------------------------------------------------- Output directory number
  std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));
  std::string outputDir2 = outputDir + "out" + dirNum;

  Output<LT> writeboundary(grid, bulkNodes, outputDir2, myRank, nProcs);
  writeboundary.add_file("boundary");
  ScalarField gamma(1, grid.size());
  VectorField<LT> normals(1, grid.size());
  for (auto &bn: bndInterp) {
    const int nodeNo = bn.nodeNo;
    gamma(0, nodeNo) = bn.gamma;
    normals(0, 0, nodeNo) = bn.normal[0];
    normals(0, 1, nodeNo) = bn.normal[1];
    
  }
  writeboundary.add_scalar_variables({"gamma"}, {gamma});
  writeboundary.add_vector_variables({"normals"}, {normals});
  writeboundary.write();
  // MPI_Finalize();
  // return 0;

  //=====================================================================================
  //
  //                                 DEFINE RHEOLOGY
  //
  //=====================================================================================
  Rans<LT> rans(input, myRank);

  //=====================================================================================
  //
  //                               MACROSCOPIC FIELDS
  //
  //=====================================================================================

  //                                   Viscosity
  //------------------------------------------------------------------------------------- Viscosity
  ScalarField viscosity(1, grid.size());

  ScalarField gammaDot(1, grid.size());
  ScalarField srcK(1, grid.size());
  ScalarField srcE(1, grid.size());
  //                                    Density
  // //------------------------------------------------------------------------------------- Density
  ScalarField rho(1, grid.size());

  //                           Initiate density from file
  //------------------------------------------------------------------------------------- Initiate density from file
  // vtklb.toAttribute("init_rho");
  for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    rho(0, n) = 1.0; //vtklb.getScalarAttribute<lbBase_t>();
  }

  //                                   Velocity
  // //------------------------------------------------------------------------------------- Velocity
  VectorField<LT> vel(1, grid.size());

  VectorField<LT> ForceField(1, grid.size());

  //                               Initiate velocity
  //------------------------------------------------------------------------------------- Initiate velocity
  for (auto nodeNo : bulkNodes)
  {
    for (int d = 0; d < LT::nD; ++d)
      vel(0, d, nodeNo) = 0.0;
  }

  //                                   rhoK field
  //------------------------------------------------------------------------------------- rhoK field
  ScalarField rhoK(1, grid.size());

  //                         Initiate rhoK field from file
  //------------------------------------------------------------------------------------- Initiate rhoK field from file
  // vtklb.toAttribute("init_rhoK");
  for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    rhoK(0, n) = kInlet; // vtklb.getScalarAttribute<lbBase_t>();
  }

  //                                 rhoEpsilon field
  //------------------------------------------------------------------------------------- rhoEpsilon field
  ScalarField rhoEpsilon(1, grid.size());

  //                        Initiate rhoEpsilon field from file
  //------------------------------------------------------------------------------------- Initiate rhoEpsilon field from file
  // vtklb.toAttribute("init_rhoEpsilon");
  for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    rhoEpsilon(0, n) = epsilonInlet; //  vtklb.getScalarAttribute<lbBase_t>();
  }

  //=====================================================================================
  //
  //                                      SETUP BOUNDARY
  //
  //=====================================================================================
  HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

  //=====================================================================================
  //
  //                                    INITIATE LB FIELDS
  //
  //=====================================================================================

  //                                        flow LB field
  //------------------------------------------------------------------------------------- flow LB field
  LbField<LT> f(1, grid.size());
  LbField<LT> fTmp(1, grid.size());

  //                                        rho k LB field
  //------------------------------------------------------------------------------------- rho k LB field
  LbField<LT> g(1, grid.size());
  LbField<LT> gTmp(1, grid.size());

  //                                     rho epsilon LB field
  //------------------------------------------------------------------------------------- rho epsilon LB field
  LbField<LT> h(1, grid.size());
  LbField<LT> hTmp(1, grid.size());

  //                                initiate LB distributions
  //------------------------------------------------------------------------------------- initiate LB distributions
  for (auto nodeNo : bulkNodes)
  {
    for (int q = 0; q < LT::nQ; ++q)
    {
      f(0, q, nodeNo) = LT::w[q] * rho(0, nodeNo);
      g(0, q, nodeNo) = LT::w[q] * rhoK(0, nodeNo);
      h(0, q, nodeNo) = LT::w[q] * rhoEpsilon(0, nodeNo);
    }
  }

  //=====================================================================================
  //
  //                                      OUTPUT VTK
  //
  //=====================================================================================

  Output<LT> output(grid, bulkNodes, outputDir2, myRank, nProcs);
  output.add_file("test_lb_run2");
  output.add_scalar_variables({"viscosity", "rho", "rhoK", "rhoEpsilon", "gammaDotTilde", "srcK", "srcEpsilon"}, {viscosity, rho, rhoK, rhoEpsilon, gammaDot, srcK, srcE});
  output.add_vector_variables({"vel"}, {vel});

  Output<LT, int> geo(grid.pos(), outputDir2, myRank, nProcs, "geo", nodes.geo(grid, vtklb));
  geo.write();

  //=====================================================================================
  //
  //                                      MAIN LOOP
  //
  //=====================================================================================

  //                              Start time step iterations
  //------------------------------------------------------------------------------------- Start time step iterations
  // begin clock
  //------------------------------------------------------------------------------------- begin clock
  auto startRun = std::chrono::high_resolution_clock::now();
  for (int i = 0; i <= nIterations; i++)
  {

    const lbBase_t ramp{0.5 * (1 - std::cos(3.14159 * std::min(i, rampTimesteps) / rampTimesteps))};

    //                        Calculate macroscopic values for all nodes
    //------------------------------------------------------------------------------------- Calculate macroscopic values for all nodes

    /*
    for (auto nodeNo: bulkNodes) {
    }

    //                               Communicate rho fields
    //------------------------------------------------------------------------------------- Communicate rho fields
    //(FILL IN COMMUNICATION HERE)
    */

    //------------------------------------------------------------------------------------- Regularized mean over neighbors
    for (auto nodeNo : bulkNodes) {
      int xpos = grid.pos(nodeNo, 0);
      int xposMax = vtklb.getGlobaDimensions(0) - 3; 
      if ( (xpos > 0) && (xpos < xposMax) && nodes.isFluidBoundary(nodeNo)  && nodes.isMyRank(nodeNo)) {
        const auto fNode = f(0, nodeNo);
        const auto gNode = g(0, nodeNo);
        const auto hNode = h(0, nodeNo);

        std::vector<lbBase_t> fMean(LT::nQ, 0.0);
        std::vector<lbBase_t> gMean(LT::nQ, 0.0);
        std::vector<lbBase_t> hMean(LT::nQ, 0.0);

        lbBase_t norm = 0.0;
        for (auto n: grid.neighbor(nodeNo)) { 
          const int x = grid.pos(n, 0);
          if (nodes.isFluid(n) && nodes.isMyRank(n) && (x > 0) && (x < xposMax)) {
            norm += 1.0;
            for (int q=0; q < LT::nQ; ++q) {
              fMean[q] += f(0, q, n);
              gMean[q] += g(0, q, n);
              hMean[q] += h(0, q, n);
            }
          }
        }
        for (int q=0; q < LT::nQ; ++q) {
          fTmp(0, q, nodeNo) = fMean[q]/norm;
          gTmp(0, q, nodeNo) = gMean[q]/norm;
          hTmp(0, q, nodeNo) = hMean[q]/norm;
        }
      }
      else {
        for (int q=0; q < LT::nQ; ++q) {
          fTmp(0, q, nodeNo) = f(0, q, nodeNo);
          gTmp(0, q, nodeNo) = g(0, q, nodeNo);
          hTmp(0, q, nodeNo) = h(0, q, nodeNo);
        }

      }

    }

    f.swapData(fTmp); // flow LBfield
    g.swapData(gTmp); // rhoK LBfield
    h.swapData(hTmp); // rhoEpsilon LBfield

    for (auto nodeNo : bulkNodes)
    {
      //                           Copy of local LB distribution
      //------------------------------------------------------------------------------------- Copy of local LB distribution
      const auto fNode = fRegularized<LT>(f(0, nodeNo), 0); //f(0, nodeNo);
      const auto gNode = fRegularized<LT>(g(0, nodeNo), 0); //f(0, nodeNo);
      const auto hNode = fRegularized<LT>(h(0, nodeNo), 0); //f(0, nodeNo);
      /* const auto fNode = f(0, nodeNo);
      const auto gNode = g(0, nodeNo);
      const auto hNode = h(0, nodeNo);*/

      //                                    Macroscopic values
      //------------------------------------------------------------------------------------- Macroscopic values
      const lbBase_t rhoNode = calcRho<LT>(fNode);
      const auto velNode = calcVel<LT>(fNode, rhoNode, bodyForce(0, 0));
      const lbBase_t u2 = LT::dot(velNode, velNode);
      const auto cu = LT::cDotAll(velNode);

      rans.apply(fNode, rhoNode, velNode, u2, cu, bodyForce(0, 0), 0.0, gNode, hNode, calcRho<LT>(gNode), calcRho<LT>(hNode));

      const lbBase_t rhoKNode = rans.rhoK();
      const lbBase_t rhoENode = rans.rhoE();

      const lbBase_t tauNode = rans.tau();
      //const lbBase_t tauNode = 0.5 + 0.001; //rans.tau();
      const lbBase_t tauKNode = rans.tauK();
      const lbBase_t tauENode = rans.tauE();

      //                         Save density and velocity for printing
      //------------------------------------------------------------------------------------- Save density and velocity for printing
      rho(0, nodeNo) = rhoNode;
      vel.set(0, nodeNo) = velNode;
      rhoK(0, nodeNo) = rhoKNode;
      rhoEpsilon(0, nodeNo) = rhoENode;
      viscosity(0, nodeNo) = rho(0, nodeNo) * LT::c2 * (tauNode - 0.5);

      gammaDot(0, nodeNo) = rans.gammaDot();
      srcK(0, nodeNo) = rans.sourceK();
      srcE(0, nodeNo) = rans.sourceE();

      //                                    Macroscopic values
      //------------------------------------------------------------------------------------- Macroscopic values
      // const lbBase_t rhoNode = calcRho<LT>(fNode);
      // const auto velNode = calcVel<LT>(fNode, rhoNode, bodyForce(0, 0));
      // const lbBase_t rhoKNode = calcRho<LT>(gNode);
      // const lbBase_t rhoEpsilonNode = calcRho<LT>(hNode);

      //                         Save density and velocity for printing
      //------------------------------------------------------------------------------------- Save density and velocity for printing
      // rho(0, nodeNo) = rhoNode;
      // vel.set(0, nodeNo) = velNode;
      // rhoK(0, nodeNo) = rhoKNode;
      // rhoEpsilon(0, nodeNo) = rhoEpsilonNode;

      //                                    BGK-collision term
      //------------------------------------------------------------------------------------- BGK-collision term
      // const lbBase_t u2 = LT::dot(velNode, velNode);
      // const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);

      // auto omegaBGK = rans.omegaBGK(fNode, rhoNode, velNode, u2, cu, bodyForce(0, 0), 0);

      //                                    LB interaction Terms
      //------------------------------------------------------------------------------------- LB interaction Terms
      ForceField.set(0, nodeNo) = bodyForce(0, 0);

      const lbBase_t uF = LT::dot(velNode, bodyForce(0, 0));
      const auto cF = LT::cDotAll(bodyForce(0, 0));

      /*
      //                                           SRT
      //------------------------------------------------------------------------------------- TRT

      //                                      Omegas: Flow
      //------------------------------------------------------------------------------------- Omegas: Flow

      const auto deltaOmegaF = calcDeltaOmegaF<LT>(tauNode, cu, uF, cF);

      const auto feqNode = calcfeq<LT>(rhoNode, u2, cu);
      const auto omegaBGK = calcOmegaBGK_TEST<LT>(fNode, feqNode, tauNode);

      //                                      Omegas: K & E
      //------------------------------------------------------------------------------------- Omegas: K & E

      const auto dOmegaFK = calcDeltaOmegaFDiff<LT>(tauKNode, rhoKNode/rhoNode, cu, uF, cF);
      const auto dOmegaFE = calcDeltaOmegaFDiff<LT>(tauENode, rhoENode/rhoNode, cu, uF, cF);

      const auto dOmegaSourceK = calcDeltaOmegaR<LT>(tauKNode, cu, rans.sourceK());
      const auto dOmegaSourceE = calcDeltaOmegaR<LT>(tauKNode, cu, rans.sourceE());

      const auto geqNode = calcfeq<LT>(rhoKNode, u2, cu);
      const auto omegaBGK_K = calcOmegaBGK_TEST<LT>(gNode, geqNode, tauKNode);

      const auto heqNode = calcfeq<LT>(rhoENode, u2, cu);
      const auto omegaBGK_E = calcOmegaBGK_TEST<LT>(hNode, heqNode, tauENode);
      */

      //                                           TRT
      //------------------------------------------------------------------------------------- TRT

      //                                      Omegas: Flow
      //------------------------------------------------------------------------------------- Omegas: Flow

      const auto deltaOmegaF = calcDeltaOmegaFTRT<LT>(tauNode, 1.0, cu, uF, cF);
      const auto omegaBGK = calcOmegaBGKTRT<LT>(fNode, tauNode, 1.0, rhoNode, u2, cu);

      //                                      Omegas: K & E
      //------------------------------------------------------------------------------------- Omegas: K & E

/*      const auto dOmegaFK = calcDeltaOmegaFDiffTRT<LT>(tauKNode, 1.0, rhoKNode / rhoNode, cu, uF, cF);
      const auto dOmegaFE = calcDeltaOmegaFDiffTRT<LT>(tauENode, 1.0, rhoENode / rhoNode, cu, uF, cF);

      const auto dOmegaSourceK = calcDeltaOmegaRTRT<LT>(tauKNode, 1.0, cu, rans.sourceK());
      const auto dOmegaSourceE = calcDeltaOmegaRTRT<LT>(tauKNode, 1.0, cu, rans.sourceE());

      const auto omegaBGK_K = calcOmegaBGKTRT<LT>(gNode, tauKNode, 1.0, rhoKNode, u2, cu);
      const auto omegaBGK_E = calcOmegaBGKTRT<LT>(hNode, tauENode, 1.0, rhoENode, u2, cu);
*/
      const auto dOmegaFK = calcDeltaOmegaFDiffTRT<LT>(1.0, tauKNode, rhoKNode / rhoNode, cu, uF, cF);
      const auto dOmegaFE = calcDeltaOmegaFDiffTRT<LT>(1.0, tauENode, rhoENode / rhoNode, cu, uF, cF);

      const auto dOmegaSourceK = calcDeltaOmegaRTRT<LT>(1.0, tauKNode, cu, rans.sourceK());
      const auto dOmegaSourceE = calcDeltaOmegaRTRT<LT>(1.0, tauKNode, cu, rans.sourceE());

      const auto omegaBGK_K = calcOmegaBGKTRT<LT>(gNode, 1.0, tauKNode, rhoKNode, u2, cu);
      const auto omegaBGK_E = calcOmegaBGKTRT<LT>(hNode, 1.0, tauENode, rhoENode, u2, cu);


      //                               Collision and propagation
      //------------------------------------------------------------------------------------- Collision and propagation
      fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);
      gTmp.propagateTo(0, nodeNo, gNode + omegaBGK_K + dOmegaFK + dOmegaSourceK, grid);
      hTmp.propagateTo(0, nodeNo, hNode + omegaBGK_E + dOmegaFE + dOmegaSourceE, grid);

    } //------------------------------------------------------------------------------------- End bulkNodes




    rans.solidBnd(f, fTmp, rho, vel, viscosity, g, gTmp, rhoK, h, hTmp, rhoEpsilon, bndInterp, nodes,grid);

    //                                   Swap data_ from fTmp to f etc.
    //------------------------------------------------------------------------------------- Swap data_ from fTmp to f etc.
    f.swapData(fTmp); // flow LBfield
    g.swapData(gTmp); // rhoK LBfield
    h.swapData(hTmp); // rhoEpsilon LBfield

    //=====================================================================================
    //
    //                                      BOUNDARY CONDITIONS
    //
    //=====================================================================================

    //                                            MPI
    //------------------------------------------------------------------------------------- MPI
    mpiBoundary.communicateLbField(0, f, grid);
    mpiBoundary.communicateLbField(0, g, grid);
    mpiBoundary.communicateLbField(0, h, grid);

    //                                   Half way bounce back
    //------------------------------------------------------------------------------------- Half way bounce back
    // bounceBackBnd.apply(f, grid);
    // bounceBackBnd.apply(g, grid);
    // bounceBackBnd.apply(h, grid);

    rans.zouHeFixedVelocityLeftBnd(inletBoundaryNodes, f, u_ref * ramp);
    // rans.zouHeFixedValueLeftBnd(inletBoundaryNodes, g, rho, kInlet * ramp + (1 - ramp) * 1e-4);
    // rans.zouHeFixedValueLeftBnd(inletBoundaryNodes, h, rho, epsilonInlet * ramp + (1 - ramp) * 1.6e-8);
    rans.zouHeFixedValueLeftBnd(inletBoundaryNodes, g, rho, kInlet );
    rans.zouHeFixedValueLeftBnd(inletBoundaryNodes, h, rho, epsilonInlet );


    // rans.zouHeFixedValueRightBnd(outletBoundaryNodes, f, 1.0, ForceField, grid);
    // rans.zouHeOpenRightBnd(outletBoundaryNodes, g, rhoK, ForceField, grid);
    // rans.zouHeOpenRightBnd(outletBoundaryNodes, h, rhoEpsilon, ForceField, grid);

    rans.copyDistBnd(outletBoundaryNodes, f, 4, grid);
    rans.copyDistBnd(outletBoundaryNodes, g, 4, grid);
    rans.copyDistBnd(outletBoundaryNodes, h, 4, grid);


    //=====================================================================================
    //
    //                                      WRITE TO FILE
    //
    //=====================================================================================
    if (int(0.2 * nItrWrite) > 0)
    {
      if (((i % int(0.2 * nItrWrite)) == 0) && myRank == 0)
      {

        std::cout << "Iteration: " << i << std::endl;
      }
    }
    else if (myRank == 0)
    {
      std::cout << "Iteration: " << i << std::endl;
    }

    if (((i % nItrWrite) == 0))
    {

      output.write(i);
      if (myRank == 0)
      {
        std::cout << "PLOT AT ITERATION : " << i << std::endl;
      }
    }

  } //------------------------------------------------------------------------------------- End time step iterations

  // end clock
  //------------------------------------------------------------------------------------- end clock
  auto endRun = std::chrono::high_resolution_clock::now();
  if (myRank == 0)
  {
    double runTime = std::chrono::duration_cast<std::chrono::nanoseconds>(endRun - startRun).count();
    runTime *= 1e-9; // from nanoseconds to seconds
    std::cout << "Run time = " << std::fixed << runTime << std::setprecision(3) << " sec" << std::endl;
  }

  MPI_Finalize();

  return 0;
}
