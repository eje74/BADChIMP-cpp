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
#define LT D3Q19
#define VTK_CELL VTK::pixel

//=====================================================================================
//
//                                 FIND BOUNDARY LINKS
//
//=====================================================================================
//                                SOLID-FLUID BOUNDARY
//------------------------------------------------------------------------------------- FIND SOLID-FLUID BOUNDARY
template<typename DXQY>
std::vector<std::vector<int>> findSolidFluidLinks(const Nodes<DXQY> &nodes, const Grid<DXQY> &grid) 
{
    std::vector<std::vector<int>> ret; // List of node numbers to all solid boundary nodes for myRank process
    for (int n = 1; n < nodes.size(); n++) 
    { // Loop over all grid nodes excpet the default node (node number = 0)
      int isBoundary = (nodes.getTag(n) >> 3) & 1;
      if (isBoundary && nodes.isMyRank(n))
      {
        for (int q = 0; q < DXQY::nQNonZero_; ++q) 
        {
          auto nn = grid.neighbor(q, n); // node neighbor
          int nnTag = nodes.getTag(nn) & 3; // Returns either 0, 1, 2, 3 based on the two smallest bits
          if (nnTag == 0)  // Add to list 
          {
            auto q_rev = DXQY::reverseDirection(q);
            ret.emplace_back(std::initializer_list<int>{n, q_rev, nn, q});
          }
        }
      }
    }

    return ret;
}

//                                PRESSURE-FLUID BOUNDARY
//------------------------------------------------------------------------------------- FIND PRESSURE-FLUID BOUNDARY
template<typename DXQY>
std::vector<std::vector<int>> findPressureFluidLinks(const Nodes<DXQY> &nodes, const Grid<DXQY> &grid) 
{
    std::vector<std::vector<int>> ret; // List of node numbers to all pressure boundary nodes for myRank process
    for (int n = 1; n < nodes.size(); n++) 
    { // Loop over all grid nodes excpet the default node (node number = 0)
      int isBoundary = (nodes.getTag(n) >> 4) & 1;
      if (isBoundary && nodes.isMyRank(n))
      {
        for (int q = 0; q < DXQY::nQNonZero_; ++q) 
        {
          auto nn = grid.neighbor(q, n); // node neighbor
          int nnTag = nodes.getTag(nn) & 3; // Returns either 0, 1, 2, 3 based on the two smallest bits
          if (nnTag == 3)  // Add to list if pressure links
          {
            auto q_rev = DXQY::reverseDirection(q);
            ret.emplace_back(std::initializer_list<int>{n, q_rev, nn, q});
          }
        }
      }
    }

    return ret;
}

//=====================================================================================
//
//                                    BOUNDARY CONDITIONS
//
//=====================================================================================
//                                SOLID-FLUID BOUNDARY
//------------------------------------------------------------------------------------- SOLID-FLUID BOUNDARY
template<typename DXQY>
void applySolidFluidBoundary(LbField<DXQY> &f, const std::vector<std::vector<int>> &boundaryLinks)
{
  for (auto link: boundaryLinks)
  {
    int nodeFluid = link[0];
    int qUnknown = link[1];
    int nodeWall = link[2];
    int qKnown = link[3];

    f(0, qUnknown, nodeFluid) = f(0, qKnown, nodeWall);
  }
}

//                                SOLID-FLUID BOUNDARY
//------------------------------------------------------------------------------------- SOLID-FLUID BOUNDARY
template<typename DXQY>
void applyPressureFluidBoundary(LbField<DXQY> &f, const VectorField<DXQY> &vel, const lbBase_t rho_w, const std::vector<std::vector<int>> &boundaryLinks)
/* Here we are using anti bounce back.
 * eq. 5.53 (p. 200) in Krugers book, but without velocity interpolation
 */
{
  for (auto link: boundaryLinks)
  {
    int nodeFluid = link[0];
    int qUnknown = link[1];
    int nodeWall = link[2];
    int qKnown = link[3];

    lbBase_t u2 = DXQY::dot(vel(0, nodeFluid), vel(0, nodeFluid));
    lbBase_t cu = DXQY::cDotRef(qKnown, vel(0, nodeFluid));
    lbBase_t w = DXQY::w[qKnown];

    f(0, qUnknown, nodeFluid) = -f(0, qKnown, nodeWall) + 2*w*rho_w*(1 + 0.5*(DXQY::c4Inv*cu*cu - DXQY::c2Inv*u2));
  }
}


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
  // std::string mpiDir = chimpDir + "input/mpi/";
  std::string inputDir = chimpDir + "input/";
  std::string outputDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/Python/CSSR/RelPerm/VtkGeo/";
  std::string mpiDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/Python/CSSR/RelPerm/VtkGeo/";
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
  // Set node tags
  vtklb.toAttribute("nodetags");
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    int tag = vtklb.getScalarAttribute<int>();
    nodes.setTag(tag, n); 
  }
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
  //                               Relaxation time
  //------------------------------------------------------------------------------------- Relaxation time
  lbBase_t tau = input["fluid"]["tau"];
  //                                    Body force
  //------------------------------------------------------------------------------------- Body force
  VectorField<LT> bodyForce(1, 1, input["fluid"]["bodyforce"]);
  //=====================================================================================
  //
  //                               WRITE PARAMETERS TO SCREEN
  //
  //=====================================================================================
  if (myRank==0) 
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
  //=====================================================================================
  //
  //                                MACROSCOPIC FIELDS
  //
  //=====================================================================================

  //                                     Fluid Flow
  //===================================================================================== Fluid Flow
  //                                      Density
  //------------------------------------------------------------------------------------- Density
  ScalarField rho(1, grid.size());
    //                              Initiate density
  //------------------------------------------------------------------------------------- Initiate density
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) 
  {
    rho(0, n) = 1.0;
  }
  //                                      Sources
  //------------------------------------------------------------------------------------- Sources
  // ScalarField Rfield(nFluidFields, grid.size());
  // ScalarField Qfield(1, grid.size());
  //                                       Velocity
  //------------------------------------------------------------------------------------- Velocity
  VectorField<LT> vel(1, grid.size());
  //                                   Initiate velocity
  //------------------------------------------------------------------------------------- Initiate velocity
  for (auto nodeNo: bulkNodes) {    
    vel.set(0, nodeNo) = 0;
    // for (int dim=0; dim<LT::nD; ++dim){
    //   vel(0, dim, nodeNo) = 0;
    }
  // //                                       Force Field
  // //------------------------------------------------------------------------------------- Force Field
  // VectorField<LT> ForceField(1, grid.size());
  // //                                   Initiate Force Field
  // //------------------------------------------------------------------------------------- Initiate Force Field
  // for (auto fieldNo=0; fieldNo < ForceField.num_fields(); ++fieldNo) {
  //   for (auto nodeNo: bulkNodes) {
  //     ForceField.set(fieldNo, nodeNo) = 0;
  //   }
  // }
  //=====================================================================================
  //
  //                                SETUP BOUNDARY
  //
  //=====================================================================================
  //                                  Solid fluid
  //------------------------------------------------------------------------------------- Solid-fluid
  auto solidFluidLinks = findSolidFluidLinks(nodes, grid);
  //                                  pressure - fluid
  //------------------------------------------------------------------------------------- Pressure fluid
  auto pressureFluidLinks = findPressureFluidLinks(nodes, grid);
  //=====================================================================================
  //
  //                                LB FIELDS
  //
  //===================================================================================== Lb fields
  //                                 Declaration
  //------------------------------------------------------------------------------------- declaration 
  LbField<LT> f(1, grid.size()); 
  LbField<LT> fTmp(1, grid.size());
  //                           Initiate lb distributions
  //------------------------------------------------------------------------------------- Initiate lb distributions
  for (auto nodeNo: bulkNodes) 
  {    
    auto u2 = LT::dot(vel(0, nodeNo), vel(0, nodeNo));
    auto cu = LT::cDotAll(vel(0, nodeNo));
    f.set(0, nodeNo) = calcfeq<LT>(rho(0, nodeNo), u2, cu);
    fTmp.set(0, nodeNo) = 0;
  }
  //=====================================================================================
  //
  //                                  OUTPUT VTK
  //
  //=====================================================================================
  Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs); 
  output.add_file("lb_run");
  output.add_scalar_variables(
    {"rho"}, 
		{ rho});
  output.add_vector_variables(
    {"vel"}, 
		{ vel});
  //######################################################################################
  //
  //                                  MAIN LOOP
  //
  //######################################################################################
  for (int i = 0; i <= nIterations; i++) 
  {
    //   int rampTimesteps = 5000; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////NB!    
    //   const lbBase_t ramp{ 0.5 * (1-std::cos(PI*std::min(i, rampTimesteps)/rampTimesteps)) };
    //                                   Main calculation loop
    //------------------------------------------------------------------------------------- 
    for (auto nodeNo: bulkNodes) 
    {
  //     Qfield(0, nodeNo)=0.0;
  //     ForceField.set(0, nodeNo)=0.0;   
  //     ForceField.set(0, nodeNo) += bodyForce(0, 0);
      //                            Copy of local velocity distribution
      //------------------------------------------------------------------------------------- Copy of local velocity distribution
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
  //     //const auto feqNode = calcfeq<LT>(rhoNode, u2, cu);
  //     //const auto omegaBGKTot = newtonian.omegaBGK(tau, fNode, rhoNode, velNode, u2, cu, ForceField(0, nodeNo), Qfield(0, nodeNo));      
      const std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);            
      //                           Calculate the Guo-force correction
      //------------------------------------------------------------------------------------- Calculate the Guo-force correction
      const lbBase_t uF = LT::dot(velNode, bodyForce(0, 0));
      const std::valarray<lbBase_t> cF = LT::cDotAll(bodyForce(0, 0));
      const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);
  //     //                           Calculate the mass-source correction
  //     //------------------------------------------------------------------------------------- Calculate the mass-source correction
  //     const std::valarray<lbBase_t> deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tau, cu, u2, Qfield(0, nodeNo));
      //                               Collision and propagation
      //------------------------------------------------------------------------------------- Collision and propagation
      fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);
      // fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF + deltaOmegaQ0, grid);
    } //------------------------------------------------------------------------------------- End for bulkNodes
    //                           Swap lb fields
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
    //                             Solid-fluid boundary
    //------------------------------------------------------------------------------------- solid-fluid
    applySolidFluidBoundary(f, solidFluidLinks);    
    //                             Pressure-fluid boundary
    //------------------------------------------------------------------------------------- pressure-fluid
    applyPressureFluidBoundary(f, vel, 1.0, pressureFluidLinks);
    //=====================================================================================
    //
    //                                 WRITE TO FILE
    //
    //=====================================================================================
    if ( ((i % nItrWrite) == 0)  ) {
      output.write(i);
      if (myRank==0) {
	      std::cout << "PLOT AT ITERATION: " << i << std::endl;
      }
   }  
  } //----------------------------------------------------------------------------------------  End for nIterations
  
  MPI_Finalize();
  
  return 0;
}
