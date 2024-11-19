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

template<typename DXQY>
std::vector<int> findPressureFluidNodes(const Nodes<DXQY> &nodes, const Grid<DXQY> &grid) 
{
    std::vector<int> ret; // List of node numbers to all pressure boundary nodes for myRank process
    for (int n = 1; n < nodes.size(); n++) 
    { // Loop over all grid nodes excpet the default node (node number = 0)
      int isBoundary = (nodes.getTag(n) >> 4) & 1;
      if (isBoundary && nodes.isMyRank(n))
      {
        ret.push_back(n);
      }
    }

    return ret;
}

//                                FLUID-FLUID BOUNDARY
//------------------------------------------------------------------------------------- FIND PRESSURE-FLUID BOUNDARY
template<typename DXQY>
std::vector<std::vector<int>> findFluidFluidLinks(const Nodes<DXQY> &nodes, const Grid<DXQY>&grid)
{
  std::vector<std::vector<int>> ret;
  // Loop over all grid nodes except the default node (node number = 0)  
  for (int n = 1; n < nodes.size(); n++) 
  { 
    int isFluidBoundary = (nodes.getTag(n) >> 2) & 1;
    if (isFluidBoundary && nodes.isMyRank(n))
    {
      int phase = nodes.getTag(n) & 3;
      // Here use the fluid phase 1 as base, so that we do not double count the
      // fluid-fluid boundary links.
      if (phase == 1)
      {
        for (int q = 0; q < DXQY::nQNonZero_; ++q) 
        {
          auto nn = grid.neighbor(q, n); // node neighbor
          int nnPhase = nodes.getTag(nn) & 3; // Returns either 0, 1, 2, 3 based on the two smallest bits
          if (nnPhase == 2)  // Add to list if pressure links
          {
            auto q_rev = DXQY::reverseDirection(q);
            ret.emplace_back(std::initializer_list<int>{n, q_rev, nn, q});
          }
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
//                                FLUID-FLUID BOUNDARY
//------------------------------------------------------------------------------------- PRESSURE-FLUID BOUNDARY
template<typename DXQY>
void applyFluidFluidBoundary(LbField<DXQY> &f, const VectorField<DXQY> &vel, const VectorField<DXQY> &norm, const std::vector<std::vector<int>> &boundaryLinks, Nodes<DXQY> &nodes)
{
  for (auto link: boundaryLinks)
  {
    int nodeFluid1 = link[0];
    int qUnknown1 = link[1];
    int nodeFluid2 = link[2];
    int qKnown1 = link[3];

    lbBase_t f1 = f(0, qUnknown1, nodeFluid1);
    std::valarray<lbBase_t> u_w = 0.5*(vel(0, nodeFluid1) + vel(0, nodeFluid2));
    std::valarray<lbBase_t> n_vec = 0.5*(norm(0, nodeFluid1) + norm(0, nodeFluid2));
    lbBase_t length = std::sqrt( n_vec[0]*n_vec[0] +  n_vec[1]*n_vec[1] + n_vec[2]*n_vec[2] );
    n_vec /= (length + 1e-14);
    lbBase_t cn = DXQY::cDotRef(qUnknown1, n_vec);
    lbBase_t cu = DXQY::cDotRef(qUnknown1, u_w);
    lbBase_t un = DXQY::dot(u_w, n_vec);
    lbBase_t w = DXQY::w[qUnknown1];
    int isbnd1 = 1 - ((nodes.getTag(nodeFluid1) >> 3) & 1);
    int isbnd2 = 1 - ((nodes.getTag(nodeFluid2) >> 3) & 1);
    lbBase_t dfu = 6*w*(cu - un*cn)*isbnd1*isbnd2;
    f(0, qUnknown1, nodeFluid1) = f(0, qKnown1, nodeFluid2) + dfu;
    f(0, qKnown1, nodeFluid2) = f1 - dfu;
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

  //                              Node tag attributes
  //------------------------------------------------------------------------------------- nodetags attribute
  //       0: solid
  //       1: fluid phase 1
  //       2: fluid phase 2 
  //       3: pressure boundary' ghost node (treated as a solid node)
  //       4: fluid-fluid interface
  //       8: fluid-solid interface
  //      16: pressure boundary
  ScalarField tagsTmp(1, grid.size());
  vtklb.toAttribute("nodetags");
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    int tag = vtklb.getScalarAttribute<int>();
    nodes.setTag(tag, n); 
    tagsTmp(0, n) = tag;
  }
  //                              domain attributes
  //------------------------------------------------------------------------------------- domain labels
  ScalarField domains(1, grid.size());
  vtklb.toAttribute("domains");
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    int domainNumber = vtklb.getScalarAttribute<int>();
    domains(0, n) = domainNumber;
  }
  //                              force indicator
  //------------------------------------------------------------------------------------- force indicator
  ScalarField forceOn(1, grid.size());
  vtklb.toAttribute("force");
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    int val = vtklb.getScalarAttribute<int>();
    forceOn(0, n) = val;
  }
  //                              interior domains
  //------------------------------------------------------------------------------------- interior domains
  ScalarField interiorDomains(1, grid.size());
  std::vector<int> interiorDomainsLabel(grid.size(), 0);
  std::vector<lbBase_t> addMassSource(grid.size(), 0.0);
  vtklb.toAttribute("interior_domains");
  int localDomainLabelMax = 0;
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    int val = vtklb.getScalarAttribute<int>();
    interiorDomains(0, n) = val;
    interiorDomainsLabel[n] = val;
    if (nodes.isMyRank(n))
      if (val > localDomainLabelMax)
        localDomainLabelMax = val;
  }
  //                              normals
  //------------------------------------------------------------------------------------- normals
  VectorField<LT> normals(1, grid.size());
  vtklb.toAttribute("normal_x");
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    lbBase_t val = vtklb.getScalarAttribute<lbBase_t>();
    normals(0, 0, n) = val;
  }
  vtklb.toAttribute("normal_y");
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    lbBase_t val = vtklb.getScalarAttribute<lbBase_t>();
    normals(0, 1, n) = val;
  }
  vtklb.toAttribute("normal_z");
  for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    lbBase_t val = vtklb.getScalarAttribute<lbBase_t>();
    normals(0, 2, n) = val;
  }

  //                              setup communication for mass conservation
  //------------------------------------------------------------------------------------- init mass conservation
  int globalDomainLabelMax;
  MPI_Allreduce(&localDomainLabelMax, &globalDomainLabelMax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  std::cout << "Rank: " << myRank << " local max domain label = " << localDomainLabelMax << "  global domain max = " << globalDomainLabelMax << std::endl;
  // size of the interior of the "interior domains", that is, at
  // the nodes where mass is added or subtracted.
  std::vector<lbBase_t> massSourceScaleFactor(globalDomainLabelMax + 1, 0);
  std::vector<lbBase_t> massChangeLocal(globalDomainLabelMax + 1, 0);
  std::vector<lbBase_t> massChange(globalDomainLabelMax + 1, 0);
  {
    std::vector<lbBase_t> tmp(globalDomainLabelMax + 1, 0.0);
    for (int n = 1; n < grid.size(); ++n)
    {
      if (nodes.isMyRank(n)) 
      { 
        int label = interiorDomainsLabel[n];
        int tag = nodes.getTag(n);
        if ( (label > 0) &&  (tag < 3) ) 
        {
          tmp[label] += 1;
          addMassSource[n] = 1.0;
        } 
      }
    }
    // MPI Communicate number of mass sources per interior domain
    MPI_Allreduce(tmp.data(), massSourceScaleFactor.data(), globalDomainLabelMax+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    for (int i = 1; i < globalDomainLabelMax + 1; ++i) 
    {
      if ( massSourceScaleFactor[i] == 0.0) {
        std::cout << "ERROR: Interior domain " << i << " has no interior nodes!" << std::endl;
        exit(1);
      } 
      massSourceScaleFactor[i] = 1.0/massSourceScaleFactor[i];
    }
  }
  //                              setup calculation of mass flux
  //------------------------------------------------------------------------------------- calculation of mass flux
  auto pressureFluidNodes = findPressureFluidNodes(nodes, grid);
 
  Output<LT> outputTest(grid, bulkNodes, outputDir, myRank, nProcs);
  outputTest.add_file("geo");
  outputTest.add_scalar_variables(
    {"tags", "domains", "force", "interior_domains"},
    {tagsTmp, domains, forceOn, interiorDomains}
  );
  outputTest.add_vector_variables(
    {"normals"},
    {normals}
  );
  outputTest.write(0);


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
  //                                  fluid - fluid
  //------------------------------------------------------------------------------------- fluid fluid
  auto fluidFluidLinks = findFluidFluidLinks(nodes, grid);
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
  //                           Check convergence of rel.perm
  //------------------------------------------------------------------------------------- Check convergence of rel.perm
  std::vector<lbBase_t> oldMassFlux(2, 0.0);

  //######################################################################################
  //
  //                                  MAIN LOOP
  //
  //######################################################################################
  for (int i = 0; i <= nIterations; i++) 
  {
    //   int rampTimesteps = 5000; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////NB!    
    //   const lbBase_t ramp{ 0.5 * (1-std::cos(PI*std::min(i, rampTimesteps)/rampTimesteps)) };
    //                                   Mass conservation
    //------------------------------------------------------------------------------------- mass conservation
    // Sett local to zeros
    std::fill(massChangeLocal.begin(), massChangeLocal.end(), 0.0);
    for (auto nodeNo: bulkNodes) {
      const std::valarray<lbBase_t> fNode = f(0, nodeNo); 
      const lbBase_t rhoNode = calcRho<LT>(fNode);
      // Source term
      int label = interiorDomainsLabel[nodeNo];
      massChangeLocal[label] += 1.0 - rhoNode;
    } 
    MPI_Allreduce(massChangeLocal.data(), massChange.data(), massChange.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // Sett local to zeros
    std::fill(massChangeLocal.begin(), massChangeLocal.end(), 0.0);    
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
      lbBase_t rhoNode = calcRho<LT>(fNode);
      // Source term
      int label = interiorDomainsLabel[nodeNo];
      lbBase_t qMassConservation = 0.9*2*massSourceScaleFactor[label]*massChange[label]*addMassSource[nodeNo];
      // Add source term
      rhoNode += 0.5*qMassConservation;
      massChangeLocal[label] += 1.0 - rhoNode;
      const std::valarray<lbBase_t> forceNode = bodyForce(0, 0)*forceOn(0, nodeNo);
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
      //                           Calculate the mass-source correction
      //------------------------------------------------------------------------------------- Calculate the mass-source correction
      const std::valarray<lbBase_t> deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tau, cu, u2, qMassConservation);
      //                               Collision and propagation
      //------------------------------------------------------------------------------------- Collision and propagation
      // fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);
      fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF + deltaOmegaQ0, grid);
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
    mpiBoundary.communciateVectorField_TEST(vel);
    //                             Solid-fluid boundary
    //------------------------------------------------------------------------------------- solid-fluid
    applySolidFluidBoundary(f, solidFluidLinks);    
    //                             Pressure-fluid boundary
    //------------------------------------------------------------------------------------- pressure-fluid
    applyPressureFluidBoundary(f, vel, 1.0, pressureFluidLinks);
    //                             fluid-fluid boundary
    //------------------------------------------------------------------------------------- fluid-fluid
    applyFluidFluidBoundary(f, vel, normals, fluidFluidLinks, nodes);
    //=====================================================================================
    //
    //                                 WRITE TO FILE
    //
    //=====================================================================================
    if ( ((i % nItrWrite) == 0)  ) {
      output.write(i);
      std::vector<lbBase_t> massFluxLocal(2, 0.0);
      for (auto n: pressureFluidNodes)
      {
        int fluidPhase = (nodes.getTag(n) & 3) - 1;
        if ( (fluidPhase!= 0) && (fluidPhase != 1) ) {
          std::cout << "Fluid phase = " << fluidPhase << " in write mass flux" << std::endl;
          MPI_Finalize();
          exit(1);
        }
        massFluxLocal[fluidPhase] += vel(0, 2, n)*rho(0, n);
      }
      std::vector<lbBase_t>  massFluxGlobal(2, 0.0);
      MPI_Allreduce(massFluxLocal.data(), massFluxGlobal.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (myRank==0) {
        lbBase_t q1 = 0.5*massFluxGlobal[0];
        lbBase_t q2 = 0.5*massFluxGlobal[1];
        lbBase_t q1_change = (q1 - oldMassFlux[0])/(q1 + 1e-15);
        lbBase_t q2_change = (q2 - oldMassFlux[1])/(q2 + 1e-15);
	      std::cout << "PLOT AT ITERATION: " << i << std::endl;
        std::cout << "q1 = " << q1 << " (" << q1_change << ")" << std::endl;
        std::cout << "q2 = " << q2 << " (" << q2_change << ")" << std::endl;
        oldMassFlux[0] = q1;
        oldMassFlux[1] = q2;
      }
      //                           Update mass source
      //------------------------------------------------------------------------------------- Update mass source
      MPI_Allreduce(massChangeLocal.data(), massChange.data(), massChange.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      // Sett local to zeros
      std::fill(massChangeLocal.begin(), massChangeLocal.end(), 0.0);
   }  
  } //----------------------------------------------------------------------------------------  End for nIterations
  
  MPI_Finalize();
  
  return 0;
}
