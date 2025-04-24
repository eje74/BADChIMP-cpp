
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
//------------------------------------------------------------------------------------- Link structure
struct Boundarylink
{
  int fluidNode, fluidDir;
  int solidNode, solidDir;
};

//------------------------------------------------------------------------------------- Count boundary links
template<typename DXQY>
int boundaryLinksCount(
  int tag, 
  const std::vector<int> &boundaryTags,
  const std::vector<int> &bulkNodes, 
  const Nodes<DXQY> &nodes, 
  const Grid<DXQY> &grid)
  /* Returns the number of boundary links with the given tag */
{
  int cnt = 0;
  for (auto nodeNo: bulkNodes) 
  {
    for (int q=0; q < DXQY::nQ; ++q) 
    {
      int nodeNeig = grid.neighbor(q, nodeNo);
      if (nodes.isSolid(nodeNeig) && boundaryTags[nodeNeig] == tag)
        cnt += 1; 
    }
  }
  return cnt;
}

//------------------------------------------------------------------------------------- create boundary link vector
template<typename DXQY>
std::vector<Boundarylink> makeBoundaryLinks(
  int tag, 
  const std::vector<int> &boundaryTags,
  const std::vector<int> &bulkNodes, 
  const Nodes<DXQY> &nodes, 
  const Grid<DXQY> &grid)
  /* Returns a vector of BoundaryLinks containing all links from fluid 
     to solids nodes with the given tag */
{
  int linksNum = boundaryLinksCount(tag, boundaryTags, bulkNodes, nodes, grid);
  std::vector<Boundarylink> ret(linksNum);
  int cnt = 0;
  for (auto nodeNo: bulkNodes) 
  {
    for (int q=0; q < DXQY::nQ; ++q) 
    {
      int nodeNeig = grid.neighbor(q, nodeNo);
      if (nodes.isSolid(nodeNeig) && boundaryTags[nodeNeig] == tag)
      {
        ret[cnt].solidNode = nodeNeig;
        ret[cnt].solidDir = q;
        ret[cnt].fluidNode = nodeNo;
        ret[cnt].fluidDir = DXQY::reverseDirection(q);
        cnt += 1;
      } 
    }
  }
  return ret;
}
//------------------------------------------------------------------------------------- solid-fluid
template<typename DXQY>
void bounceBackApply(
  LbField<DXQY> &f, 
  const std::vector<Boundarylink> &links)
{
  for (auto link: links)
    f(0, link.fluidDir, link.fluidNode) = f(0, link.solidDir, link.solidNode);
}
//------------------------------------------------------------------------------------- pressure-fluid
template<typename DXQY>
void antiBounceBackApply(
  lbBase_t rho,
  LbField<DXQY> &f, 
  const VectorField<DXQY> &vel,
  const std::vector<Boundarylink> &links)
{
  for (auto link: links) 
  {
    const std::valarray<lbBase_t> velNode = vel(0, link.fluidNode);
    int q = link.solidDir;
    lbBase_t cu = DXQY::cDotRef(q, velNode);
    lbBase_t u2 = DXQY::dot(velNode, velNode);
    lbBase_t val = 2.0*DXQY::w[q]*rho*(1.0 + 0.5*(DXQY::c4Inv*cu*cu - DXQY::c2Inv*u2));

    f(0, link.fluidDir, link.fluidNode) = -f(0, link.solidDir, link.solidNode) + val;
  }
}

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

  //===================================================================================== Setup output file
  std::string outputFile = outputDir + "sumvel.dat";

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
  //------------------------------------------------------------------------------------- Filename base
  std::string filenamebase = input["filenames"]["basename"]; 
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
    std::cout << "  basename = " << filenamebase << std::endl;
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
  std::vector<int> geoTags(grid.size());
  vtklb.toAttribute("geo_tag");
  for (int nodeNo=vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo)
  {
    geoTags[nodeNo] = vtklb.getScalarAttribute<int>();
  }
  //------------------------------------------------------------------------------------- test-tag
  //------------------------------------------------------------------------------------- solid-fluid
  auto wallBndLinks = makeBoundaryLinks(0, geoTags, bulkNodes, nodes, grid);
  //------------------------------------------------------------------------------------- pressure-fluid
  auto pressureBndLinks = makeBoundaryLinks(1, geoTags, bulkNodes, nodes, grid);


  ScalarField tagNeig(2, grid.size());
  for (auto nodeNo: bulkNodes)
  {
    int cnt1 = 0;
    int cnt2 = 0;
    for (auto neigNo: grid.neighbor(nodeNo)) {
      if ( geoTags[neigNo] == 1 ) cnt1 += 1;
      if ( geoTags[neigNo] == 2 ) cnt2 += 1; 
    }
    tagNeig(0, nodeNo) = cnt1;
    tagNeig(1, nodeNo) = cnt2; 

  }

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
  // ==================================================================================== Output
  // ------------------------------------------------------------------------------------ vtk
  Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs);
  output.add_file("lb_run");
  output.add_scalar_variables(
      {"rho", "tag_neig"},
      {rho, tagNeig});
  output.add_vector_variables(
      {"vel"},
      {vel});
  // ------------------------------------------------------------------------------------ Sum velocites
  MPI_Status status;
  MPI_File fh;
  if (myRank == 0)
    MPI_File_open(MPI_COMM_SELF, outputFile.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  // ###################################################################################### MAIN LOOP
  std::time_t start, end;

  if (myRank == 0)
  {
    std::time(&start);
  }

  for (int i = 0; i <= nIterations; i++)
  {
    //------------------------------------------------------------------------------------- Rampup
    int rampTimesteps = 1000; 
    const lbBase_t ramp{ 0.5 * (1-std::cos(3.14159*std::min(i, rampTimesteps)/rampTimesteps)) };
    std::valarray<lbBase_t> forceNode = bodyForce(0, 0);
    forceNode[2] = ramp*forceNode[2];

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
    // bounceBackBnd.apply(f, grid);
    bounceBackApply(f, wallBndLinks);
    //------------------------------------------------------------------------------------- pressure-fluid
    antiBounceBackApply(1.0, f, vel, pressureBndLinks);
    //------------------------------------------------------------------------------------- fluid-fluid

    //===================================================================================== Write
    if (((i % nItrWrite) == 0))
    {
    //------------------------------------------------------------------------------------- write to screen
    if (myRank == 0)
      {
        std::cout << "ITERATION:  " << i << std::endl;
      }
      /// output.write(i);
      //----------------------------------------------------------------------------------- sum velocity
      std::vector<lbBase_t> sumVelLocal(LT::nD, 0.0);
      for (auto nodeNo: bulkNodes) {
        for (int nd=0; nd<LT::nD; ++nd)
          sumVelLocal[nd] += vel(0, nd, nodeNo);
      }
      std::vector<lbBase_t> sumVelGlobal(LT::nD, 0.0);
      MPI_Allreduce(sumVelLocal.data(), sumVelGlobal.data(), LT::nD, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //------------------------------------------------------------------------------------- write sum velocity
      const int precision = 8;
      if (myRank == 0) {
        std::ostringstream oss;
        oss << std::to_string(i) << ", " << std::setprecision(precision) << sumVelGlobal[0];
        for (int d = 1; d < 3; ++d)
          oss << ", " << std::setprecision(precision) << sumVelGlobal[d];
        oss << "\n";
        std::string ws = oss.str();
        MPI_File_write(fh, ws.c_str(), ws.size(), MPI_CHAR, &status);
      }

    }
  } //----------------------------------------------------------------------------------------  End for nIterations
  if (myRank == 0)
  {
    std::time(&end);
    double runTime = double(end - start);
    std::cout << "Run time = " << runTime << "\n";
  }

  //------------------------------------------------------------------------------------- write fields to file
  output.write(nIterations);

  if (myRank == 0)  
    MPI_File_close(&fh);
  MPI_Finalize();

  return 0;
}
