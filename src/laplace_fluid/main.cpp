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
#include "./LBpressureboundary.h"

// SET THE LATTICE TYPE
//#define LT D2Q9
#define LT D3Q19

//#define VTK_CELL VTK::pixel

void zouHePressureBoundary(const std::vector<int> &bndNodes, LbField<LT> &f, const lbBase_t rho, const VectorField<LT> &force, const Grid<LT> &grid)
{
  for (const auto &nodeNo : bndNodes)
  {
    const std::valarray<lbBase_t> fNode = f(0, nodeNo);
    const std::valarray<lbBase_t> forceNode = force(0, nodeNo);
    const lbBase_t rho_ux = rho - 0.5 * forceNode[0] - (fNode[2] + fNode[6] + fNode[8] + 2 * (fNode[3] + fNode[4] + fNode[5]));
    f(0, 0, nodeNo) = fNode[4] + (2. / 3.) * rho_ux - (1. / 3.) * forceNode[0];
    f(0, 1, nodeNo) = fNode[5] + 0.5 * (fNode[6] - fNode[2]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] + (1. / 4.) * forceNode[1];
    f(0, 7, nodeNo) = fNode[3] + 0.5 * (fNode[2] - fNode[6]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] - (1. / 4.) * forceNode[1];

    // LT::qSum(f(0,nodeNo))
  }
}

void zouHePressureBoundaryRight(const std::vector<int> &bndNodes, LbField<LT> &f, const lbBase_t rho, const VectorField<LT> &force, const Grid<LT> &grid)
{
  for (const auto &nodeNo : bndNodes)
  {
    const std::valarray<lbBase_t> fNode = f(0, nodeNo);
    const std::valarray<lbBase_t> forceNode = force(0, nodeNo);
    const lbBase_t rho_ux = rho - 0.5 * forceNode[0] - (fNode[2] + fNode[6] + fNode[8] + 2 * (fNode[0] + fNode[1] + fNode[7]));
    f(0, 4, nodeNo) = fNode[0] + (2. / 3.) * rho_ux - (1. / 3.) * forceNode[0];
    f(0, 3, nodeNo) = fNode[7] + 0.5 * (fNode[6] - fNode[2]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] + (1. / 4.) * forceNode[1];
    f(0, 5, nodeNo) = fNode[1] + 0.5 * (fNode[2] - fNode[6]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] - (1. / 4.) * forceNode[1];

    // LT::qSum(f(0,nodeNo))
  }
}

int copyFile(std::string fn_src, std::string fn_dst)
{
  std::ifstream src(fn_src, std::ios::binary);
  std::ofstream dst(fn_dst, std::ios::binary);
  dst << src.rdbuf();

  return src && dst;
}

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
  std::string chimpDir = "../";
  // std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
  std::string mpiDir = chimpDir + "input/mpi/";
  std::string inputDir = chimpDir + "input/";
  std::string laplacePressureDir = mpiDir + "uturn/";

  Input input(inputDir + "input.dat");

  // const bool laplacePressureRun = static_cast<int>( input["iterations"]["laplaceforcing"]);

  /*    int dirNum = static_cast<int>( input["iterations"]["dirnum"]);
      std::string dirNum_s = dirNum < 10 ? "0" : "";
      dirNum_s += std::to_string(dirNum);
      std::string outputDir = chimpDir + "oscillatingoutput" + dirNum_s + "/";

      // ---------------------------------------------------------------------- Copy main and input
      copyFile(chimpDir + "src/laplace_fluid/main.cpp", outputDir + "main.cpp_back_up");
      copyFile(inputDir + "input.dat", outputDir + "intput.dat_back_up");
  */

  // const lbBase_t dP = input["fluid"]["dP"];
  // const lbBase_t freq = input["fluid"]["f"];

  const bool laplacePressureRun = true;
  std::string outputDir = chimpDir + "output/";

  // ***********************
  // SETUP GRID AND GEOMETRY
  // ***********************
  LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
  Grid<LT> grid(vtklb);
  Nodes<LT> nodes(vtklb, grid);
  BndMpi<LT> mpiBoundary(vtklb, nodes, grid);

  /* //obselete
  // Set bulk nodes
  vtklb.toAttribute("pressure_boundary");
  std::vector<int> mark(grid.size());
  for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
      const int val = vtklb.getScalarAttribute<int>();
      if (val == -1) {
          mark[nodeNo] = 1;
      }
      else {
          mark[nodeNo] = 0;
      }
  }
  */

  int maxBoundaryIndicator = addPressureBoundary("pressure_boundary", vtklb, nodes, grid);
  std::cout << "maxBoundaryIndicator = " << maxBoundaryIndicator << std::endl;

  // Set bulk nodes
  std::vector<int> bulkNodes = findBulkNodes(nodes);

  // *************
  // READ FROM INPUT
  // *************
  // Number of iterations
  int nIterations = static_cast<int>(input["iterations"]["max"]);
  // Write interval
  int nItrWrite = static_cast<int>(input["iterations"]["write"]);
  // Relaxation time
  lbBase_t tau = input["fluid"]["tau"];

  lbBase_t dP_abs = input["fluid"]["dP_abs"];
  std::valarray<lbBase_t> dP = dP_abs * inputAsValarray<lbBase_t>(input["fluid"]["dP_rel"]);
  lbBase_t PressurePeriod = input["fluid"]["PressurePeriod"];
  // lbBase_t dP = 0.0;

  //                                    Output directory number
  //------------------------------------------------------------------------------------- Output directory number
  std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));
  std::string outputDir2 = outputDir + "/outFl" + dirNum;

  // ******************
  // MACROSCOPIC FIELDS
  // ******************
  // Density
  ScalarField rho(1, grid.size());
  // Initiate density from file
  // vtklb.toAttribute("init_rho");
  for (int n = vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n)
  {
    rho(0, n) = 1.0; // vtklb.getScalarAttribute<lbBase_t>();
  }

  //--------------------------------------------------------------------------------- laplace force & laplace pressure

  VectorField<LT> laplaceForce(maxBoundaryIndicator, grid.size());
  ScalarField laplacePressure(maxBoundaryIndicator, grid.size());
  for (int n = 0; n < maxBoundaryIndicator; ++n)
  {
    VectorField<LT> jTmp(1, grid.size());
    ScalarField psiTmp(1, grid.size());
    std::string filename = "laplace_pressure_rank_" + std::to_string(myRank) + "_fieldnum_" + std::to_string(n);
    jTmp.readFromFile(mpiDir + filename);
    psiTmp.readFromFile(mpiDir + filename);
    for (auto &nodeNo : bulkNodes)
    {
      laplaceForce.set(n, nodeNo) = jTmp(0, nodeNo);
      laplacePressure(n, nodeNo) = psiTmp(0, nodeNo);
    }
  }

  ScalarField pressure(1, grid.size());
  ScalarField pPert(1, grid.size());

  // Velocity
  VectorField<LT> vel(1, grid.size());
  // Initiate velocity
  for (auto nodeNo : bulkNodes)
  {
    for (int d = 0; d < LT::nD; ++d)
      vel(0, d, nodeNo) = 0.0;
  }

  // Force
  VectorField<LT> force(1, grid.size());

  // ******************
  // SETUP BOUNDARY
  // ******************
  // HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

  // ---------------------------------------------------------------------------------- Read signed distance
  // need to be defined before pressure.
  // auto sd = readSignedDistance("signed_distance", vtklb, nodes, grid);

  fluidPressureBoundary<LT> pBnd(findFluidBndNodes(nodes), nodes, grid);

  std::vector<std::string> strCartDirs{"x", "y", "z"};
  std::valarray<lbBase_t> nullVector{0.0, 0.0, 0.0};

  VectorField<LT> pBndNorms(1, maxBoundaryIndicator);

  for (int i = 0; i < maxBoundaryIndicator; i++)
  {
    pBndNorms.set(0, i) = nullVector;
  }

  for (int i = 0; i < LT::nD; i++)
  {
    vtklb.toAttribute("boundary_normal_" + strCartDirs[i]);
    for (int n = 0; n < vtklb.numSubsetEntries(); ++n)
    {
      const auto ent = vtklb.getSubsetAttribure<lbBase_t>();
      int tag = nodes.getTag(ent.nodeNo);
      if (tag > 0)
      {
        pBndNorms(0, i, tag - 1) += ent.val;
      }
    }
  }

  for (int i = 0; i < maxBoundaryIndicator; i++)
  {
    lbBase_t tmpLength = 0;
    tmpLength = std::sqrt(LT::dot(pBndNorms(0, i), pBndNorms(0, i)));
    // std::cout << "Rank "<< myRank <<" Boundary "<< i << "( ";
    for (int dir = 0; dir < LT::nD; dir++)
    {
      pBndNorms(0, dir, i) /= tmpLength;
      // std::cout << pBndNorms(0, dir, i)<< " ";
    }
    // std::cout << ")" << std::endl;
  }
  pBndNorms(0, 0, 0) = 0.340088;
  pBndNorms(0, 1, 0) = -0.352946;
  pBndNorms(0, 2, 0) = -0.871647;

  for (int i = 0; i < maxBoundaryIndicator; i++)
  {
    std::cout << "Rank " << myRank << " Boundary " << i << "( ";
    for (int dir = 0; dir < LT::nD; dir++)
    {
      std::cout << pBndNorms(0, dir, i) << " ";
    }
    std::cout << ")" << std::endl;
  }

  // *********
  // LB FIELDS
  // *********
  LbField<LT> f(1, grid.size());
  LbField<LT> fTmp(1, grid.size());
  // initiate lb distributions
  for (auto nodeNo : bulkNodes)
  {
    for (int q = 0; q < LT::nQ; ++q)
    {
      f(0, q, nodeNo) = LT::w[q] * rho(0, nodeNo);
    }
  }

  ScalarField tagsField(1, grid.size());
  ScalarField nodeTypeField(1, grid.size());
  ScalarField nodeNoField(1, grid.size());

  for (auto nodeNo : bulkNodes)
  {
    tagsField(0, nodeNo) = 0;
    int pNode = 0;
    int sNode = 0;
    for (int q = 0; q < LT::nQ; q++)
    {
      int neighNo = grid.neighbor(q, nodeNo);
      if (nodes.getTag(neighNo) == -1 && nodes.isSolid(neighNo))
      {
        pNode = 1;
      }
      else if (nodes.getTag(neighNo) != -1 && nodes.isSolid(neighNo))
      {
        sNode = 1;
      }
    }
    if (pNode && sNode)
      tagsField(0, nodeNo) = 1;
  }

  for (auto nodeNo : bulkNodes)
  {
    nodeTypeField(0, nodeNo) = 0;
    if (tagsField(0, nodeNo) > 0)
    {
      for (int q = 0; q < LT::nQ - 1; q++)
      {
        int neighNo = grid.neighbor(q, nodeNo);
        if (tagsField(0, neighNo) > 0)
          nodeTypeField(0, nodeNo) += 1;
      }
    }
  }

  for (auto nodeNo : bulkNodes)
  {
    if (tagsField(0, nodeNo) > 0)
    {
      int numNeigh = 0;
      if (nodes.getTag(nodeNo) == 2 || nodes.getTag(nodeNo) == 4)
      {
        numNeigh = 1;
      }
      if (nodeTypeField(0, nodeNo) <= numNeigh)
      {
        for (int q = 0; q < LT::nQ - 1; q++)
        {
          int neighNo = grid.neighbor(q, nodeNo);
          if (nodes.isSolid(neighNo))
          {
            nodes.addNodeTag(-1, neighNo);
          }
        }
      }
    }
  }

  for (auto nodeNo : bulkNodes)
  {
    // tagsField(0, nodeNo) = nodes.getTag(nodeNo);

    // nodeTypeField(0, nodeNo) = nodes.getType(nodeNo);
    nodeNoField(0, nodeNo) = nodeNo;
  }

  // -------------------------------------------------------------------------------------- Checking node neighborhood 
  // field 0: solid neighbors
  // field 1: -1 tag neighbors
  ScalarField numSolidNeig(2, grid.size());
  VectorField<LT> boundaryMeanDir(2, grid.size());
  ScalarField numSolidNeigBnd(2, grid.size());
  VectorField<LT> boundaryMeanDirBnd(2, grid.size());
  for (auto nodeNo: bulkNodes) {
    for (int n=0; n<2; ++n) {
      numSolidNeig(n, nodeNo) = 0;
      boundaryMeanDir.set(n, nodeNo) = 0;
    }
    std::valarray<lbBase_t> cSum0(0.0, LT::nD);
    std::valarray<lbBase_t> cSum1(0.0, LT::nD);
    //auto neighbors = grid.neighbor(nodeNo);
    for (int q=0; q < LT::nQ; ++q) {
      const int neigNo = grid.neighbor(q, nodeNo);
      if (neigNo == 0) {
        std::cout << "ERROR: Neighbor is DEFAULT NODE!" << std::endl;
        exit(0);
      }
      if (nodes.isSolid(neigNo)) {
        numSolidNeig(0, nodeNo) += 1;
        cSum0 += LT::cValarray(q);
        if (nodes.getTag(neigNo) == -1) {
          numSolidNeig(1, nodeNo) += 1;
          cSum1 += LT::cValarray(q);
        }
      } else {
        if (nodes.getTag(neigNo) == -1) {
          std::cout << "There exists nodes that are tagged as -1 and NOT solid\n";
        }
      }
    }
    boundaryMeanDir.set(0, nodeNo) = cSum0;
    boundaryMeanDir.set(1, nodeNo) = cSum1;
  }

  

  // **********
  // OUTPUT VTK
  // **********
  Output<LT> output(grid, bulkNodes, outputDir2, myRank, nProcs);

  if (laplacePressureRun)
    output.add_file("lb_run_laplace_fluid");
  else
    output.add_file("lb_run_fluid");

  output.add_scalar_variables({"rho", "pressure", "p_perturb", "nodeType", "BndTags", "nodeNoField", "numSolidNeig", "numSolidNeigBnd"}, {rho, pressure, pPert, nodeTypeField, tagsField, nodeNoField, numSolidNeig, numSolidNeigBnd});
  output.add_vector_variables({"vel", "force", "boundaryMeanDir", "boundaryMeanDirBnd"}, {vel, force, boundaryMeanDir, boundaryMeanDirBnd});

  // *********
  // MAIN LOOP
  // *********

  for (int i = 0; i <= nIterations; i++)
  {

    // const lbBase_t ramp{ 0.5 * (1-std::cos(3.14159*std::min(i, 10000)/10000.0)) };
    const lbBase_t ramp{(0.6 + 0.4 * std::sin(2 * 3.14159 * i / PressurePeriod))};

    for (auto nodeNo : bulkNodes)
    {
      // Copy of local velocity diestirubtion
      const std::valarray<lbBase_t> fNode = f(0, nodeNo);

      // Set force

      // std::valarray<lbBase_t> forceNode = 0.0*laplaceForce(0, nodeNo);
      std::valarray<lbBase_t> forceNode(0.0, LT::nD);

      for (int n = 0; n < maxBoundaryIndicator; ++n)
      {
        forceNode += dP[n] * ramp * laplaceForce(n, nodeNo);
      }

      // const lbBase_t dimX = vtklb.getGlobaDimensions(0) ;
      // forceNode[0] = deltaPressure/dimX;
      // forceNode[0] = deltaPressure/(dimX-2);

      // forceNode[0] = deltaPressure/(dimX-3);
      // forceNode[1] = 0.0;

      force.set(0, nodeNo) = forceNode;

      // Macroscopic values
      const lbBase_t rhoNode = calcRho<LT>(fNode);
      const auto velNode = calcVel<LT>(fNode, rhoNode, forceNode);
      pressure(0, nodeNo) = rhoNode * LT::c2;
      pPert(0, nodeNo) = 0;
      for (int n = 0; n < maxBoundaryIndicator; ++n)
      {
        pressure(0, nodeNo) += (laplacePressureRun ? dP[n] * ramp * laplacePressure(n, nodeNo) : 0.0);
        pPert(0, nodeNo) += (laplacePressureRun ? dP[n] * ramp * laplacePressure(n, nodeNo) : 0.0);
      }

      // Save density and velocity for printing
      rho(0, nodeNo) = rhoNode;
      vel.set(0, nodeNo) = velNode;

      // BGK-collision term
      const lbBase_t u2 = LT::dot(velNode, velNode);
      const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
      const std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);

      // Calculate the Guo-force correction
      const lbBase_t uF = LT::dot(velNode, forceNode);
      const std::valarray<lbBase_t> cF = LT::cDotAll(forceNode);
      const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);

      // Collision and propagation
      fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);

    } // End nodes

    // Swap data_ from fTmp to f;
    f.swapData(fTmp); // LBfield

    // *******************
    // BOUNDARY CONDITIONS
    // *******************
    // Mpi
    mpiBoundary.communicateLbField(0, f, grid);
    // Half way bounce back
    // bounceBackBnd.apply(f, grid);
    // Pressure
    // std::vector <lbBase_t> rho_bnd {1.0 + dP[0], 1.0 + dP[1], 1.0 + dP[2], 1.0 + dP[3], 1.0 + dP[4], 1.0 + dP[5], 1.0 + dP[6]};
    std::vector<lbBase_t> rho_bnd(maxBoundaryIndicator, 1.0);

    // pBnd.apply(0, f, pBndNorms, rho_bnd, vel, nodes, grid);
    pBnd.applyTest(0, f, pBndNorms, rho_bnd, vel, nodes, grid, numSolidNeigBnd, boundaryMeanDirBnd);

    /*
    if (laplacePressureRun) {
              zouHePressureBoundary(boundaryNodesp1, f, 1.0, force, grid);
              zouHePressureBoundary(boundaryNodesp2, f, 1.0, force, grid);
          }
          else {
              zouHePressureBoundary(boundaryNodesp1, f, 1.0 + 0.5*deltaPressure*LT::c2Inv, force, grid);
              zouHePressureBoundary(boundaryNodesp2, f, 1.0 - 0.5*deltaPressure*LT::c2Inv, force, grid);

          }
    */

    // *************
    // WRITE TO FILE
    // *************
    if (((i % nItrWrite) == 0))
    {
      output.write(i);
      if (myRank == 0)
      {
        std::cout << "PLOT AT ITERATION : " << i << std::endl;
      }
    }

  } // End iterations

  MPI_Finalize();

  return 0;
}
