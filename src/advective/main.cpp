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

// SET THE LATTICE TYPE
#define LT D2Q9
template <typename T>
using Output = LBOutputUnstructured<LT, T, VTK::BINARY, VTK::voxel>;
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
    //std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
    std::string chimpDir = "./";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

    // ***********************
    // SETUP GRID AND GEOMETRY
    // ***********************
    Input input(inputDir + "input.dat");
    LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
    Grid<LT> grid(vtklb);
    Nodes<LT> nodes(vtklb, grid);
    BndMpi<LT> mpiBoundary(vtklb, nodes, grid);
    // Set bulk nodes
    std::vector<int> bulkNodes = findBulkNodes(nodes);
    
    // *************
    // READ FROM INPUT
    // *************
    // Number of iterations
    // int nIterations = static_cast<int>( input["iterations"]["max"]);
    int nIterations = input["iterations"]["max"];
    // Write interval
    // int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    int nItrWrite = input["iterations"]["write"];
    // Relaxation time
    lbBase_t tau = input["fluid"]["viscosity"]*LT::c2Inv + 0.5;
    // Body force
    VectorField<LT> bodyForceInit(1, 1);
    // bodyForceInit.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);
    bodyForceInit.set(0, 0) = input["fluid"]["bodyforce"];

    lbBase_t const momXInit = input["fluid"]["momx"];
    // std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));
    std::string dirNum = input["out"]["directoryNum"];
    std::string outputDir2 = outputDir + "/out" + dirNum;
    
    // *************
    // DEFINE RHEOLOGY
    // *************
    Newtonian<LT> newtonian(tau);
    
    // ******************
    // MACROSCOPIC FIELDS
    // ******************
    // Viscosity
    ScalarField viscosity(1, grid.size());
    // gammaDot
    ScalarField gammaDot(1, grid.size());
    // epsilonDot
    ScalarField epsilonDot(1, grid.size());
    // E00
    ScalarField E00(1, grid.size());
    // E01
    ScalarField E01(1, grid.size());
    // E00_2
    ScalarField E00_2(1, grid.size());
    // Density
    ScalarField rho(1, grid.size());
    // Tracer field
    ScalarField phi(2, grid.size());
    // Initiate density from file
    vtklb.toAttribute("init_rho");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        rho(0, n) = vtklb.getScalarAttribute<lbBase_t>();
    }
    // Initiate phi from file
    vtklb.toAttribute("init_phi");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        phi(0, n) = vtklb.getScalarAttribute<lbBase_t>();
    }
    // Initiate force from file
    VectorField<LT> bodyForce(1, grid.size());
    vtklb.toAttribute("force");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
      bodyForce(0, 0, n) = vtklb.getScalarAttribute<lbBase_t>();
      bodyForce(0, 0, n) *= bodyForceInit(0, 0, 0);
      bodyForce(0, 1, n) = bodyForceInit(0, 1, 0);
    }
    // Velocity
    VectorField<LT> vel(1, grid.size());
    // Initiate velocity
    for (auto nodeNo: bulkNodes) {
        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
	vel(0, 0, nodeNo) = momXInit;
    }

    // ******************
    // SETUP BOUNDARY
    // ******************
    HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

    // *********
    // LB FIELDS
    // *********
    LbField<LT> f(1, grid.size()); 
    LbField<LT> fTmp(1, grid.size());
    LbField<LT> g(2, grid.size()); 
    LbField<LT> gTmp(2, grid.size());
    // initiate lb distributions
    for (auto nodeNo: bulkNodes) {
      auto rhoNode = rho(0, nodeNo);
      auto velNode = vel(0, nodeNo);
      const lbBase_t u2 = LT::dot(velNode, velNode);
      const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
        for (int q = 0; q < LT::nQ; ++q) {
          f(0, q, nodeNo) = rhoNode*LT::w[q]*(1 + (LT::c2Inv*cu[q] + LT::c4Inv0_5*(cu[q]*cu[q] - LT::c2*u2)) ) ;
          g(0, q, nodeNo) = phi(0, nodeNo)*f(0, q, nodeNo);
          g(1, q, nodeNo) = phi(0, nodeNo)*rhoNode*LT::w[q]*(1 + (LT::c2Inv*cu[q] + LT::c4Inv0_5*(cu[q]*cu[q] - LT::c2*u2)) -0.5*(LT::cNorm[q]*LT::cNorm[q]*LT::c2Inv - LT::nD));
        }
    }

    // **********
    // OUTPUT VTK
    // **********
    Output<double> output(grid, bulkNodes, outputDir2, myRank, nProcs);
    output.add_file("lb_run");
    output.add_variables({"rho", "phi", "vel", "viscosity", "gammaDot", "epsilonDot", "E00", "E01", "E00_2"},
                         { rho,   phi,   vel,   viscosity,   gammaDot,   epsilonDot,   E00,   E01,   E00_2});

    // VTK::Output<VTK_CELL, double> output(VTK::BINARY, grid.getNodePos(bulkNodes), outputDir2, myRank, nProcs);
    // output.add_file("lb_run");
    // output.add_variable("rho", 1, rho.get_data(), rho.get_field_index(0, bulkNodes));
    // output.add_variable("phi", 1, phi.get_data(), phi.get_field_index(0, bulkNodes));
    // output.add_variable("phi1", 1, phi.get_data(), phi.get_field_index(1, bulkNodes));
    // output.add_variable("vel", LT::nD, vel.get_data(), vel.get_field_index(0, bulkNodes));
    // output.add_variable("viscosity", 1, viscosity.get_data(), viscosity.get_field_index(0, bulkNodes));
    // output.add_variable("gammaDot", 1, gammaDot.get_data(), gammaDot.get_field_index(0, bulkNodes));
    // output.add_variable("epsilonDot", 1, epsilonDot.get_data(), epsilonDot.get_field_index(0, bulkNodes));
    // output.add_variable("E00", 1, E00.get_data(), E00.get_field_index(0, bulkNodes));
    // output.add_variable("E01", 1, E01.get_data(), E01.get_field_index(0, bulkNodes));
    // output.add_variable("E00_2", 1, E00_2.get_data(), E00_2.get_field_index(0, bulkNodes));

    
    // *********
    // MAIN LOOP
    // *********
    for (int i = 0; i <= nIterations; i++) {
        // Calculate macroscopic values for all nodes
      /*
      for (auto nodeNo: bulkNodes) {
            
        }
      */
        // Communicate rho fields
        
        for (auto nodeNo: bulkNodes) {
            // Copy of local velocity diestirubtion
            const std::valarray<lbBase_t> fNode = f(0, nodeNo);
	    const std::valarray<lbBase_t> gNode = g(0, nodeNo);
	    const std::valarray<lbBase_t> g1Node = g(1, nodeNo);
	    
            // Macroscopic values
            const lbBase_t rhoNode = calcRho<LT>(fNode);
            const auto velNode = calcVel<LT>(fNode, rhoNode, bodyForce(0, nodeNo));

	    const lbBase_t phiNode = calcRho<LT>(gNode);
	    const auto MiNode = LT::qSumC(gNode);

	    const lbBase_t phi1Node = calcRho<LT>(g1Node);

            // Save density and velocity for printing
            rho(0, nodeNo) = rhoNode;
            vel.set(0, nodeNo) = velNode;
	    phi(0, nodeNo) = phiNode;
	    phi(1, nodeNo) = phi1Node;
	    
            // BGK-collision term
            const lbBase_t u2 = LT::dot(velNode, velNode);
            const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
            auto omegaBGK = newtonian.omegaBGK(fNode, rhoNode, velNode, u2, cu, bodyForce(0, nodeNo), 0);

            // Calculate the Guo-force correction
            const lbBase_t uF = LT::dot(velNode, bodyForce(0, nodeNo));
            const std::valarray<lbBase_t> cF = LT::cDotAll(bodyForce(0, nodeNo));
            tau = newtonian.tau();
            viscosity(0, nodeNo) = newtonian.viscosity();
	    gammaDot(0, nodeNo) = newtonian.gammaDot();
	    epsilonDot(0, nodeNo) = newtonian.epsilonDot();
	    E00(0, nodeNo) = newtonian.E00();
	    E01(0, nodeNo) = newtonian.E01();
	    
	    E00_2(0, nodeNo) = newtonian.E00_2();
            const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);

	    
	    const std::valarray<lbBase_t> omegaAdvMom = LT::cDotAll(phiNode*velNode - MiNode);
	    std::valarray<lbBase_t> omegaPhi(LT::nQ);
	    std::valarray<lbBase_t> omegaPhi1(LT::nQ);
	    for (int q = 0; q < LT::nQ; ++q){
	      omegaPhi[q]=2*LT::w[q]*LT::c2Inv*omegaAdvMom[q];
	      omegaPhi1[q]=-/*1.999**/(g1Node[q]-phi1Node*rhoNode*LT::w[q]*(1 + (LT::c2Inv*cu[q] + LT::c4Inv0_5*(cu[q]*cu[q] - LT::c2*u2)) -0.5*(LT::cNorm[q]*LT::cNorm[q]*LT::c2Inv - LT::nD)));
	    }
		
            // Collision and propagation
            fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);
	    gTmp.propagateTo(0, nodeNo, gNode + omegaPhi, grid);
	    gTmp.propagateTo(1, nodeNo, g1Node + omegaPhi1, grid);

        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield
	g.swapData(gTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(0, f, grid);
	mpiBoundary.communicateLbField(g, grid);
	
        // Half way bounce back
        bounceBackBnd.apply(f, grid);
	bounceBackBnd.apply(g, grid);
        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0)  ) {
	  /*
	  for (auto nn: bulkNodes) {
                velIO(0, 0, nn) = vel(0, 0, nn);
                velIO(0, 1, nn) = vel(0, 1, nn);
                velIO(0, 2, nn) = 0;
            }
            output.write("lb_run", i);
	  */
	  output.write(i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
            }
        }

    } // End iterations

    MPI_Finalize();

    return 0;
}
