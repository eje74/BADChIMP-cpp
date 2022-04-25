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

//                                SET THE LATTICE TYPE
//------------------------------------------------------------------------------------- SET THE LATTICE TYPE
#define LT D2Q9

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
    //std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

    
    //=====================================================================================
    //
    //                             SETUP GRID AND GEOMETRY
    //
    //=====================================================================================

    Input input(inputDir + "input.dat");
    LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
    Grid<LT> grid(vtklb);
    Nodes<LT> nodes(vtklb, grid);
    BndMpi<LT> mpiBoundary(vtklb, nodes, grid);

    //                                 Set bulk nodes
    //------------------------------------------------------------------------------------- Set bulk nodes
    std::vector<int> bulkNodes = findBulkNodes(nodes);

    
    //                                 Set solid boundary
    //------------------------------------------------------------------------------------- Set solid boundary
    //std::vector<int> solidBoundaryNodes = findSolidBndNodes(nodes);
  
  
    //=====================================================================================
    //
    //                                 READ FROM INPUT
    //
    //=====================================================================================
    
    
    //                               Number of iterations
    //------------------------------------------------------------------------------------- Number of iterations
    // int nIterations = static_cast<int>( input["iterations"]["max"]);
    int nIterations = input["iterations"]["max"];

    //                                  Write interval
    //------------------------------------------------------------------------------------- Write interval
    // int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    int nItrWrite = input["iterations"]["write"];

    //                                  Relaxation time
    //------------------------------------------------------------------------------------- Relaxation time                          // FJERNE?
    lbBase_t tau = LT::c2Inv * input["fluid"]["viscosity"] + 0.5;

    //                                    Body force
    //------------------------------------------------------------------------------------- Body force
    VectorField<LT> bodyForce(1, 1);
    // bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);
    bodyForce.set(0, 0) = input["fluid"]["bodyforce"];

    //                                    Output directory number
    //------------------------------------------------------------------------------------- Output directory number
    std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));  
    std::string outputDir2 = outputDir + "/out" + dirNum;

    
    //=====================================================================================
    //
    //                                 DEFINE RHEOLOGY 
    //
    //=====================================================================================

    
    
    Rans<LT> rans(input);
    
  
    //=====================================================================================
    //
    //                               MACROSCOPIC FIELDS 
    //
    //=====================================================================================

    //                                   Viscosity
    //------------------------------------------------------------------------------------- Viscosity
    ScalarField viscosity(1, grid.size());

    //                                    Density
    //------------------------------------------------------------------------------------- Density
    ScalarField rho(1, grid.size());

    //                           Initiate density from file
    //------------------------------------------------------------------------------------- Initiate density from file
    vtklb.toAttribute("init_rho");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        rho(0, n) = vtklb.getScalarAttribute<lbBase_t>();
    }
    
    //                                   Velocity
    //------------------------------------------------------------------------------------- Velocity
    VectorField<LT> vel(1, grid.size());

    //                               Initiate velocity
    //------------------------------------------------------------------------------------- Initiate velocity
    for (auto nodeNo: bulkNodes) {
        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }

    //                                   rhoK field
    //------------------------------------------------------------------------------------- rhoK field
    ScalarField rhoK(1, grid.size());
    
    //                         Initiate rhoK field from file
    //------------------------------------------------------------------------------------- Initiate rhoK field from file
    vtklb.toAttribute("init_rhoK");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        rhoK(0, n) = vtklb.getScalarAttribute<lbBase_t>();
    }
    
    //                                 rhoEpsilon field
    //------------------------------------------------------------------------------------- rhoEpsilon field
    ScalarField rhoEpsilon(1, grid.size());
    
    //                        Initiate rhoEpsilon field from file
    //------------------------------------------------------------------------------------- Initiate rhoEpsilon field from file
    vtklb.toAttribute("init_rhoEpsilon");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        rhoEpsilon(0, n) = vtklb.getScalarAttribute<lbBase_t>();
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
    for (auto nodeNo: bulkNodes) {
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = LT::w[q]*rho(0, nodeNo);
	    g(0, q, nodeNo) = LT::w[q]*rhoK(0, nodeNo);
	    h(0, q, nodeNo) = LT::w[q]*rhoEpsilon(0, nodeNo);
        }
    }

    //=====================================================================================
    //
    //                                      OUTPUT VTK
    //
    //=====================================================================================

    Output<LT> output(grid, bulkNodes, outputDir2, myRank, nProcs); 
    output.add_file("lb_run");
    output.add_scalar_variables({"viscosity", "rho", "rhoK", "rhoEpsilon"}, {viscosity, rho, rhoK, rhoEpsilon});
    output.add_vector_variables({"vel"}, {vel});

    Output<LT,int> geo(grid.pos(), outputDir2, myRank, nProcs, "geo", nodes.geo(grid, vtklb));
    geo.write();
    

  
    //=====================================================================================
    //
    //                                      MAIN LOOP
    //
    //=====================================================================================

    //                              Start time step iterations
    //------------------------------------------------------------------------------------- Start time step iterations
    for (int i = 0; i <= nIterations; i++) {
      
      //                        Calculate macroscopic values for all nodes
      //------------------------------------------------------------------------------------- Calculate macroscopic values for all nodes
      for (auto nodeNo: bulkNodes) {   
      }

      //                               Communicate rho fields
      //------------------------------------------------------------------------------------- Communicate rho fields
      //(FILL IN COMMUNICATION HERE)

	
        for (auto nodeNo: bulkNodes) {
            //                           Copy of local LB distribution
	    //------------------------------------------------------------------------------------- Copy of local LB distribution
            const auto fNode = f(0, nodeNo);
	    const auto gNode = g(0, nodeNo);
	    const auto hNode = h(0, nodeNo);

	    //                                    Macroscopic values
	    //------------------------------------------------------------------------------------- Macroscopic values
	    const lbBase_t rhoNode = calcRho<LT>(fNode);
            const auto velNode = calcVel<LT>(fNode, rhoNode, bodyForce(0, 0));
	    const lbBase_t u2 = LT::dot(velNode, velNode);
	    const auto cu = LT::cDotAll(velNode);

	    rans.apply(fNode, rhoNode, velNode, u2, cu, bodyForce(0, 0), 0.0, gNode, hNode);

	    const lbBase_t rhoKNode = rans.rhoK();
	    const lbBase_t rhoENode = rans.rhoE();

	    //                         Save density and velocity for printing
	    //------------------------------------------------------------------------------------- Save density and velocity for printing
	    rho(0, nodeNo) = rhoNode;
            vel.set(0, nodeNo) = velNode;
	    rhoK(0, nodeNo) = rhoKNode;
	    rhoEpsilon(0, nodeNo) = rhoENode;
	    viscosity(0, nodeNo) = LT::c2*(tau-0.5);
	    








	    
            //                                    Macroscopic values
	    //------------------------------------------------------------------------------------- Macroscopic values
            //const lbBase_t rhoNode = calcRho<LT>(fNode);
            //const auto velNode = calcVel<LT>(fNode, rhoNode, bodyForce(0, 0));
	    //const lbBase_t rhoKNode = calcRho<LT>(gNode);
	    //const lbBase_t rhoEpsilonNode = calcRho<LT>(hNode);

            //                         Save density and velocity for printing
	    //------------------------------------------------------------------------------------- Save density and velocity for printing
            //rho(0, nodeNo) = rhoNode;
            //vel.set(0, nodeNo) = velNode;
	    //rhoK(0, nodeNo) = rhoKNode;
	    //rhoEpsilon(0, nodeNo) = rhoEpsilonNode;
                
            //                                    BGK-collision term
	    //------------------------------------------------------------------------------------- BGK-collision term
            //const lbBase_t u2 = LT::dot(velNode, velNode);
            //const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);

	    //auto omegaBGK = rans.omegaBGK(fNode, rhoNode, velNode, u2, cu, bodyForce(0, 0), 0);

            

	    //                                    LB interaction Terms
	    //------------------------------------------------------------------------------------- LB interaction Terms
            const lbBase_t uF = LT::dot(velNode, bodyForce(0, 0));
            const auto cF = LT::cDotAll(bodyForce(0, 0));
            tau = 1.0; //rans.tau();

	    
	    std::cout<<"Force = ("<< bodyForce(0, 0, 0) << ", " << bodyForce(0, 0, 1) << ")" <<std::endl;

	    //------------------------------------------------------------------------------------- Omegas: Flow 
	    const auto deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);
	    const auto feqNode = calcfeq<LT>(rhoNode, u2, cu);
	    const auto omegaBGK = calcOmegaBGK_TEST<LT>(fNode, feqNode, tau);

	    //------------------------------------------------------------------------------------- Omegas: K & E
	    const lbBase_t tauKNode =  1.0;//rans.tauK();
	    const lbBase_t tauENode =  1.0;//rans.tauE();
	    
	    const auto dOmegaFK = calcDeltaOmegaFDiff<LT>(tauKNode, rhoKNode/rhoNode, cu, uF, cF);
	    const auto dOmegaFE = calcDeltaOmegaFDiff<LT>(tauENode, rhoENode/rhoNode, cu, uF, cF);

	    const auto dOmegaSourceK = calcDeltaOmegaR<LT>(tauKNode, cu, rans.sourceK());
	    const auto dOmegaSourceE = calcDeltaOmegaR<LT>(tauKNode, cu, rans.sourceE());

	    const auto geqNode = calcfeq<LT>(rhoKNode, u2, cu);
	    const auto omegaBGK_K = calcOmegaBGK_TEST<LT>(gNode, geqNode, tauKNode);

	    const auto heqNode = calcfeq<LT>(rhoENode, u2, cu);
	    const auto omegaBGK_E = calcOmegaBGK_TEST<LT>(hNode, heqNode, tauENode);
	    
	    
	    //                               Collision and propagation
	    //------------------------------------------------------------------------------------- Collision and propagation
            fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);
	    gTmp.propagateTo(0, nodeNo, gNode + omegaBGK_K + dOmegaFK + dOmegaSourceK, grid);
	    hTmp.propagateTo(0, nodeNo, hNode + omegaBGK_E + dOmegaFE + dOmegaSourceE, grid);

        }//------------------------------------------------------------------------------------- End bulkNodes

	
        //                                   Swap data_ from fTmp to f etc.
	//------------------------------------------------------------------------------------- Swap data_ from fTmp to f etc.
        f.swapData(fTmp);  // flow LBfield
	g.swapData(gTmp);  // rhoK LBfield
	h.swapData(hTmp);  // rhoEpsilon LBfield
	

 
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
        bounceBackBnd.apply(f, grid);
	bounceBackBnd.apply(g, grid);
	bounceBackBnd.apply(h, grid);

	//=====================================================================================
	//
	//                                      WRITE TO FILE
	//
	//=====================================================================================
        if ( ((i % nItrWrite) == 0)  ) {
            
            output.write(i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
            }
        }

    }//------------------------------------------------------------------------------------- End time step iterations
    

    MPI_Finalize();

    return 0;
}
