// //////////////////////////////////////////////
//
//  SUBGRID BOUNDARY CONDTIONS
//
// TO RUN PROGRAM: type "mpirun -np <#procs> badchimpp" in command
// line in main directory
//
// //////////////////////////////////////////////

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>

#include "../LBSOLVER"
#include "../IO"

#include "LBsubgridboundary.h"
#include "LBsubgridboundaryd.h"
// #include "LBboundarybasic.h"
#include "LBregularbasic.h"
#include "LBbouncebackbasic.h"
#include "LBregularboundarybasic.h"
#include "LBregularboundarybasicII.h"
#include "LBregularboundarybasicIII.h"

// SET THE LATTICE TYPE
#define LT D2Q9


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
    std::string chimpDir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

    // ***********************
    // SETUP GRID AND GEOMETRY
    // ***********************
    std::cout << " SETUP " << std::endl;
    Input input(inputDir + "input.dat");
    std::cout << " ETTER INPUT " << std::endl;
    
    LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
    Grid<LT> grid(vtklb);
    Nodes<LT> nodes(vtklb, grid);
    BndMpi<LT> mpiBoundary(vtklb, nodes, grid);


    // SETUP BULK NODES
    std::vector<int> bulkNodes = findBulkNodes(nodes);

    // *************
    // SET LB VALUES
    // *************

    // Number of iterations
    int nIterations = static_cast<int>( input["iterations"]["max"]);
    // Write interval
    int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    // Relaxation time
    lbBase_t tau = input["fluid"]["tau"][0];
    // Vector source
    VectorField<LT> bodyForce(1, 1);
    bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);
    // Driver force
    std::valarray<lbBase_t> force = bodyForce(0, 0);

    //output directory
    std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));
    std::string outDir2 = outputDir+"out"+dirNum;


    // **************
    // GEOMETRY
    // **************
    // -- New boundary condition
    

    ScalarField surfaceDistance(1, grid.size());vtklb.toAttribute("q");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        surfaceDistance(0, n) = vtklb.getScalarAttribute<double>();
    }

    VectorField<LT> surfaceNormal(1, grid.size());
    vtklb.toAttribute("nx");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        surfaceNormal(0, 0, n) = vtklb.getScalarAttribute<double>();
    }
    vtklb.toAttribute("ny");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        surfaceNormal(0, 1, n) = vtklb.getScalarAttribute<double>();
    }

    VectorField<LT> surfaceTangent(1, grid.size());
    vtklb.toAttribute("tx");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        surfaceTangent(0, 0, n) = vtklb.getScalarAttribute<double>();
    }
    vtklb.toAttribute("ty");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        surfaceTangent(0, 1, n) = vtklb.getScalarAttribute<double>();
    }


/* EJE
   ScalarField ux_analytic(1, grid.size());
    vtklb.toAttribute("u_analytic");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        ux_analytic(0, n) = vtklb.getScalarAttribute<double>();
    }

    ScalarField pi_xx_analytic(1, grid.size());
    ScalarField pi_xy_analytic(1, grid.size());
    vtklb.toAttribute("Pi_neq_xx");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        pi_xx_analytic(0, n) = vtklb.getScalarAttribute<double>();
    }
    vtklb.toAttribute("Pi_neq_xy");
    for (int n=vtklb.beginNodeNo(); n<vtklb.endNodeNo(); ++n) {
        pi_xy_analytic(0, n) = vtklb.getScalarAttribute<double>();
    } 
*/

    // ******************
    // MACROSCOPIC FIELDS
    // ******************
    // Density
    ScalarField rho(1, grid.size());
    // Velocity
    VectorField<LT> vel(1, grid.size());
    // Tensor values
    ScalarField pi_xx(1, grid.size());
    ScalarField pi_xy(1, grid.size());
    ScalarField pi_yy(1, grid.size());
    // Initiate values
    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 1.0;
        pi_xx(0, nodeNo) = 0;// EJE pi_xx_analytic(0, nodeNo);
        pi_xy(0, nodeNo) = 0;// EJE pi_xy_analytic(0, nodeNo);
        pi_yy(0, nodeNo) = 0;

        // EJE vel(0,0,nodeNo) = ux_analytic(0, nodeNo);
        vel(0,0,nodeNo) = 0.01;
        for (int d=1; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }

    // ******************
    // SETUP BOUNDARY
    // ******************
    // WORKING
//    OneNodeSubGridBndDyn<LT> fluidWallBnd(findFluidBndNodes(nodes), nodes, grid, surfaceDistance, surfaceNormal, surfaceTangent, rho, bodyForce, tau);


//    OneNodeSubGridBnd<LT> fluidWallBnd(findFluidBndNodes(nodes), nodes, grid, qAttribute, surfaceNormal, surfaceTangent, rho, bodyForce, tau);
//    HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);
//    BounceBackBasic<LT> bndBasic(findFluidBndNodes(nodes), nodes, grid);
    RegularBoundaryBasicIII<LT> fluidNoSlipBnd(findFluidBndNodes(nodes), surfaceDistance, surfaceNormal, surfaceTangent, rho, bodyForce, tau, nodes, grid);
    ScalarField boundaryPlot(1, grid.size());
    
    for (int nodeNo = 0; nodeNo < grid.size(); ++nodeNo) {
        boundaryPlot(0, nodeNo) = 0.0;
    }

    // *********
    // LB FIELDS
    // *********
    LbField<LT> f(1, grid.size());  // LBfield
    LbField<LT> fTmp(1, grid.size());  // LBfield
    // initieate values

    for (auto nodeNo: bulkNodes) {
        auto cF = LT::cDotAll(bodyForce(0,0));
        auto cu = LT::cDotAll(vel(0, nodeNo));
        auto uu = vel(0, 0, nodeNo)*vel(0, 0, nodeNo);
        for (int alpha = 0; alpha < LT::nQ; ++alpha) {
            auto c = LT::c(alpha);
            f(0, alpha, nodeNo)  = LT::w[alpha] * rho(0, nodeNo);  
            f(0, alpha, nodeNo) += LT::w[alpha] * LT::c2Inv * rho(0, nodeNo) * cu[alpha];
            f(0, alpha, nodeNo) += LT::w[alpha] * LT::c4Inv0_5 *(cu[alpha] * cu[alpha] - LT::c2*uu) * rho(0, nodeNo);
// EJE            f(0, alpha, nodeNo) +=     LT::w[alpha] * LT::c4Inv0_5 *(c[0] * c[0] - LT::c2) * pi_xx_analytic(0, nodeNo);            
// EJE            f(0, alpha, nodeNo) += 2 * LT::w[alpha] * LT::c4Inv0_5 *c[0] * c[1] * pi_xy_analytic(0, nodeNo);
            f(0, alpha, nodeNo) -= 0.5*LT::c2Inv*LT::w[alpha]*cF[alpha];
        }
    }

    // **********
    // OUTPUT VTK
    // **********
    auto node_pos = grid.getNodePos(bulkNodes); // Need a named variable as Outputs constructor takes a reference as input
    auto global_dimensions = vtklb.getGlobaDimensions();
    // Setup output file
    Output output(global_dimensions, outDir2, myRank, nProcs, node_pos);
    output.add_file("lb_run");
    output.add_file("boundary");
    // Add density to the output file
    VectorField<D3Q19> velIO(1, grid.size());
    VectorField<D3Q19> velAnalyticIO(1, grid.size());    
    output["lb_run"].add_variable("rho", rho.get_data(), rho.get_field_index(0, bulkNodes), 1);
    output["lb_run"].add_variable("Pi_xx", pi_xx.get_data(), pi_xx.get_field_index(0, bulkNodes), 1);
    output["lb_run"].add_variable("Pi_xy", pi_xy.get_data(), pi_xy.get_field_index(0, bulkNodes), 1);
    output["lb_run"].add_variable("Pi_yy", pi_yy.get_data(), pi_yy.get_field_index(0, bulkNodes), 1);
    //output["lb_run"].add_variable("vel", velIO.get_data(), vel.get_field_index(0, bulkNodes), LT::nD);
    output["lb_run"].add_variable("vel", velIO.get_data(), velIO.get_field_index(0, bulkNodes), 3);
    // EJE output["lb_run"].add_variable("vel_analytic", velAnalyticIO.get_data(), velAnalyticIO.get_field_index(0, bulkNodes), 3);
    // Print geometry 
    outputGeometry("lb_geo", outDir2, myRank, nProcs, nodes, grid, vtklb);
    // Write error
    // std::ofstream fileOutputError(outDir2 + "/ux_error.dat");  // EJE
    
    // *********
    // MAIN LOOP
    // *********
    // For all time steps
    const std::clock_t beginTime = std::clock();
    for (int i = 0; i <= nIterations; i++) {
        // ***************
        // GLOBAL COUNTERS
        // ***************
        double rhoSumLocal = 0;
        // For all bulk nodes
        for (auto nodeNo: bulkNodes) {
            // Copy of local velocity diestirubtion
            std::valarray<lbBase_t> fNode = f(0, nodeNo);

            // MACROSCOPIC VALUES
            lbBase_t rhoNode = calcRho<LT>(fNode);
            // std::valarray<lbBase_t> velNode = calcVel<LT>(fNode, rhoNode, force);
            auto velNode = calcVel<LT>(fNode, rhoNode, force);
            // velNode = calcVel(fNode, rhoNode, force); // Kan bruke using

            rho(0, nodeNo) = rhoNode - 1;
            rhoSumLocal += rhoNode-1;
            for (int d = 0; d < LT::nD; ++d)
                vel(0, d, nodeNo) = velNode[d];


            // BGK COLLISION TERM
            // SRT
            lbBase_t u2 = LT::dot(velNode, velNode);
            std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
            std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);
            lbBase_t uF = LT::dot(velNode, force);
            std::valarray<lbBase_t> cF = LT::cDotAll(force);
            std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);

            // COLLISION AND PROPAGATION
            fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);

        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(0, f, grid);
        // Half way bounce back
        // WORKING
        // fluidWallBnd.applyNonSlip(0, f, nodes, grid, surfaceDistance, surfaceNormal, surfaceTangent, bodyForce, tau);

        
        fluidNoSlipBnd.apply(0, f, surfaceDistance, surfaceNormal, surfaceTangent, bodyForce, tau, nodes, grid);
//        fluidWallBnd.apply(0, f, nodes, grid);
//        bounceBackBnd.apply(0, f, grid);
//        bndBasic.apply(0, f, grid);
        // *************
        // WRITE TO FILE
        // *************
        double rhoSumGlobal;
        MPI_Allreduce(&rhoSumLocal, &rhoSumGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        // Write error to tile
        

       /* EJE if ( (((i) % (nItrWrite/100)) == 0) && (i > 0) ) 
        {    

            lbBase_t u_error = 0;
            lbBase_t pi_error = 0;
            for ( auto nodeNo: bulkNodes ) 
            {
                // Error in velocity
                auto u = vel(0, nodeNo);
                auto u2 = LT::dot(u, u);
                u_error += std::abs( u[0] - ux_analytic(0, nodeNo) );
                for ( int d = 1; d < LT::nD; ++d)
                {
                    u_error += std::abs( u[d] );
                }
                // Error in strain rate tensor
                auto rhoNode = rho(0, nodeNo) + 1;
                auto cu = LT::cDotAll(u);
                auto fNode = f(0, nodeNo);
                lbBase_t pi[LT::nD][LT::nD]; 
                pi[0][0] = 0;
                pi[0][1] = 0;
                pi[1][0] = 0;
                pi[1][1] = 0;

                for (int alpha = 0; alpha < LT::nQ; ++alpha)
                {
                    auto cu_alpha = cu[alpha];
                    auto w = LT::w[alpha];
                    auto f_eq_alpha = w * rhoNode * ( 1 + LT::c2Inv*cu_alpha + LT::c4Inv0_5 * (cu_alpha*cu_alpha - LT::c2 * u2) );
                    auto f_alpha = fNode[alpha];
                    auto c = LT::c(alpha);
                    for (int ii = 0; ii < LT::nD; ++ii)
                    {
                        auto ci = c[ii];
                        pi[ii][ii] += (ci*ci - LT::c2) * (f_alpha - f_eq_alpha);
                        for (int jj = ii+1; jj < LT::nD; ++jj)
                        {
                            auto cj = c[jj];
                            pi[ii][jj] += ci * cj * (f_alpha - f_eq_alpha);
                        }
                    }
                }
                pi_xx(0, nodeNo) = pi[0][0];
                pi_xy(0, nodeNo) = pi[0][1];
                pi_yy(0, nodeNo) = pi[1][1];

            }
            u_error /= bulkNodes.size();
            fileOutputError << i << " " << u_error << "\n";

        } */
        
        if ( ((i % nItrWrite) == 0) && (i > 0) ) {
            for (auto nn: bulkNodes) {
                velIO(0, 0, nn) = vel(0, 0, nn);
                velIO(0, 1, nn) = vel(0, 1, nn);
                velIO(0, 2, nn) = 0;
// EJE                velAnalyticIO(0, 0, nn) = ux_analytic(0, nn);
// EJE                velAnalyticIO(0, 1, nn) = 0;
// EJE                velAnalyticIO(0, 2, nn) = 0;
            }
            // Stress tensor
            
            output.write("lb_run", i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << " ( " << float( std::clock () - beginTime ) /  CLOCKS_PER_SEC << " sec)" << std::endl;
                std::cout << "Error in mass:" << rhoSumGlobal << "  " << std::endl;
                // Setup plot over line
                int nx = 1;
                for (int ny = 1; ny < 19; ++ny ) {
                    std::vector<int> pos { nx, ny};
                    int nodeNo = grid.nodeNo(pos);
//                   std::cout << "vel = " << vel(0, 0, nodeNo) << ", " << vel(0, 1, nodeNo) << std::endl;
                }
            }
        }

    } // End iterations
    // EJE fileOutputError.close();

    MPI_Finalize();

    return 0;
}
