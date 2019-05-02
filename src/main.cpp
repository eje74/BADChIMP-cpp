// //////////////////////////////////////////////
//
// TWO PHASE SOLVER
//
//
// Using the color gradient method (Rothman-Keller
//  type) with the Reis Phillips surface pertubation
//  and / Latva-Kokko recolor step. (As suggested
//  by Leclaire et al and Yu-Hang Fu et al)
//
// This is an explicit two-phase implementation. A
// full multiphase implementation will be added later.
//
// TO RUN PROGRAM: type bin/runner on command
// line in main directory
//
// //////////////////////////////////////////////

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "LBlatticetypes.h"
#include "LBgrid.h"
#include "LBfield.h"
#include "LBboundary.h"
#include "LBhalfwaybb.h"
#include "LBbulk.h"
#include "LBgeometry.h"
#include "LBmacroscopic.h"
#include "LBcollision.h"
#include "LBiteration.h"
#include "LBinitiatefield.h"
#include "LBsnippets.h"

#include "LBbndmpi.h"

#include "Input.h"
#include "Output.h"
#include "Mpi_class.h"
#include "Geo.h"
#include "Field.h"


// SET THE LATTICE TYPE
// #define LT D2Q9
#define LT D3Q19




int main()
// JLV int main(int argc, char *argv[])
{
    std::cout << "Begin test Two phase new" << std::endl;

    std::string mpiDir = "/home/ejette/Programs/GITHUB/badchimpp/input/mpi/";
    std::string inputDir = "/home/ejette/Programs/GITHUB/badchimpp/input/";

    // std::string mpiDir = "/home/olau/Programs/Git/BADChIMP-cpp/input/mpi/";
    // std::string inputDir = "/home/olau/Programs/Git/BADChIMP-cpp/input/";

    // read input files
    //Input input("input.dat"); //input.print();
    Input input(inputDir + "input.dat");

    // initialize MPI
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // JLV
    //    Mpi mpi; //(&argc, &argv, input["mpi"]["procs"]);
    //    mpi.start(&argc, &argv, "node_labels.mpi", "rank.mpi");
    //    mpi.print();
    // JLV



    // READ BOUNDARY FILES WITH MPI INFORMATION
    MpiFile<LT> rankFile(mpiDir + "rank.mpi");
    MpiFile<LT> localFile(mpiDir + "rank_" + std::to_string(myRank) + "_labels.mpi");
    MpiFile<LT> globalFile(mpiDir + "node_labels.mpi");


    // SETUP GRID
    Grid<LT> grid  = Grid<LT>::makeObject(localFile, rankFile);


    // SETUP MPI BOUNDARY
    BndMpi<LT> mpiBoundary(myRank);

    mpiBoundary.setupBndMpi(localFile, globalFile, rankFile, grid);


    // SETUP BOUNCE BACK BOUNDARY (fluid boundary)
    HalfWayBounceBack<LT> bbBnd = makeFluidBoundary<HalfWayBounceBack>(myRank, grid);

    // SETUP SOLID BOUNDARY
    std::vector<int> solidBnd = findSolidBndNodes(myRank, grid);

    // SETUP BULK NODES
    std::vector<int> bulkNodes = findBulkNodes(myRank, grid);

    // SETUP CONST PRESSURE AND FLUID SINK
    std::vector<int> constDensNodes, sourceNodes;
    for (auto bulkNode: bulkNodes) {
        int y = grid.pos(bulkNode, 1);
        if (grid.getRank(bulkNode) == myRank) {
            if ( y < 3) // sink
                sourceNodes.push_back(bulkNode);
            if ( y > (rankFile.dim(1) - 2*1 - 4) ) // Const pressure
                constDensNodes.push_back(bulkNode);
        }
    }


    // READ INPUT FILE
    // lbBase_t force[3] = {0.0, 0.0, 0.}; //Denne burde vi få fra input

    // Scalar source
    ScalarField Q(2, grid.size());
    for (int n = 0; n < Q.size(); ++n) {
        Q(0, n) = 0.0;
        Q(1, n) = 0.0;
    }

    // Vector source
    VectorField<LT> bodyForce(1, 1);
    std::vector<lbBase_t> tmpVec = input["fluid"]["bodyforce"];
    for (int d = 0; d < LT::nD; ++d)
        bodyForce(0, d, 0) = tmpVec[static_cast<std::size_t>(d)];

    int nIterations = static_cast<int>( input["iterations"]["max"]);

    lbBase_t tau0 = input["fluid"]["tau"][0];
    lbBase_t tau1 = input["fluid"]["tau"][1];
    lbBase_t sigma = input["fluid"]["sigma"];
    lbBase_t beta = input["fluid"]["beta"];

    // SET DERIVED VARAIBLES
    lbBase_t nu0Inv = 1.0 / (LT::c2 * (tau0 - 0.5));
    lbBase_t nu1Inv = 1.0 / (LT::c2 * (tau1 - 0.5));

    // SETUP LB FIELDS
    LbField<LT> f(2, grid.size());  // LBfield
    LbField<LT> fTmp(2, grid.size());  // LBfield
    // STL
    //LB_field f_(2, nNodes, LT::nQ);
    //LB_field fTmp_(2, nNodes, LT::nQ);

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(2, grid.size()); // LBfield
    VectorField<LT> vel(1, grid.size()); // LBfield
    ScalarField cgField(1, grid.size()); // LBfield
    // STL
    //Scalar_field rho(2, nNodes); // LBfield
    //Vector_field vel_(1, nNodes, LT::nD); // LBfield
    //Scalar_field cgField_(1, nNodes); // LBfield



    // FILL MACROSCOPIC FIELDS
    //   Fluid densities and velocity
    std::srand(8549388);
    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 0.5; // 1.0*(1.0 * std::rand()) / (RAND_MAX * 1.0);
        rho(1, nodeNo) = 1 - rho(0, nodeNo);

        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }
    //   Solid boundary (Wettability)
    for (auto nodeNo: solidBnd) {
        rho(0, nodeNo) = 1.0;
        rho(1, nodeNo) = 1.0 - rho(0, nodeNo);
    }


    // INITIATE LB FIELDS
    // -- phase 0
    //f(0,bulk).initializeLB(rho(0,bulk), vel(0,bulk))
    initiateLbField(0, 0, 0, bulkNodes, rho, vel, f);  // LBinitiatefield
    // -- phase 1
    initiateLbField(1, 1, 0, bulkNodes, rho, vel, f);  // LBinitiatefield



    //    Output output("out", myRank, nProcs-1, geo2, 0);
//    output.add_file("fluid");
    // 0+1 to skip first dummy node
//    output["fluid"].add_variables({"rho0","rho1"}, {&rho(0,0),&rho(1,0)}, {sizeof(rho(0,0)),sizeof(rho(1,0))}, {1,1}, {rho.num_fields(),rho.num_fields()});
    //output.set_time(0);
    //output.write("fluid","chem");
    //output.write("all");
//    output.write("fluid",0);



    // -----------------MAIN LOOP------------------
    /* Comments to main loop:
     * Calculation of cu is kept outside of calcOmega and calcDeltaOmega
     *  to avoid recalculation  of the dot-product.
     * It is therefore natural to give all scalar products as inputs to the
     *  calcOmega and calcDeltaOmega.
     * Note that we have a oneIteration function in 'LBiteration.h'. This seemed
     *   to affect the speed an gave a slow down of 10-15 %.
     */

    for (int i = 0; i < nIterations; i++) {

        for (auto nodeNo : bulkNodes) {  // Change nElements to nNodes?
            // UPDATE MACROSCOPIC DENSITIES
            lbBase_t rho0Node, rho1Node;
            // Calculate rho for each phase
            calcRho<LT>(&f(0,0,nodeNo), rho0Node);  // LBmacroscopic
            rho(0, nodeNo) = rho0Node; // save to global field
            calcRho<LT>(&f(1,0,nodeNo), rho1Node);  // LBmacroscopic
            rho(1, nodeNo) = rho1Node; // save to global field

            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }  // End for all bulk nodes

        // Need to set the density when when calculating the densities
        // Predefine two regions:
        //   1) Pressure is constant
        //   2) Constant sink
        // Strategy - When mass is removed, then remove both oil and water
        //          - If mass is added at the top then then only add oil.

        for (auto nodeNo: constDensNodes) {
            lbBase_t q0, q1;
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);
            setConstDensity(q0, q1, rho0Node, rho1Node, 1.0);
            Q(0, nodeNo) = q0;
            Q(1, nodeNo) = q1;
            rho0Node += 0.5*q0;
            rho1Node += 0.5*q1;
            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }

        for (auto nodeNo: sourceNodes) {
            lbBase_t q0, q1;
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);
            setConstSource(q0, q1, rho0Node, rho1Node, -1e-4);
            Q(0, nodeNo) = q0;
            Q(1, nodeNo) = q1;
            rho0Node += 0.5*q0;
            rho1Node += 0.5*q1;
            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }


        for (auto nodeNo: solidBnd) { // Change getNumNodes to nNodes ?
            const lbBase_t rho0Node = rho(0, nodeNo);
            const lbBase_t rho1Node = rho(1, nodeNo);
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }


        //  MPI: COMMUNCATE SCALAR 'cgField'
        mpiBoundary.communciateScalarField(cgField);


        for (auto nodeNo: bulkNodes) { // Change nElements to nNodes?

            // UPDATE MACROSCOPIC VARIABLES
            lbBase_t rhoNode, rho0Node, rho1Node;
            //lbBase_t forceNode[LT::nD] = {0.0, 0.0};
            lbBase_t forceNode[LT::nD]; // = {0.0, 0.0, 0.0};
            lbBase_t velNode[LT::nD];
            lbBase_t fTot[LT::nQ];

            // Set the local total lb distribution
            for (int q = 0; q < LT::nQ; ++q)
                fTot[q] = f(0, q, nodeNo) + f(1, q, nodeNo);

            // Set densities
            rho0Node = rho(0, nodeNo);
            rho1Node = rho(1, nodeNo);
            // Total density
            rhoNode = rho0Node + rho1Node;
            // Total velocity and force
            forceNode[0] = 0.0;//(rho0Node - rho1Node) / rhoNode * bodyForce(0, 0, 0);
            forceNode[1] = rho0Node/rhoNode * bodyForce(0, 1, 0);
            forceNode[2] = 0.0;


            calcVel<LT>(fTot, rhoNode, velNode, forceNode);  // LBmacroscopic
            vel(0, 0, nodeNo) = velNode[0];
            vel(0, 1, nodeNo) = velNode[1];


            // Correct mass density for mass source
            lbBase_t q0Node = Q(0, nodeNo);
            lbBase_t q1Node = Q(1, nodeNo);
            rho0Node += 0.5*q0Node;
            rho1Node += 0.5*q1Node;
            rho(0, nodeNo) = rho0Node;
            rho(1, nodeNo) = rho1Node;


            // CALCULATE BGK COLLISION TERM
            // Mean collision time /rho_tot/\nu_tot = \sum_s \rho_s/\nu_s
            lbBase_t tau;
            tau = LT::c2Inv * rhoNode / (rho0Node*nu0Inv + rho1Node*nu1Inv) + 0.5;

            lbBase_t omegaBGK[LT::nQ]; // Holds the collision values
            lbBase_t uu, cu[LT::nQ];  // pre-calculated values
            uu = LT::dot(velNode, velNode);  // Square of the velocity
            LT::cDotAll(velNode, cu);  // velocity dotted with lattice vectors
            calcOmegaBGK<LT>(fTot, tau, rhoNode, uu, cu, omegaBGK);  // LBcollision


            // CALCULATE FORCE CORRECTION TERM
            lbBase_t deltaOmegaF[LT::nQ];
            lbBase_t  uF, cF[LT::nQ];
            uF = LT::dot(velNode, forceNode);
            LT::cDotAll(forceNode, cF);

            calcDeltaOmegaF<LT>(tau, cu, uF, cF, deltaOmegaF);  // LBcollision

            // CALCULATE MASS SOURCE CORRECTION TERM
            lbBase_t deltaOmegaQ0[LT::nQ], deltaOmegaQ1[LT::nQ];
            calcDeltaOmegaQ<LT>(tau, cu, uu, q0Node, deltaOmegaQ0);
            calcDeltaOmegaQ<LT>(tau, cu, uu, q1Node, deltaOmegaQ1);


            // -- calculate the normelaized color gradient
            lbBase_t CGNorm, CG2;
            lbBase_t colorGradNode[LT::nD];

            lbBase_t scalarTmp[LT::nQ];
            for (int q = 0; q < LT::nQ; ++q) {
                int neigNode = grid.neighbor(q, nodeNo);
                scalarTmp[q] = cgField(0, neigNode);
            }
            LT::grad(scalarTmp, colorGradNode);

            CG2 = LT::dot(colorGradNode, colorGradNode);
            CGNorm = sqrt(CG2);

            // Normalization of colorGrad
            for (int d = 0; d < LT::nD; ++d)
                colorGradNode[d] /= CGNorm + (CGNorm < lbBaseEps);

            // CALCULATE SURFACE TENSION PERTURBATION
            lbBase_t deltaOmegaST[LT::nQ];
            lbBase_t bwCos[LT::nQ];
            lbBase_t cCGNorm[LT::nQ];
            LT::cDotAll(colorGradNode, cCGNorm);
            lbBase_t AF0_5 = 2.25 * CGNorm * sigma / tau;
            lbBase_t rhoFacBeta = beta * rho0Node * rho1Node / rhoNode;

            for (int q = 0; q < LT::nQNonZero_; ++q) {
                deltaOmegaST[q] = AF0_5 * (LT::w[q] * cCGNorm[q]*cCGNorm[q] - LT::B[q] );
                bwCos[q] = rhoFacBeta * LT::w[q] * cCGNorm[q] /  LT::cNorm[q];
            }
            deltaOmegaST[LT::nQNonZero_] = -AF0_5 * LT::B[LT::nQNonZero_];
            bwCos[LT::nQNonZero_] = 0.0; // This should be zero by default


            // COLLISION AND PROPAGATION
            lbBase_t c0, c1;
            c0 = (rho0Node/rhoNode);  // Concentration of phase 0
            c1 = (rho1Node/rhoNode);  // Concentration of phase 1
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = c0 * (fTot[q] + omegaBGK[q] + deltaOmegaF[q] +  deltaOmegaST[q]) +  bwCos[q] + deltaOmegaQ0[q];
                fTmp(1, q,  grid.neighbor(q, nodeNo)) = c1 * (fTot[q] + omegaBGK[q] + deltaOmegaF[q] +  deltaOmegaST[q]) -  bwCos[q] + deltaOmegaQ1[q];
            }

        } // End nodes


        // PRINT

       if ( (i % static_cast<int>(input["iterations"]["write"])) == 0) {
           std::cout << "PLOT AT ITERATION : " << i << std::endl;
           std::string tmpName("/home/olau/Programs/Git/BADChIMP-cpp/output/rho_val_");
           tmpName += std::to_string(myRank) + "_" + std::to_string(i);
           tmpName += ".dat";
           std::ofstream ofs;
           ofs.open(tmpName);
            for (auto nodeNo: bulkNodes) {
                // std::cout << "(" << grid.pos(nodeNo, 0) << ", " << grid.pos(nodeNo, 1) << ") : (" << nodeNo << ", " << grid.getRank(nodeNo) << ") : ";
                // std::cout << rho(0, nodeNo) << " + " << rho(1, nodeNo) <<  " = " <<  rho(0, nodeNo) + rho(1, nodeNo) << std::endl;
                ofs << std::setprecision(23) << grid.pos(nodeNo, 0) << " " << grid.pos(nodeNo, 1) << " " << grid.pos(nodeNo, 2) << " " << rho(0, nodeNo) << " " << rho(1, nodeNo) << std::endl;
            }
            ofs.close();
        }

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // MPI Boundary
        // Hente verdier hente fra ghost
        // Sette i bulk
        mpiBoundary.communicateLbField(grid, f, 0);
        mpiBoundary.communicateLbField(grid, f, 1);


/*        if (myRank == 0) {
            for (int n = 1; n < grid.size(); ++n) {
                test[grid.pos(n, 1) + 1][grid.pos(n, 0) + 1] = f(0, 0, n);
            }
            for (int y = 0; y < 9; ++y) {
                for (int x = 0; x < 14; ++x)
                    std::cout << std::setw(3) << test[y][x];
                std::cout << std::endl;
            }
        } */



        // BOUNDARY CONDITIONS
        bbBnd.apply(0, f, grid);  // LBboundary
        bbBnd.apply(1, f, grid);


    } // End iterations
    // -----------------END MAIN LOOP------------------

    //mpi.end();

    MPI_Finalize();

    return 0;


}

