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
// TO RUN PROGRAM: type "mpirun -np <#procs> bdchmp" in command
// line in main directory
//
// //////////////////////////////////////////////

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "../lbsolver/LBbndmpi.h"
#include "../lbsolver/LBboundary.h"
#include "../lbsolver/LBcollision.h"
#include "../lbsolver/LBcollision2phase.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBgeometry.h"
#include "../lbsolver/LBgrid.h"
#include "../lbsolver/LBhalfwaybb.h"
#include "../lbsolver/LBinitiatefield.h"
#include "../lbsolver/LBmacroscopic.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBsnippets.h"
#include "../lbsolver/LButilities.h"

//#include "../io/Input.h"
//#include "../io/Output.h"
#include "../IO.h"

#include "../lbsolver/LBvtk.h"

#include<algorithm> // std::max

// SET THE LATTICE TYPE
//#define LT D2Q9
//#define VTK_CELL VTK::pixel
#define LT D3Q19
#define VTK_CELL VTK::voxel

int main()
{
    // SETUP MPI
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // ********************************
    // SETUP THE INPUT AND OUTPUT PATHS
    // ********************************
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

    // SETUP BOUNCE BACK BOUNDARY (fluid boundary)
    std::cout << "bbBnd" << std::endl;
    HalfWayBounceBack<LT> bbBnd(findBulkNodes(nodes), nodes, grid); // = makeFluidBoundary<HalfWayBounceBack>(nodes, grid);

    // SETUP SOLID BOUNDARY
    std::vector<int> solidBnd = findSolidBndNodes(nodes);

    // SETUP BULK NODES
    std::vector<int> bulkNodes = findBulkNodes(nodes);

    // ****************
    // SETUP MASS SOURCES
    // ****************
    // Source markers:
    std::vector<int> sourceMarker(grid.size());
    //   read from file
    vtklb.toAttribute("source");
    std::vector<int> constDensNodes, sourceNodes;
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        int val = vtklb.getScalar<int>();
        sourceMarker[nodeNo] = val;
        if (val == 1) {// sink
    	  sourceNodes.push_back(nodeNo);
        } else if (val == 2) {// Const pressure
	      constDensNodes.push_back(nodeNo);
    	}
    }


    // Vector source
    VectorField<LT> bodyForce(1, 1);
    // bodyForce.set(0, 0) = inputAsValarray<lbBase_t>(input["fluid"]["bodyforce"]);
    bodyForce.set(0, 0) = input["fluid"]["bodyforce"];

    //int nIterations = static_cast<int>( input["iterations"]["max"]);
    int nIterations = input["iterations"]["max"];

    lbBase_t tau0 = input["fluid"]["tau"][0];
    lbBase_t tau1 = input["fluid"]["tau"][1];
    lbBase_t sigma = input["fluid"]["sigma"];
    lbBase_t beta = input["fluid"]["beta"];

    // SET DERIVED VARIABLES
    lbBase_t nu0Inv = 1.0 / (LT::c2 * (tau0 - 0.5));
    lbBase_t nu1Inv = 1.0 / (LT::c2 * (tau1 - 0.5));

    // Scalar source
    ScalarField Q(2, grid.size());
    for (int n = 0; n < Q.size(); ++n) {
        Q(0, n) = 0.0;
        Q(1, n) = 0.0;
    }
    // SETUP LB FIELDS
    LbField<LT> f(2, grid.size());  // LBfield
    LbField<LT> fTmp(2, grid.size());  // LBfield

    // SETUP MACROSCOPIC FIELDS
    ScalarField rho(2, grid.size()); // LBfield
    VectorField<LT> vel(1, grid.size()); // LBfield
    ScalarField cgField(1, grid.size()); // LBfield


    // FILL MACROSCOPIC FIELDS
    //   Fluid densities and velocity
    std::srand(8549388);

     // ****************
    // SETUP MASS DENSITIES
    // ****************
    // Density markers:
    // read from file
    vtklb.toAttribute("rho0");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        float val = vtklb.getScalar<float>();
        rho(0, nodeNo) = val;
	for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }

    vtklb.toAttribute("rho1");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        float val = vtklb.getScalar<float>();
        rho(1, nodeNo) = val;
	for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }
    vtklb.toAttribute("wettability");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
      float val = vtklb.getScalar<float>();
      if(nodes.isSolidBoundary(nodeNo)){
	rho(0, nodeNo) =val;
	rho(1, nodeNo) =1-val;
      }
    }

    /*
    for (auto nodeNo: bulkNodes) {
        rho(0, nodeNo) = 1.0; //(1.0 * std::rand()) / (RAND_MAX * 1.0);
        rho(1, nodeNo) = 1 - rho(0, nodeNo);

        for (int d=0; d < LT::nD; ++d)
            vel(0, d, nodeNo) = 0.0;
    }
    */
    /*
    //   Solid boundary (Wettability)
    for (auto nodeNo: solidBnd) {
        rho(0, nodeNo) = 1.0;
        rho(1, nodeNo) = 1.0 - rho(0, nodeNo);
    }
    */

    int numNodes = bulkNodes.size();
    int numNodesGlobal;
    MPI_Allreduce(&numNodes, &numNodesGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);



    // INITIATE LB FIELDS
    // -- phase 0
    initiateLbField(0, 0, 0, bulkNodes, rho, vel, f);  // LBinitiatefield
    // -- phase 1
    initiateLbField(1, 1, 0, bulkNodes, rho, vel, f);  // LBinitiatefield

    //std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));
    //std::string dirNum = input["out"]["directoryNum"]; 
    //std::string outDir2 = outputDir+"out"+dirNum;
    outputDir += input["outdir"];

    // **********
    // OUTPUT VTK
    // **********
    Output<LT> output(grid.pos(bulkNodes), outputDir, myRank, nProcs);
    output.add_file("fluid");
    output.add_variables({{"rho", &rho}, {"vel", &vel}}, bulkNodes);
    // output.add_variables({"rho0", "rho1"}, rho, bulkNodes);
    // output.add_variables({"vel"}, vel, bulkNodes);
    //output.write(0);

    // std::vector<int> geo(grid.size(), -1);
    // for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
    //     geo[nodeNo] = nodes.isSolid(nodeNo) ? 1 : 0;        
    // }
    Output<LT,int> geoout(grid.pos(), outputDir, myRank, nProcs, "geo", nodes.geo(grid, vtklb), 1);
    geoout.write();
    //outputGeometry("geo", outputDir, myRank, nProcs, nodes, grid, vtklb);
    // JLV

    // -----------------MAIN LOOP------------------
    /* Comments to main loop:
     * Calculation of cu is kept outside of calcOmega and calcDeltaOmega
     *  to avoid recalculation  of the dot-product.
     * It is therefore natural to give all scalar products as inputs to the
     *  calcOmega and calcDeltaOmega.
     * Note that we have a oneIteration function in 'LBiteration.h'. This seemed
     *   to affect the speed an gave a slow down of 10-15 %.
     */

    for (int i = 0; i <= nIterations; i++) {

        for (auto nodeNo : bulkNodes) {
            // UPDATE MACROSCOPIC DENSITIES
            // Calculate rho for each phase
            lbBase_t rho0Node = rho(0, nodeNo) = calcRho<LT>(f(0,nodeNo));  // LBmacroscopic
            lbBase_t rho1Node = rho(1, nodeNo) = calcRho<LT>(f(1,nodeNo));  // LBmacroscopic

            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }  // End for all bulk nodes

        // Need to set the density when when calculating the densities
        // Predefine two regions:
        //   1) Pressure is constant
        //   2) Constant sink
        // Strategy - When mass is removed, then remove both oil and water
        //          - If mass is added at the top then then only add oil.

	/*
	for (auto nodeNo: constDensNodes) {
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);

            setConstDensity(Q(0, nodeNo), Q(1, nodeNo), rho0Node, rho1Node, 1.0);
            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }

        for (auto nodeNo: sourceNodes) {
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);
	    lbBase_t fluxNode = -1e-3;
            setConstSource(Q(0, nodeNo), Q(1, nodeNo), rho0Node, rho1Node, fluxNode);
	    if((rho0Node+rho1Node)<0.8){
	      Q(0, nodeNo)=Q(1, nodeNo)=0.0;
	    }

            // Calculate color gradient kernel
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }
	*/


        for (auto nodeNo: solidBnd) {
            const lbBase_t rho0Node = rho(0, nodeNo);
            const lbBase_t rho1Node = rho(1, nodeNo);
            cgField(0, nodeNo) = (rho0Node - rho1Node)/(rho0Node + rho1Node);
        }

        //  MPI: COMMUNCATE SCALAR 'cgField'
        mpiBoundary.communciateScalarField(cgField);


	//----------------------------------Flux---------------------------------------

	lbBase_t meanfcX=0.0;
	lbBase_t meanfcXGlobal;
	for (auto nodeNo : bulkNodes) {
	  std::valarray<lbBase_t> fTot = f(0, nodeNo) + f(1, nodeNo);
	  std::valarray<lbBase_t> sumfc = LT::qSumC(fTot);
	  meanfcX+=sumfc[0];
	}
	MPI_Allreduce(&meanfcX, &meanfcXGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	meanfcXGlobal/=numNodesGlobal;

	lbBase_t setMomX=input["fluid"]["momx"]; //1e-3;
	//if (myRank==0) {
	//  std::cout << "setMomX: " << setMomX << std::endl;
	//}
	lbBase_t fluxForceX=2*(setMomX-meanfcXGlobal);
	bodyForce(0,0,0)=fluxForceX;
	//----------------------------------end Flux---------------------------------------


        for (auto nodeNo: bulkNodes) {

            // Set the local total lb distribution
            // Use std:valarray insatead of auto as valarray uses expression templates that do not work well with auto
            std::valarray<lbBase_t> fTot = f(0, nodeNo) + f(1, nodeNo);

            // UPDATE MACROSCOPIC VARIABLES
            // -- densities
            lbBase_t rho0Node = rho(0, nodeNo);
            lbBase_t rho1Node = rho(1, nodeNo);
            // -- total density
            lbBase_t rhoNode = rho0Node + rho1Node;
            // -- force
            std::valarray<lbBase_t> forceNode = bodyForce(0, 0);//setForceGravity(rho0Node, rho1Node, bodyForce, 0);
            // -- velocity
            std::valarray<lbBase_t> velNode = calcVel<LT>(fTot, rhoNode, forceNode);  // LBmacroscopic

            vel.set(0, nodeNo) = velNode;

            // Correct mass density for mass source
            lbBase_t q0Node = Q(0, nodeNo);
            lbBase_t q1Node = Q(1, nodeNo);
            rho(0, nodeNo) = rho0Node += 0.5*q0Node;
            rho(1, nodeNo) = rho1Node += 0.5*q1Node;


            // CALCULATE BGK COLLISION TERM
            // Mean collision time /rho_tot/\nu_tot = \sum_s \rho_s/\nu_s
            lbBase_t tau = LT::c2Inv * rhoNode / (rho0Node*nu0Inv + rho1Node*nu1Inv) + 0.5;

            lbBase_t uu = LT::dot(velNode, velNode);  // Square of the velocity
            std::valarray<lbBase_t> cu = LT::cDotAll(velNode);  // velocity dotted with lattice vectors
            std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fTot, tau, rhoNode, uu, cu);  // LBcollision


            // CALCULATE FORCE CORRECTION TERM
            lbBase_t  uF = LT::dot(velNode, forceNode);
            std::valarray<lbBase_t>  cF = LT::cDotAll(forceNode);
            std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);  // LBcollision

            // CALCULATE MASS SOURCE CORRECTION TERM
            std::valarray<lbBase_t> deltaOmegaQ0 = calcDeltaOmegaQ<LT>(tau, cu, uu, q0Node);
            std::valarray<lbBase_t> deltaOmegaQ1 = calcDeltaOmegaQ<LT>(tau, cu, uu, q1Node);

            // CALCULATE SURFACE TENSION PERTURBATION
            // -- calculate the normalized color gradient
            std::valarray<lbBase_t> colorGradNode = grad(cgField, 0, nodeNo, grid);
            lbBase_t CGNorm = vecNorm<LT>(colorGradNode);
            colorGradNode *= 1.0/(CGNorm + (CGNorm < lbBaseEps));

            std::valarray<lbBase_t> cCGNorm = LT::cDotAll(colorGradNode);

            std::valarray<lbBase_t> deltaOmegaST = calcDeltaOmegaST<LT>(tau, sigma, CGNorm, cCGNorm);
            std::valarray<lbBase_t> deltaOmegaRC = calcDeltaOmegaRC<LT>(beta, rho0Node, rho1Node, rhoNode, cCGNorm);

            // COLLISION AND PROPAGATION
            lbBase_t c0, c1;
            c0 = (rho0Node/rhoNode);  // Concentration of phase 0
            c1 = (rho1Node/rhoNode);  // Concentration of phase 1
            for (int q = 0; q < LT::nQ; ++q) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = c0 * (fTot[q] + omegaBGK[q] + deltaOmegaF[q] +  deltaOmegaST[q]) +  deltaOmegaRC[q]  + deltaOmegaQ0[q];
                fTmp(1, q,  grid.neighbor(q, nodeNo)) = c1 * (fTot[q] + omegaBGK[q] + deltaOmegaF[q] +  deltaOmegaST[q]) -  deltaOmegaRC[q]  + deltaOmegaQ1[q];
            }

        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield

        // MPI Boundary
        // Hente verdier hente fra ghost
        // Sette i bulk


        mpiBoundary.communicateLbField(0, f, grid);
        mpiBoundary.communicateLbField(1, f, grid);


        // BOUNDARY CONDITIONS
        bbBnd.apply(0, f, grid);  // LBboundary
        bbBnd.apply(1, f, grid);




        // PRINT
        if ( (i % static_cast<int>(input["iterations"]["write"])) == 0)
        {



            if (myRank==0)
                std::cout << "PLOT AT ITERATION : " << i << std::endl;

/*            std::string tmpName(outDir);
            tmpName += std::to_string(myRank) + "_" + std::to_string(i);
            tmpName += ".dat";
            std::ofstream ofs;
            ofs.open(tmpName);
            if (!ofs) {
                std::cout << "Error: could not open file: " << tmpName << std::endl;
                return 1;
            }

            for (auto nodeNo: bulkNodes) {
                ofs << std::setprecision(23) << grid.pos(nodeNo, 0) << " " << grid.pos(nodeNo, 1) << " " << grid.pos(nodeNo, 2) << " " << rho(0, nodeNo) << " " << rho(1, nodeNo)
                    << " " << vel(0, 0, nodeNo) << " " << vel(0, 1, nodeNo) << " " << vel(0, 2, nodeNo) << " " << nodes.getType(nodeNo) << std::endl;
            }
            ofs.close(); */

            // JLV
            output.write(i);

	    if (myRank==0){
	      std::ofstream ofs;
	      //std::string tmpName(outDir2+"/force.dat");
	      std::string tmpName(outputDir+"/force.dat");
	      ofs.open(tmpName, std::ios::app);
	      if (!ofs) {
                std::cout << "Error: could not open file: " << tmpName << std::endl;
                return 1;
	      }
	      ofs << i << " " << std::setprecision(23) << bodyForce(0, 0, 0) << std::endl;

	      ofs.close();
	    }
	    /*
	    for (auto nodeNo: bulkNodes) {
	      if(rho(1,nodeNo)>1+1e-8 || rho(1,nodeNo)<0.0){
		std::cout<< rho(1,nodeNo) << " ";
	      }

	    }
	    */
            // JLV
        }

    } // End iterations
    // -----------------END MAIN LOOP------------------

    //mpi.end();

    MPI_Finalize();

    return 0;
}
