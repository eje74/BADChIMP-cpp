//=====================================================================================
//                      Q U E M A D A   R H E O L O G Y 
//
//                      A N N E L U S   G E O M E T R Y
//
//=====================================================================================

#include "../LBSOLVER.h"
#include "../IO.h"
#include "LBsixdof.h"
#include "LBrehology.h"
// SET THE LATTICE TYPE
#define LT D3Q19

int main()
{
    //---------------------------------------------------------------------------------
    //                                   Setup 
    //--------------------------------------------------------------------------------- 

    //                                   Setup 
    //--------------------------------------------------------------------------------- mpi
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    //                                   Setup 
    //--------------------------------------------------------------------------------- input output paths
    std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

    //                                   Setup 
    //--------------------------------------------------------------------------------- grid and geometry objects
    Input input(inputDir + "input.dat");
    LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
    Grid<LT> grid(vtklb);
    Nodes<LT> nodes(vtklb, grid);
    BndMpi<LT> mpiBoundary(vtklb, nodes, grid);
    //--------------------------------------------------------------------------------- Set bulk nodes
    std::vector<int> bulkNodes = findBulkNodes(nodes);


    //                                   Setup 
    //--------------------------------------------------------------------------------- read input-file
    // Number of iterations
    int nIterations = input["iterations"]["max"];
    // Write interval
    int nItrWrite = input["iterations"]["write"];
    // Relaxation time
    //lbBase_t tau = input["fluid"]["tau"];
    lbBase_t tauSym = 1.0;
    lbBase_t tauAnti = 1.0;

    //---------------------------------------------------------------------------------
    //                               Physical system 
    //--------------------------------------------------------------------------------- 

    //                               Physical system 
    //--------------------------------------------------------------------------------- Rheology object
    QuemadaLawRheology<LT> quemada(inputDir + "test.dat", 0.002361508);
    //                               Physical system 
    //--------------------------------------------------------------------------------- indicator field
    ScalarField phi(1, grid.size());
    vtklb.toAttribute("phi");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) 
    {
        const lbBase_t val = vtklb.getScalarAttribute<lbBase_t>();
        phi(0, nodeNo) = val;
    }


    //--------------------------------------------------------------------------------- position relative to origo
    VectorField<LT> pos(1, grid.size());
    vtklb.toAttribute("x0");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) 
    {
        const lbBase_t val = vtklb.getScalarAttribute<lbBase_t>();
        pos(0, 0, nodeNo) = val;
        pos(0, 2, nodeNo) = 0;
    }
    vtklb.toAttribute("y0");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) 
    {
        const lbBase_t val = vtklb.getScalarAttribute<lbBase_t>();
        pos(0, 1, nodeNo) = val;
    }

    //--------------------------------------------------------------------------------- boundary indicator
    ScalarField boundaryIndicator(1, grid.size());
    vtklb.toAttribute("boundary_indicator");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) 
    {
        const int val = vtklb.getScalarAttribute<int>();
        boundaryIndicator(0, nodeNo) = val;
    }

    //--------------------------------------------------------------------------------- boundary velocity
    VectorField<LT> boundaryVelocity(1, grid.size());
    vtklb.toAttribute("vel_bnd_x");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) 
    {
        const lbBase_t val = vtklb.getScalarAttribute<lbBase_t>();
        boundaryVelocity(0, 0, nodeNo) = val;
    }

    vtklb.toAttribute("vel_bnd_y");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) 
    {
        const lbBase_t val = vtklb.getScalarAttribute<lbBase_t>();
        boundaryVelocity(0, 1, nodeNo) = val;
        boundaryVelocity(0, 2, nodeNo) = 0.0;
    }




    //--------------------------------------------------------------------------------- pipe centre movement
    std::valarray<lbBase_t> pipePos(0.0, LT::nD);
    std::valarray<lbBase_t> pipeVel(0.0, LT::nD);
    std::valarray<lbBase_t> pipeAks(0.0, LT::nD);
    //--------------------------------------------------------------------------------- pipe rotation
    lbBase_t omega = 0;
    lbBase_t alpha = 0.0;
    lbBase_t omegaInner = 0.01;
    lbBase_t alphaInner = 0;
    lbBase_t omegaOuter = 0;
    lbBase_t rInner = 32;
    lbBase_t rOuter = 54;

    //--------------------------------------------------------------------------------- viscosity
    ScalarField viscosity(1, grid.size());
    //--------------------------------------------------------------------------------- rho
    ScalarField rho(1, grid.size());
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        rho(0, n) = 1.0;
    } 
    //--------------------------------------------------------------------------------- velocity
    VectorField<LT> vel(1, grid.size());
    VectorField<LT> velAna(1, grid.size());
    // Initiate velocity
    for (auto nodeNo: bulkNodes) {
        for (int d=0; d < LT::nD; ++d) 
        {
            vel(0, d, nodeNo) = 0.0;
            velAna(0, d, nodeNo) = 0.0;
        }
    }
    //--------------------------------------------------------------------------------- force
    VectorField<LT> force(1, grid.size());

    //                               Physical system 
    //--------------------------------------------------------------------------------- boundary conditions
    InterpolatedBounceBackBoundary<LT> ippBndOuter(findFluidBndNodes(nodes, boundaryIndicator, 1.0), nodes, grid);
    InterpolatedBounceBackBoundary<LT> ippBndInner(findFluidBndNodes(nodes, boundaryIndicator, 2.0), nodes, grid);
    // InterpolatedBounceBackBoundary<LT> ippBnd(findFluidBndNodes(nodes), nodes, grid);
    VectorField<LT> solidFluidForce(1, grid.size());

    //--------------------------------------------------------------------------------- lb fields
    LbField<LT> f(1, grid.size()); 
    LbField<LT> fTmp(1, grid.size());
    // initiate lb distributions
    //--------------------------------------------------------------------------------- lb fields - conservation of mass
    lbBase_t rhoTotInit = 0;
    lbBase_t rhoTot = 0;
    for (auto nodeNo: bulkNodes) {
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = LT::w[q]*rho(0, nodeNo);
        }
        rhoTotInit += rho(0,nodeNo);
    }
    rhoTot = rhoTotInit;

    //                               Output  
    ScalarField delta(1, grid.size());
    //--------------------------------------------------------------------------------- vtk output
    Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs);
    output.add_file("lb_run_annulus_immersed");
    //output.add_file("lb_run");
    output.add_scalar_variables({"phi", "boundary_indicator", "rho", "viscosity", "delta"}, {phi, boundaryIndicator,rho, viscosity, delta});
    output.add_vector_variables({"vel", "forceBoundary", "solidforce", "pos", "velAna"}, {vel, force, solidFluidForce, pos, velAna});
    

    //==================================================================================
    //                                  M A I N   L O O P  
    //==================================================================================
    std::ofstream writeDeltaP;
    std::ofstream writeSurfaceForce;
    std::ofstream writeTorque;
    if (myRank == 0 ) {
        writeDeltaP.open(outputDir + "deltaP.dat");
        writeSurfaceForce.open(outputDir + "surfaceforce.dat");
        writeTorque.open(outputDir + "torque.dat");
    }

    std::valarray<lbBase_t> surfaceForce(LT::nD);

    // lbBase_t u_mean = 0;
    const std::time_t beginLoop = std::time(NULL);
    for (int i = 0; i <= nIterations; i++) {
        if ( i < 1000) {
            const lbBase_t my_pi = 3.14159265358979323;
            omegaInner = 0.1*0.5*(1 - std::cos(i*my_pi/(1.0*1000)));
            alphaInner = (0.1*0.5*my_pi/1000.0)*std::sin(i*my_pi/(1.0*1000));
        } else {
            omegaInner = 0.1;
            alphaInner = 0.0;
        }
        //omegaInner *= 0.001;
        //                              main loop 
        //----------------------------------------------------------------------------- flux correction force - initiation
        /* 
        lbBase_t uxLocal = 0.0;
        for (auto nodeNo: bulkNodes) {
            auto cf = LT::qSumC(f(0, nodeNo));
            uxLocal += cf[2];
        }
        lbBase_t uxGlobal;
        MPI_Allreduce(&uxLocal, &uxGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        lbBase_t numNodesLocal = bulkNodes.size();
        lbBase_t numNodesGlobal;
        MPI_Allreduce(&numNodesLocal, &numNodesGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  

        //----------------------------------------------------------------------------- set flux
        const int Ninit = 1000;
        const lbBase_t u_low = 0.01;
        
        if (i < Ninit) {
            const lbBase_t my_pi = 3.14159265358979323;
            u_mean = 0.5*u_low*(1 - std::cos(i*my_pi/(1.0*Ninit)));
        }

        //----------------------------------------------------------------------------- set flux correction force
        lbBase_t forceX = 2*(u_mean*numNodesGlobal - uxGlobal)/numNodesGlobal;     
        if (myRank == 0) {
            writeDeltaP << u_mean << " " << forceX << std::endl;
        }
        */
        //                               main loop 
        //----------------------------------------------------------------------------- Begin node loop
        lbBase_t rhoTotInitGlobal;
        MPI_Allreduce(&rhoTotInit, &rhoTotInitGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        lbBase_t rhoTotGlobal;
        MPI_Allreduce(&rhoTot, &rhoTotGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        lbBase_t numNodesLocal = bulkNodes.size();
        lbBase_t numNodesGlobal;
        MPI_Allreduce(&numNodesLocal, &numNodesGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
        const lbBase_t dRhoNode = (rhoTotInitGlobal - rhoTotGlobal)/numNodesGlobal;
        rhoTot = 0;

        pipeVel[1] = 0.005*std::sin(3.14159*0.0003*i);
        pipeAks[1] = 0.005*3.14159*0.0003*std::cos(3.14159*0.0003*i);
        //pipeVel[1] = 0;
        pipePos = pipePos + pipeVel + 0.5*pipeAks;
        //pipePos[1] = 0;

        lbBase_t torque = 0;
        lbBase_t torqueAna = 0;
        lbBase_t forceX = 0;
        lbBase_t forceY = 0;


        for (auto nodeNo: bulkNodes) {

            

            // std::valarray<lbBase_t> fNode = f(0, nodeNo);
            // for (int q=0; q<LT::nD; ++q)
            //     fNode[q] += LT::w[q]*dRhoNode;
            //for (int q = )            
            const std::valarray<lbBase_t> fNode =  fRegularized<LT>(f(0, nodeNo), dRhoNode);

            // const std::valarray<lbBase_t> fluxCorrectionForce = {0, 0, forceX};
            /* std::valarray<lbBase_t> fluxCorrectionForce = {0, 0, 0};

            std::valarray<lbBase_t> forceNode = fluxCorrectionForce;
            force.set(0, nodeNo) = forceNode; 
            */
            //------------------------------------------------------------------------- rho 
            const lbBase_t rhoNode = calcRho<LT>(fNode);
            rhoTot += rhoNode;
            rho(0, nodeNo) = rhoNode;
            // pressure(0, nodeNo) = LT::c2*rhoNode + deltaP*laplacePressure(0, nodeNo);

            //------------------------------------------------------------------------- force
            const auto nodePos = pos(0, nodeNo);
            const lbBase_t dx = nodePos[0] - pipePos[0];
            const lbBase_t dy = nodePos[1] - pipePos[1];
            const lbBase_t rNormTrue = std::sqrt( dx*dx +  dy*dy );
            const lbBase_t rNorm = std::max(rNormTrue, 0.5); 
            const lbBase_t dist = rNormTrue - rInner; 
            boundaryIndicator(0, nodeNo) = dist;
            std::valarray<lbBase_t> tangent(LT::nD);
            tangent[0] = -dy/rNorm;
            tangent[1] =  dx/rNorm;
            tangent[2] = 0;
            const std::valarray<lbBase_t> velInner = omegaInner*tangent + pipeVel;
            const lbBase_t omega1 = omegaInner/rInner;
            const lbBase_t alpha1 = alphaInner/rInner;

            //------------------------------------------------------------------------- force rotation
            std::valarray<lbBase_t> forceRot(0.0, 3);

            forceRot[0] = dx;
            forceRot[1] = dy;
            forceRot *= -rhoNode*omega1*omega1; //*heavisideStepReg<LT>(dist);

            forceRot[0] += -2*rhoNode*omega1*(vel(0,1,nodeNo) - (rNormTrue*omega1*tangent[1] + pipeVel[1]));
            forceRot[1] +=  2*rhoNode*omega1*(vel(0,0,nodeNo) - (rNormTrue*omega1*tangent[0] + pipeVel[1]));

            forceRot[0] += -rhoNode*alpha1*dy;
            forceRot[1] +=  rhoNode*alpha1*dx;
            //------------------------------------------------------------------------- force translation
            std::valarray<lbBase_t> forceTra(0.0, 3);
            forceTra = rhoNode*pipeAks;

            forceRot += forceTra;

            forceRot *= heavisideStepReg<LT>(dist);

            //forceRot += 0*2*0.01*(rhoNode*rNormTrue*omega1*tangent - Mi - 0.5*forceRot)*heavisideStepReg<LT>(dist+2);

            const auto Mi = LT::qSumC(fNode);
            const std::valarray<lbBase_t> forceBoundary = 2*(rhoNode*velInner - Mi - 0.5*forceRot)*dirceDeltaReg(dist);
//            const auto forceBoundary = immersedBoundaryForce<LT>(dist, fNode, rhoNode, omegaInner, velInner);
            
            torque += (tangent[0]*forceBoundary[0] + tangent[1]*forceBoundary[1])*rInner;
            forceX += forceBoundary[0];
            forceY += forceBoundary[1];

            //const std::valarray<lbBase_t> velInternal = omegaInner*velRot*(rInner - 5)/(rInner);
            //const auto forceInner = immersedBoundaryForce<LT>(dist+5, fNode, rhoNode, omegaInner, velInternal);

            const std::valarray<lbBase_t> forceNode = forceBoundary + forceRot;// + forceInner;
            force.set(0, nodeNo) = forceBoundary;

            //------------------------------------------------------------------------- velocity 
            auto velNode = calcVel<LT>(fNode, rhoNode, forceNode);
            const lbBase_t speed = std::sqrt(velNode[0]*velNode[0] + velNode[1]*velNode[1] + velNode[2]*velNode[2]);
            if (speed > 0.3)
                velNode = 0.3*velNode/speed;

            vel.set(0, nodeNo) = velNode;
            //------------------------------------------------------------------------- Calculation of the analytical velocity
            lbBase_t speedR;
            const lbBase_t eta = rInner/rOuter;
            const lbBase_t my = 0; // omegaOuter/omegaInner
            const lbBase_t A = omega1*(my - eta*eta)/(1 - eta*eta);
            const lbBase_t B = omega1*rInner*rInner*(1 - my)/(1 - eta*eta);
            if (rNormTrue >= rInner) {
                speedR = A*rNormTrue + B/rNormTrue;
                torqueAna += -(1.0/3.0)*(tauSym - 0.5)*rInner*dirceDeltaReg(rNormTrue - rInner)*(A - B/(rInner*rInner));
            } else {
                speedR = omega1*rNormTrue;
                torqueAna += -(1.0/3.0)*(tauSym - 0.5)*rInner*dirceDeltaReg(rNormTrue - rInner)*(A - B/(rInner*rInner));
//                torqueAna += 0*(1.0/3.0)*(tauSym - 0.5)*rInner*dirceDeltaReg(rNormTrue - rInner)*omega1;
            }
            std::valarray<lbBase_t> velAnaNode(LT::nD);
            velAnaNode = speedR*tangent;
            velAna.set(0, nodeNo) = velAnaNode;
            
            //torqueAna += -(1.0/3.0)*(tauSym - 0.5)*rInner*dirceDeltaReg(rNormTrue - rInner)*(A - B/(rInner*rInner));

            //------------------------------------------------------------------------- force + vel + rot 
            // const std::valarray<lbBase_t> posNode = pos(0, nodeNo);
            // const auto velNode = calcVelRot<LT>(fNode, rhoNode, posNode, fluxCorrectionForce, omega, alpha);
            // vel.set(0, nodeNo) = velNode;
            // const std::valarray<lbBase_t> forceNode = pseudoForce(posNode, velNode, omega, alpha) + fluxCorrectionForce;
            // force.set(0,nodeNo) = forceNode;
            //                           main loop 
            //------------------------------------------------------------------------- BGK-collision term
            const lbBase_t u2 = LT::dot(velNode, velNode);
            const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
            //------------------------------------------------------------------------- tau
            // auto omegaBGK = carreau.omegaBGK(fNode, rhoNode, velNode, u2, cu, forceNode, 0);
            // tau = powerLaw.tau(fNode, rhoNode, velNode, u2, cu, forceNode);
            tauSym = quemada.tau(fNode, rhoNode, velNode, u2, cu, forceNode);
//            tau = 0.502788344475428;
//            tau = 5055766889508577;    
//            tauSym = 0.5055766889508577;
            tauSym += (1- tauSym)*heavisideStepReg<LT>(dist + 4.5, 2.5);
            tauAnti = 1.0;
            // viscosity(0, nodeNo) = tau;
            viscosity(0, nodeNo) = tauSym;
            // auto omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);
            const auto omegaBGK = calcOmegaBGKTRT<LT>(fNode, tauSym, tauAnti, rhoNode, u2, cu);
            //------------------------------------------------------------------------- Guo-force correction
            const lbBase_t uF = LT::dot(velNode, forceNode);
            const std::valarray<lbBase_t> cF = LT::cDotAll(forceNode);
            //tau = carreau.tau();
            // const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);
            const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaFTRT<LT>(tauSym, tauAnti, cu, uF, cF);

            //------------------------------------------------------------------------- Collision and propagation
            fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);
        } //--------------------------------------------------------------------------- End node loop
     

        // std::cout << rhoTot << " " << rhoTotInit << " " << dRhoNode << std::endl;
        //                               main loop 
        //----------------------------------------------------------------------------- Swap fTmp and f
        f.swapData(fTmp);  
        //                            bondary conditions
        //----------------------------------------------------------------------------- mpi
        mpiBoundary.communicateLbField(0, f, grid);
        //----------------------------------------------------------------------------- inlet, outlet and solid
        // bounceBackBnd.apply(f, grid);
        //ippBnd.apply(phi, 0, f, grid);
        // surfaceForce = ippBnd.apply(phi, 0, f, grid, solidFluidForce);
        ippBndOuter.applyTest(phi, rho, omegaOuter, boundaryVelocity, 0, f, nodes, grid, solidFluidForce);
        surfaceForce = ippBndInner.applyTest(phi, rho, omegaInner, boundaryVelocity, 0, f, nodes, grid, solidFluidForce);
        lbBase_t FzLocal;
        lbBase_t FzGlobal;
        MPI_Allreduce(&FzLocal, &FzGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        lbBase_t torqueGlobal;
        MPI_Allreduce(&torque, &torqueGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        lbBase_t torqueAnaGlobal;
        MPI_Allreduce(&torqueAna, &torqueAnaGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        lbBase_t forceXGlobal, forceYGlobal;
        MPI_Allreduce(&forceX, &forceXGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&forceY, &forceYGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


        if (myRank == 0) {
            writeSurfaceForce << FzGlobal << std::endl;
            writeTorque << torqueGlobal << " " << forceXGlobal << " " << forceYGlobal << std::endl;
        }
        // abbBnd.apply(0, f, vel, pind, grid);
        /*wetBoundary.applyBoundaryCondition(0, f, force, grid);
        lbBase_t deltaMass  = 0.0;
        for (const auto & nodeNo: bulkNodes) {
            deltaMass += 1 - LT::qSum(f(0, nodeNo));
        }
        wetBoundary.addMassToWallBoundary(0, f, force, deltaMass);
        */

        //                               main loop 
        //----------------------------------------------------------------------------- write to file
        if ( ((i % nItrWrite) == 0)  ) {
            output.write(i);
            
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
                std::cout << "TORQUE = " << torqueGlobal  << "  " << torqueAnaGlobal << std::endl;
            }
        }

    } //------------------------------------------------------------------------------- End iterations
    if(myRank == 0) {
        std::cout << " RUNTIME = " << std::time(NULL) - beginLoop << std::endl;
    }
    writeDeltaP.close();
    writeSurfaceForce.close();
    writeTorque.close();

    MPI_Finalize();
    return 0;
}
