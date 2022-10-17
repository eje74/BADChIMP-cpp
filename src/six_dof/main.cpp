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


// HELPER CONSTANTS AND FUNCTIONS
#define PI 3.14159265358979323

lbBase_t gVal(const lbBase_t x)
{
    return 0.5*(1-std::cos(x));
}


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
    lbBase_t tauSymOrg = 0.5005639348377272;


    //---------------------------------------------------------------------------------
    //                               Physical system 
    //--------------------------------------------------------------------------------- 
    //--------------------------------------------------------------------------------- pipe rotation
    // const lbBase_t my_pi = 3.14159265358979323;
    lbBase_t dl = 0.05 *1e-2;
    lbBase_t dt = 1.2531885282826405e-05;
    //lbBase_t omegaInner = 0.01;
    //lbBase_t alphaInner = 0;
    lbBase_t omegaOuter = 0;
    const lbBase_t rInner = 0.5*0.127/dl;
    const lbBase_t rOuter = 0.5*0.216/dl;
    //--------------------------------------------------------------------------------- Constants for CM-movment
    // Latteral
    const lbBase_t cm_amp_r = 0.9*(rOuter - rInner);
    //const lbBase_t cm_amp_z = 0.5*cm_amp_r;
    const lbBase_t tmp_w = 2*0.005/(cm_amp_r*12.0 + 0.5*cm_amp_r*5.0);
    
    const lbBase_t cm_w_r = 7*tmp_w;
    const lbBase_t cm_w_x = 2*tmp_w;
    const lbBase_t cm_w_y = 5*tmp_w;
    // Axial
    const lbBase_t cm_amp_z = 0.5*cm_amp_r*3*tmp_w;
    const lbBase_t cm_w_z = 3*tmp_w;
    //--------------------------------------------------------------------------------- Constants for rotation
    const lbBase_t rot_amp = 0.03/rInner;
    const lbBase_t rot_w = tmp_w;
    //--------------------------------------------------------------------------------- Constants for forcing
    const lbBase_t dp_amp = 0.01*8*LT::c2*(tauSymOrg-0.5);
    const lbBase_t dp_w = 0.005*tmp_w;

    std::cout << cm_amp_r << " " << tmp_w << std::endl;


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

    /* const lbBase_t A1_r = 0.021/dl;//0.06/dl;
    const lbBase_t w1_r = 1.1*dt;
    const lbBase_t omega1_r = 2*my_pi*0.5*dt;
    const lbBase_t A2_r = 0.021/dl;
    const lbBase_t w2_r = 0.7*dt;
    const lbBase_t omega2_r = 2*my_pi*0.3*dt;
    */
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
    //InterpolatedBounceBackBoundary<LT> ippBndInner(findFluidBndNodes(nodes, boundaryIndicator, 2.0), nodes, grid);
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
    output.add_scalar_variables({"phi", "boundary_indicator", "rho", "viscosity"}, {phi, boundaryIndicator,rho, viscosity});
    output.add_vector_variables({"vel", "forceBoundary", "solidforce"}, {vel, force, solidFluidForce});
    

    // HELPER
    std::valarray<lbBase_t> oldPos(0.0, 3);
    std::valarray<lbBase_t> oldVel(0.0, 3);

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
        /*if ( i < 1000) {
            omegaInner = 0.05*0.5*(1 - std::cos(i*my_pi/(1.0*1000)));
            alphaInner = (0.05*0.5*my_pi/1000.0)*std::sin(i*my_pi/(1.0*1000));
        } else {
            omegaInner = 0.05;
            alphaInner = 0.0;
        } */
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

/*        lbBase_t r_1 = A1_r*std::sin(w1_r*i);
        lbBase_t r_dot_1 = w1_r*A1_r*std::cos(w1_r*i);
        lbBase_t r_dot_dot_1 = -w1_r*w1_r*A1_r*std::sin(w1_r*i);

        lbBase_t r_2 = A2_r*std::sin(w2_r*i);
        lbBase_t r_dot_2 = w2_r*A2_r*std::cos(w2_r*i);
        lbBase_t r_dot_dot_2 = -w2_r*w2_r*A2_r*std::sin(w2_r*i);


        const std::valarray<lbBase_t> rVec1{std::cos(omega1_r*i), std::sin(omega1_r*i)};
        const std::valarray<lbBase_t> rVec2{-std::sin(omega1_r*i), std::cos(omega1_r*i)};
        
        pipePos = r_1*rVec1;
        pipeVel = r_dot_1*rVec1 + omega1_r*r_1*rVec2;
        pipeAks = (r_dot_dot_1 - omega1_r*omega1_r*r_1)*rVec1 + 2*omega1_r*r_dot_1*rVec2;


        const std::valarray<lbBase_t> r2Vec1{std::cos(omega2_r*i), std::sin(omega2_r*i)};
        const std::valarray<lbBase_t> r2Vec2{-std::sin(omega2_r*i), std::cos(omega2_r*i)};
        
        pipePos += r_2*r2Vec1;
        pipeVel += r_dot_2*r2Vec1 + omega2_r*r_2*r2Vec2;
        pipeAks += (r_dot_dot_2 - omega2_r*omega2_r*r_2)*r2Vec1 + 2*omega2_r*r_dot_2*r2Vec2;
*/


        //------------------------------------------------------------------------------------------Pipe pos
        const lbBase_t gw = gVal(cm_w_r*i);
        const lbBase_t gx = gVal(cm_w_x*i);
        const lbBase_t gy = gVal(cm_w_y*i);
        const lbBase_t gz = gVal(cm_w_z*i);

        pipePos[0] = cm_amp_r*gw*gx;
        pipePos[1] = cm_amp_r*gw*gy; 
        pipePos[2] = 0;
        
        //------------------------------------------------------------------------------------------Pipe vel
        const lbBase_t dgw = cm_w_r*0.5*std::sin(cm_w_r*i);
        const lbBase_t dgx = cm_w_x*0.5*std::sin(cm_w_x*i);
        const lbBase_t dgy = cm_w_y*0.5*std::sin(cm_w_y*i);
        const lbBase_t dgz = cm_w_z*0.5*std::sin(cm_w_z*i);

        pipeVel[0] = cm_amp_r*(gw*dgx + dgw*gx);
        pipeVel[1] = cm_amp_r*(gw*dgy + dgw*gy);
        pipeVel[2] = cm_amp_z*gz;



        //------------------------------------------------------------------------------------------Pipe aks
        const lbBase_t ddgw = cm_w_r*cm_w_r*0.5*std::cos(cm_w_r*i);
        const lbBase_t ddgx = cm_w_x*cm_w_x*0.5*std::cos(cm_w_x*i);
        const lbBase_t ddgy = cm_w_y*cm_w_y*0.5*std::cos(cm_w_y*i);
        //const lbBase_t ddgz = cm_w_z*cm_w_z*0.5*std::cos(cm_w_z*i);

        pipeAks[0] = cm_amp_r*(gw*ddgx + 2*dgw*dgx + ddgw*gx);
        pipeAks[1] = cm_amp_r*(gw*ddgy + 2*dgw*dgy + ddgw*gy);
        pipeAks[2] = cm_amp_z*dgz;


        /*if (myRank == 0) {
            std::valarray<lbBase_t> dPos = pipePos - oldPos;
            oldPos = pipePos;
            std::valarray<lbBase_t> dVel = pipeVel - oldVel;
            oldVel = pipeVel;

            std::cout << std::endl;
            std::cout << pipeVel[0] << " " << pipeVel[1] << " " << pipeVel[2] << std::endl;    
            std::cout << dPos[0] << " " << dPos[1] << " " << dPos[2] << std::endl;    
            std::cout << std::endl;
            std::cout << pipeAks[0] << " " << pipeAks[1] << " " << pipeAks[2] << std::endl;    
            std::cout << dVel[0] << " " << dVel[1] << " " << dVel[2] << std::endl;    
            std::cout << std::endl;
            std::cout << std::endl;
        } */


        //------------------------------------------------------------------------------------------Disk angular speed
        const lbBase_t omega_z = rot_amp*gVal(rot_w*i);
        const lbBase_t omega_z_sq = omega_z*omega_z;
        const lbBase_t omegaInner = omega_z*rInner;

        //------------------------------------------------------------------------------------------Disk angular velocity
        const lbBase_t alpha_z = rot_amp*rot_w*0.5*std::sin(rot_w*i);

        //------------------------------------------------------------------------------------------Axial pressure
        const lbBase_t dp = dp_amp*gVal(dp_w*i);


        // pipeVel[1] = 0.005*std::sin(3.14159*0.0003*i);
        // pipeAks[1] = 0.005*3.14159*0.0003*std::cos(3.14159*0.0003*i);
        //pipeVel[1] = 0;
        // pipePos = pipePos + pipeVel + 0.5*pipeAks;
        //pipePos[1] = 0;

        lbBase_t torqueX = 0;
        lbBase_t torqueY = 0;
        lbBase_t torqueZ = 0;
        lbBase_t torqueAna = 0;
        lbBase_t forceX = 0;
        lbBase_t forceY = 0;
        lbBase_t forceZ = 0;
        lbBase_t velZ = 0;


        for (auto nodeNo: bulkNodes) {
            //----------------------------------------------------------------------------------------------------------- regularized 
            const std::valarray<lbBase_t> fNode =  fRegularized<LT>(f(0, nodeNo), dRhoNode);
            //----------------------------------------------------------------------------------------------------------- rho 
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
            const lbBase_t heavsideInner = heavisideStepReg<LT>(dist);
            boundaryIndicator(0, nodeNo) = dist;
            std::valarray<lbBase_t> tangent(LT::nD);
            tangent[0] = -dy/rNorm;
            tangent[1] =  dx/rNorm;
            tangent[2] = 0;
            const std::valarray<lbBase_t> velInner = omegaInner*tangent + pipeVel;
            //const lbBase_t omega1 = omegaInner/rInner;
            //const lbBase_t alpha1 = alphaInner/rInner;

            //------------------------------------------------------------------------- force rotation
            std::valarray<lbBase_t> forceRot(0.0, 3);

            forceRot[0] = dx;
            forceRot[1] = dy;
            forceRot *= -rhoNode*omega_z_sq; //  omega1*omega1; //*heavisideStepReg<LT>(dist);

            //forceRot[0] += -2*rhoNode*omega1*(vel(0,1,nodeNo) - (rNormTrue*omega1*tangent[1] + pipeVel[1]));
            //forceRot[1] +=  2*rhoNode*omega1*(vel(0,0,nodeNo) - (rNormTrue*omega1*tangent[0] + pipeVel[0]));  // NB sjekk denne

            forceRot[0] += -rhoNode*alpha_z*dy;
            forceRot[1] +=  rhoNode*alpha_z*dx;
            //------------------------------------------------------------------------- force translation
            std::valarray<lbBase_t> forceTra(0.0, 3);
            forceTra = rhoNode*pipeAks;

            forceRot += forceTra;

            forceRot *= heavsideInner;

            //------------------------------------------------------------------------- pressure force
            forceRot[2] += dp*(1-heavsideInner);

            solidFluidForce.set(0, nodeNo) = forceRot;


            //forceRot += 0*2*0.01*(rhoNode*rNormTrue*omega1*tangent - Mi - 0.5*forceRot)*heavisideStepReg<LT>(dist+2);

            const auto Mi = LT::qSumC(fNode);
            const std::valarray<lbBase_t> forceBoundary = 2*(rhoNode*velInner - Mi - 0.5*forceRot)*dirceDeltaReg(dist);
//            const auto forceBoundary = immersedBoundaryForce<LT>(dist, fNode, rhoNode, omegaInner, velInner);
            torqueX += -tangent[0]*forceBoundary[2]*rInner;
            torqueY += -tangent[1]*forceBoundary[2]*rInner;
            torqueZ += (tangent[0]*forceBoundary[0] + tangent[1]*forceBoundary[1])*rInner;
            forceX += forceBoundary[0];
            forceY += forceBoundary[1];
            forceZ += forceBoundary[2];

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
            velZ += velNode[2]*(1-heavsideInner);
            //------------------------------------------------------------------------- Calculation of the analytical velocity
            /* lbBase_t speedR;
            const lbBase_t eta = rInner/rOuter;
            const lbBase_t my = 0; // omegaOuter/omegaInner
            const lbBase_t A = omega_z*(my - eta*eta)/(1 - eta*eta);
            const lbBase_t B = omega_z*rInner*rInner*(1 - my)/(1 - eta*eta);
            if (rNormTrue >= rInner) {
                speedR = A*rNormTrue + B/rNormTrue;
                torqueAna += -(1.0/3.0)*(tauSym - 0.5)*rInner*dirceDeltaReg(rNormTrue - rInner)*(A - B/(rInner*rInner));
            } else {
                speedR = omega_z*rNormTrue;
                torqueAna += -(1.0/3.0)*(tauSym - 0.5)*rInner*dirceDeltaReg(rNormTrue - rInner)*(A - B/(rInner*rInner));
//                torqueAna += 0*(1.0/3.0)*(tauSym - 0.5)*rInner*dirceDeltaReg(rNormTrue - rInner)*omega1;
            }
            std::valarray<lbBase_t> velAnaNode(LT::nD);
            velAnaNode = speedR*tangent;
            velAna.set(0, nodeNo) = velAnaNode;
            */
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
//            tauSym = quemada.tau(fNode, rhoNode, velNode, u2, cu, forceNode);
//            tau = 0.502788344475428;
//            tau = 5055766889508577;    
            tauSym = tauSymOrg;
            // tauSym += (1- tauSym)*heavisideStepReg<LT>(dist + 4.5, 2.5);
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
        // surfaceForce = ippBndInner.applyTest(phi, rho, omegaInner, boundaryVelocity, 0, f, nodes, grid, solidFluidForce);

        lbBase_t torqueXGlobal, torqueYGlobal, torqueZGlobal;
        MPI_Allreduce(&torqueX, &torqueXGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&torqueY, &torqueYGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&torqueZ, &torqueZGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        lbBase_t torqueAnaGlobal;
        MPI_Allreduce(&torqueAna, &torqueAnaGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        lbBase_t forceXGlobal, forceYGlobal, forceZGlobal ;
        MPI_Allreduce(&forceX, &forceXGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&forceY, &forceYGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&forceZ, &forceZGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        lbBase_t velZGlobla;
        MPI_Allreduce(&velZ, &velZGlobla, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


        if (myRank == 0) {
            writeSurfaceForce << i << " " << pipePos[0] << " " << pipePos[1] << " ";
            writeSurfaceForce << pipeVel[0] << " " << pipeVel[1] << " " << pipeVel[2] << " ";
            writeSurfaceForce << omega_z << " " << dp << " " << velZGlobla << " ";
            writeSurfaceForce << forceXGlobal << " " << forceYGlobal << " " << forceZGlobal << " ";
            writeSurfaceForce << torqueXGlobal << " "  << torqueYGlobal << " "  << torqueZGlobal << std::endl;
            writeTorque << torqueZGlobal << " " << forceXGlobal << " " << forceYGlobal << std::endl;
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
                std::cout << "TORQUE = " << torqueZGlobal  <<  std::endl;
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
