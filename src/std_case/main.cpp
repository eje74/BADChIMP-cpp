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
#include <cmath>  // for sqrt

// SET THE LATTICE TYPE
#define LT D2Q9
// #define LT D3Q19

#define EPS 1.0e-14

// Include array list
#include <array>
std::vector<std::array<int, 3>> makeBoundaryList(int indicator, const std::vector<int> & boundaryIndicator, const std::vector<int> & bulkNodes, const Grid<LT> &grid)
/*  fluid node, boundary node, direction 
 *  direction is direction from fluid node to boundary node 
 *
 *  geo : bulk fluid
*/
{
    std::vector<std::array<int, 3>> ret;

    for (auto nodeNo: bulkNodes) {
        if (boundaryIndicator[nodeNo] == 0) {
            for (int q = 0; q < LT::nQ; ++q) {
                int neigNo = grid.neighbor(q, nodeNo);
                if (boundaryIndicator[neigNo] == indicator) {
                    ret.push_back({nodeNo, neigNo, q});
                }
            }
        }
    }

    return ret;
}

void inletBoundary(LbField<LT> &f, const std::vector<std::array<int, 3>> &inlet, lbBase_t CInlet, lbBase_t HInlet, const ScalarField &rho, const Grid<LT> &grid)
{
    lbBase_t yNorm = 1.0 / HInlet;

    for (auto & el: inlet) {
        const int fluidNo = el[0];
        const int wallNo = el[1];
        const int alpha = el[2];
        const int alphaHat = LT::reverseDirection(alpha);
        const lbBase_t y = yNorm * (grid.pos(el[0], 1) - 0.5 + 0.5*LT::c(alpha, 1)); 
        const std::vector<lbBase_t> vel = {CInlet * y * (1.0 - y), 0.0};

        // Set the bounce back distribution 
        f(0, alphaHat, fluidNo) = f(0, alpha, wallNo) + 2.0 * LT::w[alphaHat] * rho(0, fluidNo) * LT::c2Inv * LT::cDotRef(alphaHat, vel);
    }
}

void outletBoundary(LbField<LT> &f, const std::vector<std::array<int, 3>> &outlet, const VectorField<LT> & vel ,const Grid<LT> &grid)
{
    lbBase_t rhoOutlet = 1.0;

    for (auto & el: outlet) {
        const int fluidNo = el[0];
        const int wallNo = el[1];
        const int alpha = el[2];
        const int alphaHat = LT::reverseDirection(alpha);

        const std::valarray<lbBase_t> velNode = vel(0, fluidNo);

        const lbBase_t u2 = LT::dot(velNode, velNode);
        const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);

        f(0, alphaHat, fluidNo) = 2.0 * LT::w[alphaHat] * rhoOutlet * (1.0 + LT::c4Inv0_5*(cu[alphaHat]*cu[alphaHat] - LT::c2*u2)) - f(0, alpha, wallNo);
    }
}

// *******************
// PARTICLE SIMULATION
// *******************
lbBase_t gamma_part(const std::vector<int> & posNode, const std::valarray<lbBase_t> & posPart, lbBase_t dPart)
{
    const lbBase_t eps = 1.0;
    const lbBase_t dx = posNode[0] - posPart[0];
    const lbBase_t dy = posNode[1] - posPart[1];
    lbBase_t r = std::sqrt(dx*dx + dy*dy);
    lbBase_t d = r - 0.5*dPart;

    if (d > eps ) return 0.0;
    if (d < -eps ) return 1.0;
    return 0.5 * (1 - std::sin(0.5*3.14159 * d / eps));
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
    std::string chimpDir = "./../";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

    // ***********************
    // SETUP GRID AND GEOMETRY
    // ***********************
    Input input(inputDir + "input_particle.dat");
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
    int nIterations = input["iterations"]["max"];
    // Write interval
    int nItrWrite = input["iterations"]["write"];
    // Relaxation time
    lbBase_t tau = input["fluid"]["tau"];
    // Driver force
    const lbBase_t CInlet_org = 6.0 * input["driver"]["uMean"];
    const lbBase_t HInlet = 41.0;
    const lbBase_t LQ = 0.0;


    // ******************
    // MACROSCOPIC FIELDS
    // ******************
    // Density
    ScalarField rho(1, grid.size());
    // Initiate density from file
    vtklb.toAttribute("init_rho");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        rho(0, n) = vtklb.getScalarAttribute<lbBase_t>();
    }

    // Velocity
    VectorField<LT> vel(1, grid.size());
    // Initiate velocity
    for (auto nodeNo: bulkNodes) {
        for (int d=0; d < LT::nD; ++d) {
            vel(0, d, nodeNo) = 0.0;
        }
    }

    // Bodyforce
    VectorField<LT> bodyForce(1, grid.size());
    for (int d=0; d < LT::nD; ++d)
        bodyForce(0, d, 0) = 0.0;

    // Particle indicator
    ScalarField gamma(1, grid.size());
    vtklb.toAttribute("gamma");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        gamma(0, n) = vtklb.getScalarAttribute<lbBase_t>();
        // gamma(0, n) = 0.0;
    }

    // Boundary indicator
    std::vector<int> boundaryIndicator(grid.size(), 0);
    vtklb.toAttribute("boundary");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        boundaryIndicator[n] = vtklb.getScalarAttribute<lbBase_t>();
    }

    auto inlet = makeBoundaryList(1, boundaryIndicator, bulkNodes, grid);
    auto outlet = makeBoundaryList(2, boundaryIndicator, bulkNodes, grid);
    
    // ******************
    // SETUP PARTICLE
    // ******************
    lbBase_t diamPart = 10.0;
    lbBase_t massPart = 3.141592*0.25*diamPart*diamPart;
    
    std::valarray<lbBase_t> velPartOld(LT::nD), velPartNew(LT::nD);
    std::valarray<lbBase_t> posPart(LT::nD), velPart(LT::nD), accPart(LT::nD); 
    posPart = {15.0, 15.0};
    velPartNew = {0.0, 0.0};
    velPartOld = {0.0, 0.0};
    velPart = {0.0, 0.0};




    // ******************
    // SETUP BOUNDARY
    // ******************
    HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

    // *********
    // LB FIELDS
    // *********
    LbField<LT> f(1, grid.size()); 
    LbField<LT> fTmp(1, grid.size());
    // initiate lb distributions
    for (auto nodeNo: bulkNodes) {
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = LT::w[q]*rho(0, nodeNo);
        }
    }

    // **********
    // OUTPUT VTK
    // **********
    VectorField<LT> boundaryForce(1, grid.size());
    Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs);
    output.add_file("lb_run");
    output.add_scalar_variables({"rho"}, {rho});
    output.add_vector_variables({"vel", "boundary_force"}, {vel, boundaryForce}); 
    output.add_scalar_variables({"gamma"}, {gamma});

    // **********
    // INIT WRITE
    // **********
    std::cout << "SETUP:\n";
    std::cout << "\tnumber of grid point = " << grid.size() << "\n";


    // *********
    // MAIN LOOP
    // *********

    // *** Force calculation
    std::valarray<lbBase_t> force_contour(0.0, LT::nD);
    std::valarray<lbBase_t> force_boundary(0.0, LT::nD);
    std::valarray<lbBase_t> force_star(0.0, LT::nD);
    std::valarray<lbBase_t> force_star_old(0.0, LT::nD);

    for (int i = 0; i <= nIterations; i++) {

        // // Ramp-up
        // lbBase_t CInlet = 0.001 * i * CInlet_org;
        // if (CInlet > CInlet_org) CInlet = CInlet_org;

        lbBase_t CInlet = CInlet_org;


        force_star_old = force_star;
        force_star = 0.0;
        force_contour = 0.0;
        force_boundary = 0.0;

        for (auto nodeNo: bulkNodes) {
            // Copy of local velocity diestirubtion
            const std::valarray<lbBase_t> fNode = f(0, nodeNo);
            const auto mi = LT::qSumC(fNode);
            const lbBase_t rhoNode = calcRho<LT>(fNode);


            // Set solid velocity
            gamma(0, nodeNo) = gamma_part(grid.pos(nodeNo), posPart, diamPart);
            const lbBase_t gammaNode = gamma(0, nodeNo);
            std::valarray<lbBase_t> forceNode = (1.0 - (gammaNode>0.0)) * bodyForce(0, 0);
            const std::valarray<lbBase_t> forceSolidBoundary =  gammaNode * 2 * ( velPart - mi); 
            forceNode += forceSolidBoundary;

            boundaryForce.set(0, nodeNo) = forceSolidBoundary;

            force_contour += (gammaNode > 1.0 - EPS) * mi;

            force_boundary += -(gammaNode < 1.0 - EPS) * forceSolidBoundary;

            // Macroscopic values
            const auto velNode = calcVel<LT>(fNode, rhoNode, forceNode);

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

            const auto mi_star = LT::qSumC(fNode + omegaBGK + deltaOmegaF);

            force_star += (gammaNode > 1.0 - EPS) * mi_star;    

            fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF, grid);

        } // End nodes

        f.swapData(fTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(0, f, grid);
        // Half way bounce back
        bounceBackBnd.apply(f, grid);
        // Inlet
        inletBoundary(f, inlet, CInlet, HInlet, rho, grid);
        outletBoundary(f, outlet, vel, grid);

        // *******************
        // update particle position
        // *******************
        if (i > 0) {
            accPart = (force_boundary + force_contour - force_star_old) / massPart;
            velPartNew = velPartOld + accPart;
            posPart = posPart + velPart;

            velPart = 0.5 * (velPartNew + velPartOld);
            velPartOld = velPartNew;
        }

        

        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0)  ) {
            output.write(i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
                std::cout << "\t velocity = ";
                for (auto el: velPart) {
                    std::cout << el << " ";
                }
                std::cout << std::endl;

                std::valarray<lbBase_t> vtmp = force_boundary + force_contour - force_star_old;
                std::cout << "\t force (best estimate) = ";
                for (auto el: vtmp) {
                    std::cout << el << " ";
                }
                std::cout << std::endl;
            }
        }

    } // End iterations

    MPI_Finalize();

    return 0;
}
