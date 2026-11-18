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
// #define LT D3Q19

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
    int nIterations = input["iterations"]["max"];
    // Write interval
    int nItrWrite = input["iterations"]["write"];
    // Relaxation time
    lbBase_t tau = input["fluid"]["tau"];
    // Fluid bodyforce driver
    // mean inlet velocity
    VectorField<LT> uMean(1, 1, input["driver"]["uinlet"]);
    // rampup time
    int rampupItr = input["driver"]["rampup"];

    // Data analysis
    std::vector<int> yLim{input["data"]["bottom_y"], input["data"]["top_y"]};
    std::vector<int> inletNodes(yLim[1] - yLim[0]);
    int cnt = 0;
    for (int ypos=yLim[0]; ypos < yLim[1]; ++ypos) {
        std::vector<int> pos{0, ypos};
        inletNodes[cnt] = grid.nodeNo(pos);
        const auto xy = grid.pos(inletNodes[cnt]);
        std::cout << xy[0] << " " << xy[1] << "\n";
        cnt++;
    }

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
    VectorField<LT> bf(1, grid.size());
    // Initiate velocity
    for (auto nodeNo: bulkNodes) {
        for (int d=0; d < LT::nD; ++d) {
            vel(0, d, nodeNo) = 0.0;
            bf(0, d, nodeNo) = 0.0;
        }
    }

    // Bodyforce
    VectorField<LT> bodyForce(1, grid.size());
    for (int d=0; d < LT::nD; ++d)
        bodyForce(0, d, 0) = 0.0;

    // Set mean velocity
    VectorField<LT> uMeanCurrent(1, grid.size());
    for (int d=0; d < LT::nD; ++d)
        uMeanCurrent(0, d, 0) = 0.0;

    // Particle indicator
    ScalarField gamma(1, grid.size());
    vtklb.toAttribute("gamma");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        gamma(0, n) = vtklb.getScalarAttribute<lbBase_t>();
    }

    // Particle surface indicator
    ScalarField delta(1, grid.size());
    vtklb.toAttribute("delta");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        delta(0, n) = vtklb.getScalarAttribute<lbBase_t>();
    }

    // Particle normal vector
    VectorField<LT> normVec(1, grid.size());
    vtklb.toAttribute("Nx");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        normVec(0, 0, n) = vtklb.getScalarAttribute<lbBase_t>();
    }
    vtklb.toAttribute("Ny");
    for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
        normVec(0, 1, n) = vtklb.getScalarAttribute<lbBase_t>();
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
    // initiate lb distributions
    for (auto nodeNo: bulkNodes) {
        for (int q = 0; q < LT::nQ; ++q) {
            f(0, q, nodeNo) = LT::w[q]*rho(0, nodeNo);
        }
    }

    // **********
    // OUTPUT VTK
    // **********
    Output<LT> output(grid, bulkNodes, outputDir, myRank, nProcs);
    output.add_file("lb_run");
    output.add_scalar_variables({"rho"}, {rho});
    output.add_vector_variables({"vel", "bf", "normVec"}, {vel, bf, normVec}); 
    output.add_scalar_variables({"gamma"}, {gamma});

    // **********
    // INIT WRITE
    // **********
    std::cout << "SETUP:\n";
    std::cout << "\tnumber of grid point = " << grid.size() << "\n";
    std::cout << "\tmean inlet velocity = [" << uMean(0, 0, 0) << ", " << uMean(0, 1, 0) << "]\n";
    std::cout << "\trampup time = " << rampupItr << "\n";
    std::cout << "\tylim = [" << yLim[0] << ", " << yLim[1] << "]\n";


    // *********
    // MAIN LOOP
    // *********
    for (int i = 0; i <= nIterations; i++) {
        // Ramp up mean velocity
        if (i > rampupItr) {
            for (int d=0; d < LT::nD; ++d)
                uMeanCurrent(0, d, 0) = uMean(0, d, 0);
        }
        else {
            double tmp = 0.5 * (1. - std::cos(i *3.14159 / rampupItr));
            for (int d=0; d < LT::nD; ++d)
                uMeanCurrent(0, d, 0) = tmp * uMean(0, d, 0);
        }

        // Calculate velocity forcing
        lbBase_t cnt = 0.0;
        lbBase_t gammaMxMean = 0.0;
        lbBase_t gammaMean = 0.0;
        for (auto nodeNo: bulkNodes) {
            const std::valarray<lbBase_t> fNode = f(0, nodeNo);
            const auto mi = LT::qSumC(fNode);
            const auto gammaNode = gamma(0, nodeNo);
            gammaMean += (1-gammaNode) * (1-gammaNode);
            gammaMxMean += (1-gammaNode) * mi[0];
            cnt += 1.0;
        }
        gammaMxMean /= cnt;
        gammaMean /= cnt;
        const double Fx = (2.0 / gammaMean) * (uMeanCurrent(0, 0, 0) - gammaMxMean);

        lbBase_t sumFs = 0.0;
        lbBase_t surfFa = 0.0;
        lbBase_t surfFb = 0.0;
        lbBase_t sumBulk = 0.0;

        for (auto nodeNo: bulkNodes) {
            // Copy of local velocity diestirubtion
            const std::valarray<lbBase_t> fNode = f(0, nodeNo);
            const auto mi = LT::qSumC(fNode);
            const lbBase_t rhoNode = calcRho<LT>(fNode);
            // Macroscopic values
            const lbBase_t deltaRho = rhoNode - 1.0;


            const lbBase_t Fsx = 2.0 * gamma(0, nodeNo) * (-mi[0] );
            const lbBase_t Fsy = 2.0 * gamma(0, nodeNo) * (-mi[1] );



            sumFs -= Fsx;
            surfFa -=  delta(0, nodeNo) * Fsx;
            surfFb -= (1 - gamma(0, nodeNo)) * Fsx;
            if (gamma(0, nodeNo) > 0.9999999999)
                sumBulk -= Fsx;
            bodyForce(0, 0, 0) = (1 - gamma(0, nodeNo)) * Fx + Fsx;
            bodyForce(0, 1, 0) =  Fsy;
            bf(0, 0, nodeNo) = Fsx;
            bf(0, 1, nodeNo) = Fsy;


            const auto xy = grid.pos(nodeNo);
            const lbBase_t yNorm = 1.0 / (43. - 2.);
            if (xy[0] == 0) {
                const lbBase_t y = (xy[1] - 0.5) * yNorm;
                const lbBase_t rhoU = 6 * uMeanCurrent(0, 0, 0) * y *(1-y) * rhoNode;
                bodyForce(0, 0, 0) = 0.68*2*(rhoU - mi[0]) + 0.32*bodyForce(0, 0, 0);
                bodyForce(0, 1, 0) = 0.68*2*(-mi[1]);
            }
            if (xy[0] == 1) {
                const lbBase_t y = (xy[1] - 0.5) * yNorm;
                const lbBase_t rhoU = 6 * uMeanCurrent(0, 0, 0) * y *(1-y) * rhoNode;
                bodyForce(0, 0, 0) = 0.16*2*(rhoU - mi[0]) + 0.16*bodyForce(0, 0, 0);
                bodyForce(0, 1, 0) = 0.16*2*(-mi[1]);
            }
            if (xy[0] == (220 - 1)) {
                const lbBase_t y = (xy[1] - 0.5) * yNorm;
                const lbBase_t rhoU = 6 * uMeanCurrent(0, 0, 0) * y *(1-y) * rhoNode;
                bodyForce(0, 0, 0) = 0.16*2*(rhoU - mi[0]) + 0.16*bodyForce(0, 0, 0);
                bodyForce(0, 1, 0) = 0.16*2*(-mi[1]);
            }


            const auto velNode = calcVel<LT>(fNode, rhoNode, bodyForce(0, 0));

            // Save density and velocity for printing
            rho(0, nodeNo) = rhoNode;
            vel.set(0, nodeNo) = velNode;
                
            // BGK-collision term
            const lbBase_t u2 = LT::dot(velNode, velNode);
            const std::valarray<lbBase_t> cu = LT::cDotAll(velNode);
            const std::valarray<lbBase_t> omegaBGK = calcOmegaBGK<LT>(fNode, tau, rhoNode, u2, cu);
            
            // Calculate the Guo-force correction
            const lbBase_t uF = LT::dot(velNode, bodyForce(0, 0));
            const std::valarray<lbBase_t> cF = LT::cDotAll(bodyForce(0, 0));
            const std::valarray<lbBase_t> deltaOmegaF = calcDeltaOmegaF<LT>(tau, cu, uF, cF);

            // Collision and propagation
            const lbBase_t sc = 1.0 / (1.0 - LT::w0);
            std::valarray<lbBase_t> val(LT::nQ);
            for (int q = 0; q < LT::nQ - 1; ++q) {
                val[q] += LT::w[q]*deltaRho*sc;
            }
            val[8] = 0.0;

            fTmp.propagateTo(0, nodeNo, fNode + omegaBGK + deltaOmegaF + val, grid);

        } // End nodes

        f.swapData(fTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(0, f, grid);
        // Half way bounce back
        bounceBackBnd.apply(f, grid);

        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0)  ) {
            output.write(i);
            if (myRank==0) {
                const lbBase_t Dlb = 10;
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
                std::cout << "\tu = [" << uMeanCurrent(0, 0, 0) << ", " << uMeanCurrent(0, 1, 0) << "]\n";
                std::cout << "\t sum Fs = " << sumFs << "  " <<  sumFs / (0.5* uMean(0, 0, 0) * uMean(0, 0, 0) * Dlb) << "\n";
                std::cout << "\t surf Fa = " << surfFa << "  " <<  surfFa / (0.5* uMean(0, 0, 0) * uMean(0, 0, 0) * Dlb) << "\n";
                std::cout << "\t surf Fb = " << surfFb << "  " <<  surfFb / (0.5* uMean(0, 0, 0) * uMean(0, 0, 0) * Dlb) << "\n";
                std::cout << "\t sum Bulk = " << sumBulk << "  " <<  sumBulk / (0.5* uMean(0, 0, 0) * uMean(0, 0, 0) * Dlb) << "\n";
            }
        }

    } // End iterations

    MPI_Finalize();

    return 0;
}
