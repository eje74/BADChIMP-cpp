//
// TO DO :
//        1) [OK] Define lattice information
//        2) Boundary nodes
//        3) What about ghost nodes. Do we need them?
//        3b) [OK] uansett ghost nodes info bør ligge i grid
//        4) [OK] Change function variables to longer descreptive names, (so we can easly see the difference input-variables and loop-iterators)
//        5) Se på denne: https://upload.wikimedia.org/wikipedia/commons/7/7d/OptimizingCpp.pdf
//        6) Legger bouncback etter hverandre
//         |b bHat
//
// GRID :
// /////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <vector>
#include "LBlattice.h"
#include "LBgrid.h"
#include "LBfield.h"
#include "LBboundary.h"
//#include "latticeboltzmann.h"





void setd2q9(Lattice &lattice)
{
    double w0 = 4.0 / 9.0;
    double w1 = 1.0 / 9.0;
    double w2 = 1.0 / 36.0;
    int cx[] = {1, 1, 0, -1, -1, -1,  0,  1, 0};
    int cy[] = {0, 1, 1,  1,  0, -1, -1, -1, 0};

    for (int q = 0; q < 8; q +=2) {
        lattice.setW(q, w1);
        lattice.setC(q, 0, cx[q]); lattice.setC(q, 1, cy[q]);
        lattice.setW(q+1, w2);
        lattice.setC(q+1, 0, cx[q+1]); lattice.setC(q+1, 1, cy[q+1]);
    }
    lattice.setC(8, 0, cx[8]); lattice.setC(8, 1, cy[8]);
    lattice.setW(8, w0);

}

template <typename BASETYPE>
inline BASETYPE LbEquilibirum(const int qDirection, const BASETYPE rho, const BASETYPE* vel, const Lattice &lattice)
{
    BASETYPE cu, uu;

    uu = lattice.dot(vel, vel);  // Only needed to be calulated once for each element, could be set at an input variable
    cu = lattice.cDot(qDirection, vel);

    return lattice.w(qDirection) * rho * ( 1.0 + lattice.c2Inv_ * cu + 0.5*lattice.c4Inv_ * (cu*cu - lattice.c2_*uu) );
}


bool insideDomain(int xNo, int yNo, int nX, int nY)
{
    return ( (xNo >= 0 ) && (xNo < nX) && (yNo >= 0) && (yNo < nY)   );
}

template<typename BASETYPE>
void makeGeometry(ScalarField<BASETYPE>& geo, GridRegular& grid, int nX, int nY)
{

    // Periodic
    for (int y = -1; y < nY+1; ++y) {
        for (int x = -1; x < nX+1; ++x)
            geo(0, grid.element(x, y)) = 2;
    }
    // Fluid
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
            geo(0, grid.element(x, y)) = 0;

    // Solid
    for (int x = 0; x < nX; ++x) {
        geo(0, grid.element(x, 0)) = 1;
    }

    // Fill ghost nodes with solid
    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            if ( (y==0) || (x==0) || (y == nY-1) || (x = nX-1) ) {
                if (geo(0, grid.element(x, y)) == 1) {
                    for (int q = 0; q < 8; ++q) {
                        int neigNode = grid.neighbor(q, grid.element(x, y));
                        int xNeig, yNeig;
                        grid.position(xNeig, yNeig, neigNode);
                        if (!insideDomain(xNeig, yNeig, nX, nY)) {
                            geo(0, grid.periodicNeighbor(q, grid.element(x, y))) = 1;
                        }
                    }
                }
            }
        }
    }
}

template<typename BASETYPE>
bool containsLinkType(int linkType, int nodeNo, ScalarField<BASETYPE>& geo, GridRegular& grid, Lattice& lattice)
{
    for (int q = 0; q < lattice.nQ() - 1; ++q) {
        if (geo(0, grid.neighbor(q, nodeNo)) == linkType)
                return true;
    }
    return false;
}


template<typename BASETYPE>
/* Find number of periodic nodes and bounce back */
void numberOfBoundaryNodes(int& nBounceBackNodes, int& nPeriodicNodes, ScalarField<BASETYPE>& geo, GridRegular& grid, Lattice& lattice, int nX, int nY)
{
    nBounceBackNodes = 0;
    nPeriodicNodes = 0;

    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            if (geo(0, grid.element(x, y)) == 0.0) {
                if (containsLinkType(1.0, grid.element(x, y), geo, grid, lattice))  nBounceBackNodes += 1;
                if (containsLinkType(2.0, grid.element(x, y), geo, grid, lattice))  nPeriodicNodes += 1;
            }
        }
    }
}


template<typename BASETYPE>
void addBoundaryNode(int node, int geoType, ScalarField<BASETYPE>& geo, Boundary& boundary, GridRegular& grid, Lattice& lattice)
{
    int nBeta=0, nGamma=0, nDelta=0;
    int betaList[4], gammaList[4], deltaList[4];
    bool linkPresent, reverseLinkPresent;

    if (geo(0, node) != 0)
        return;

    for (int q = 0; q < 4; q++) {
        linkPresent = false;
        reverseLinkPresent = false;

        if (geo(0, grid.neighbor(q, node)) == geoType)
            linkPresent = true;
        if (geo(0, grid.neighbor(lattice.reverseDirection(q), node)) == geoType)
            reverseLinkPresent = true;

        if ( !linkPresent && !reverseLinkPresent ) { // Gamma link
            gammaList[nGamma] = q;
            nGamma += 1;
        }
        else if ( linkPresent && reverseLinkPresent ) { // delta link
            deltaList[nDelta] = q;
            nDelta += 1;
        }
        else if ( !linkPresent && reverseLinkPresent ) { // beta = q
            betaList[nBeta] = q;
            nBeta += 1;
        }
        else if ( linkPresent && !reverseLinkPresent ) { // beta = lattice.reverse(q)
            betaList[nBeta] = lattice.reverseDirection(q);
            nBeta += 1;
        } else {
            printf("ERROR in addBoundary node");
            exit(1);
        }
    }

    boundary.addBoundaryNode( node,
                              nBeta,
                              betaList,
                              nGamma,
                              gammaList,
                              nDelta,
                              deltaList );

}

template<typename BASETYPE>
void setupBoundary(HalfWayBounceBack& bounceBackBoundary, Periodic& periodicBoundary, GridRegular& grid,  ScalarField<BASETYPE>& geo, Lattice& lattice, int nX, int nY)
{
    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            int node = grid.element(x, y);
            if (containsLinkType(1.0, node, geo, grid, lattice)) { // Bounce back link
                addBoundaryNode(node, 1.0, geo, bounceBackBoundary, grid, lattice);
            }
            if (containsLinkType(2.0, node, geo, grid, lattice)) { // Periodic link
                addBoundaryNode(node, 2.0, geo, periodicBoundary, grid, lattice);
            }
        }
    }
}


int main()
{
    std::cout << "Begin test";
    std::cout << std::endl;

    // Input data
    int nIterations;
    int nX, nY;
    int nPeriodicNodes, nBounceBackNodes;
    int nQ, nD;
    double tau, tau_inv;

    double force[2] = {1.0e-8, 0.0};
    double factor_force;
    nIterations = 1000;
    nX = 250; nY = 101;
    nQ = 9;
    nD = 2;
    tau = 1.0;
    tau_inv = 1.0 / tau;
    factor_force = (1 - 0.5 / tau);

    VectorField<double> my_vel(1, 2, 10);

    my_vel(0, 3)[1] = 4.0;

    // Simulation
    Lattice lattice(9, 2);
    setd2q9(lattice);    

    GridRegular grid(nX, nY, &lattice);
    ScalarField<int> geo(1, grid.nElements());
    makeGeometry(geo, grid, nX, nY);
    numberOfBoundaryNodes(nBounceBackNodes, nPeriodicNodes, geo, grid, lattice, nX, nY);

    LbField<double> f(1, nQ, grid.nElements()); // Bør f-field vite at det er nQ. Dette kan nok gjøre at ting blir gjort raskere
    LbField<double> fTmp(1, nQ, grid.nElements()); // Do not need to be a LbField (as f), since its just a memory holder, that is immediately swaped after propagation.
    ScalarField<double> rho(1, grid.nElements()); // Dette kan fikses for 'single rho fields' slik at man slipper å skrive rho(0, pos)
    VectorField<double> vel(1, nD, grid.nElements()); // Bør vel-field vite at det er 2D. Dette kan nok gjøre at ting blir gjort raskere

    HalfWayBounceBack halfWayBounceBack(nBounceBackNodes, 4);
    Periodic periodic(nPeriodicNodes, 4);

    // SETUP BOUNDARY
    setupBoundary(halfWayBounceBack, periodic, grid, geo, lattice, nX, nY);

    // INIT:
    double velTmp[2] = {0.0, 0.0};
    //VectorField<double> force(1, 2, 1);
    for (int y = 0; y < nY; y++ ) {
        for (int x = 0; x < nX; x++) {
            for (int q = 0; q < nQ; q++) {
                f(0, q, grid.element(x, y)) = LbEquilibirum(q, 1.0, velTmp, lattice);
            }
        }
    }

    // standard.collitionpropagate()
    // pressourboundarary.applyyBounrayCodtion(f, gird, latitice)
    // pressutboudnar.colltionpropagte()
    // - LOOP TYPE 1:
    for (int i = 0; i < nIterations; i++) {
   //     std::cout << "@ iteration " << i << std::endl;
        for (int y = 1; y < nY; y++ ) {
            for (int x = 0; x < nX; x++) {
                int nodeNo = grid.element(x, y);
                double* fNode = f(0, nodeNo);
                double velNode[2];
                double rhoNode;
                // * Macrosopics
                // * * rho and vel
                rhoNode = lattice.qSum(fNode);
                velNode[0] = (lattice.qSumC(0, fNode) + 0.5 * force[0]) / rhoNode;
                velNode[1] = (lattice.qSumC(1, fNode) + 0.5 * force[0]) / rhoNode;
                rho(0, nodeNo) = rhoNode;
                vel(0, 0, nodeNo) = velNode[0];
                vel(0, 1, nodeNo) = velNode[1];
                //f.zerothMoment(nodeNo, &rho(0,nodeNo));
                //f.firstMoment(nodeNo, vel(0, nodeNo));
                /* for (int d = 0; d < lattice.nD(); ++d) {
                    vel(0, d, nodeNo) += 0.5 * force[d];
                    vel(0, d, nodeNo) /= rho(0, nodeNo);
                } */

                // * Collision and propagation:
                for (int q = 0; q < nQ; q++) {  // Collision should provide the right hand side must be
                                                // collision(q, f(q, n), rho(n), u(n) \\ should be an array) ?
                    fTmp(0, q, grid.neighbor(q, nodeNo)) = f(0, q, nodeNo) - tau_inv * ( f(0, q, nodeNo) - LbEquilibirum(q, rhoNode, velNode, lattice) );
                    fTmp(0, q, grid.neighbor(q, nodeNo)) += lattice.w(q) * factor_force * (lattice.c2Inv_ * lattice.cDot(q, force)
                                                                                              + lattice.c4Inv_ *
                                                                                              (
                                                                                                  lattice.cDot(q, force) * lattice.cDot(q, velNode) - lattice.c2_ * lattice.dot(velNode, force))
                                                                                              );
                } // End collision and propagation
            }

        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);

        // BOUNDARY CONDITION Half way bounce back
        periodic.applyBoundaryCondition(f, grid, lattice);
        halfWayBounceBack.applyBoundaryCondition(f, grid, lattice);

    } // End iterations (LOOP TYPE 1)

    std::cout << std::setprecision(3) << vel(0, 0, grid.element(10, 10))  << std::endl;

    // Write to screen
    /* for (int y = nY-1; y >= 0; --y ) {
        for (int x = 0; x < nX; ++x) {
            std::cout << std::setprecision(3) << rho(0, grid.element(x, y)) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;

    for (int y = nY-1; y >= 0; --y ) {
        for (int x = 0; x < nX; ++x) {
            std::cout << std::setprecision(3) << vel(0, grid.element(x, y)) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;

    for (int y = nY-1; y >= 0; --y ) {
        for (int x = 0; x < nX; ++x) {
            std::cout << std::setprecision(3) << vel(1, grid.element(x, y)) << " ";
        }
        std::cout << std::endl;
    } */



    return 0;
}

