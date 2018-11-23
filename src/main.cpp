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
#include "LBd2q9.h"
#include "LBgrid.h"
#include "LBfield.h"
#include "LBboundary.h"


template <typename DXQY>
inline void LbEqAll(const lbBase_t tau_inv, const lbBase_t *f, const lbBase_t rho, const lbBase_t* cu, const lbBase_t uu, lbBase_t* fEqAll)
{
    for (int q = 0; q < DXQY::nQ; ++q)
        fEqAll[q] = f[q] + tau_inv * (DXQY::w[q] * rho * (1 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*uu) ) - f[q]);
}

template <typename DXQY>
inline void LbForceAll(const lbBase_t tau_factor, const lbBase_t* cu, const lbBase_t* cf, const lbBase_t uf, lbBase_t* forceAll)
{
    for (int q = 0; q < DXQY::nQ; ++q)
        forceAll[q] = DXQY::w[q]*tau_factor * (DXQY::c2Inv*cf[q] + DXQY::c4Inv * ( cf[q] * cu[q] - DXQY::c2 * uf));
}

bool insideDomain(int xNo, int yNo, int nX, int nY)
{
    return ( (xNo >= 0 ) && (xNo < nX) && (yNo >= 0) && (yNo < nY)   );
}

template <typename DXQY>
void makeGeometry(GeoField& geo, GridRegular<DXQY>& grid, int nX, int nY)
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

template <typename DXQY>
bool containsLinkType(int linkType, int nodeNo, GeoField& geo, GridRegular<DXQY>& grid)
{
    for (int q = 0; q < DXQY::nQ - 1; ++q) {
        if (geo(0, grid.neighbor(q, nodeNo)) == linkType)
            return true;
    }
    return false;
}



/* Find number of periodic nodes and bounce back */
template <typename DXQY>
void numberOfBoundaryNodes(int& nBounceBackNodes, int& nPeriodicNodes, GeoField& geo, GridRegular<DXQY>& grid, int nX, int nY)
{
    nBounceBackNodes = 0;
    nPeriodicNodes = 0;

    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            if (geo(0, grid.element(x, y)) == 0.0) {
                if (containsLinkType(1.0, grid.element(x, y), geo, grid))  nBounceBackNodes += 1;
                if (containsLinkType(2.0, grid.element(x, y), geo, grid))  nPeriodicNodes += 1;
            }
        }
    }
}

template <typename DXQY>
void addBoundaryNode(int node, int geoType, GeoField& geo, Boundary& boundary, GridRegular<DXQY>& grid)
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
        if (geo(0, grid.neighbor(DXQY::reverseDirection(q), node)) == geoType)
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
            betaList[nBeta] = DXQY::reverseDirection(q);
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


template <typename DXQY>
void setupBoundary(HalfWayBounceBack& bounceBackBoundary, Periodic& periodicBoundary, GridRegular<DXQY>& grid,  GeoField& geo, int nX, int nY)
{
    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            int node = grid.element(x, y);
            if (containsLinkType(1.0, node, geo, grid)) { // Bounce back link
                addBoundaryNode(node, 1.0, geo, bounceBackBoundary, grid);
            }
            if (containsLinkType(2.0, node, geo, grid)) { // Periodic link
                addBoundaryNode(node, 2.0, geo, periodicBoundary, grid);
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
    double tau, tau_inv;

    double force[2] = {1.0e-8, 0.0};
    double factor_force;
    nIterations = 10000;
    nX = 250; nY = 101;

    tau = 0.7;
    tau_inv = 1.0 / tau;
    factor_force = (1 - 0.5 / tau);


    GridRegular<D2Q9> grid(nX, nY);
    GeoField geo(1, grid.nElements());
    makeGeometry(geo, grid, nX, nY);
    numberOfBoundaryNodes(nBounceBackNodes, nPeriodicNodes, geo, grid, nX, nY);

    LbField f(1, D2Q9::nQ, grid.nElements()); // Bør f-field vite at det er nQ. Dette kan nok gjøre at ting blir gjort raskere
    LbField fTmp(1, D2Q9::nQ, grid.nElements()); // Do not need to be a LbField (as f), since its just a memory holder, that is immediately swaped after propagation.
    ScalarField rho(1, grid.nElements()); // Dette kan fikses for 'single rho fields' slik at man slipper å skrive rho(0, pos)
    VectorField vel(1, D2Q9::nD, grid.nElements()); // Bør vel-field vite at det er 2D. Dette kan nok gjøre at ting blir gjort raskere

    HalfWayBounceBack halfWayBounceBack(nBounceBackNodes, 4);
    Periodic periodic(nPeriodicNodes, 4);

    // SETUP BOUNDARY
    setupBoundary(halfWayBounceBack, periodic, grid, geo, nX, nY);

    // INIT:
    lbBase_t velTmp[2] = {0.0, 0.0};
    //VectorField<double> force(1, 2, 1);
    for (int y = 0; y < nY; y++ ) {
        for (int x = 0; x < nX; x++) {
            lbBase_t cu[D2Q9::nQ], uu;
            D2Q9::cDotAll(velTmp, cu);
            uu = D2Q9::dot(velTmp, velTmp);
            for (int q = 0; q < D2Q9::nQ; q++) {
                f(0, q, grid.element(x, y)) = D2Q9::w[q] * (1.0 + D2Q9::c2Inv*cu[q] + D2Q9::c4Inv0_5*(cu[q]*cu[q] - D2Q9::c2*uu));
            }
        }
    }


    // standard.collitionpropagate()
    // pressourboundarary.applyyBounrayCodtion(f, gird, latitice)
    // pressutboudnar.colltionpropagte()
    // - LOOP TYPE 1:
    for (int i = 0; i < nIterations; i++) {
        for (int y = 1; y < nY; y++ ) {
            for (int x = 0; x < nX; x++) {
                int nodeNo = grid.element(x, y);
                double velNode[2];
                double rhoNode;
                // * Macrosopics
                // * * rho and vel
                D2Q9::qSum(f(0, nodeNo), rhoNode);
                D2Q9::qSumC(f(0, nodeNo), velNode);
                for (int d = 0; d < D2Q9::nD; ++d) {
                    velNode[d] += 0.5 * force[d];
                    velNode[d] /= rhoNode;
                }
                rho(0, nodeNo) = rhoNode;
                vel(0, 0, nodeNo) = velNode[0];
                vel(0, 1, nodeNo) = velNode[1];

                lbBase_t uu, uF;
                uu = D2Q9::dot(velNode, velNode);
                uF = D2Q9::dot(velNode, force);

                // * Collision and propagation:
                lbBase_t  cul[D2Q9::nQ];
                D2Q9::cDotAll(velNode, cul);

                lbBase_t fEql[D2Q9::nQ];
                LbEqAll<D2Q9>(tau_inv, f(0, nodeNo), rhoNode, cul, uu, fEql);

                lbBase_t  cfl[D2Q9::nQ];
                D2Q9::cDotAll(force, cfl);

                lbBase_t forcel[D2Q9::nQ];
                LbForceAll<D2Q9>(factor_force, cul, cfl, uF, forcel);

                for (int q = 0; q < D2Q9::nQ; q++) {  // Collision should provide the right hand side must be
                    fTmp(0, q,  grid.neighbor(q, nodeNo)) = fEql[q] + forcel[q];//fTmp(0, q, grid.neighbor(q, nodeNo)) = fEql[q] + forcel[q];
                } // End collision and propagation */
            }

        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);

        // BOUNDARY CONDITION Half way bounce back
        periodic.applyBoundaryCondition(f, grid);
        halfWayBounceBack.applyBoundaryCondition(f, grid);

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

