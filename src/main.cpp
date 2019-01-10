//
// TO RUN PROGRAM: type bin/runner in command line in main directory
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
#include "LBbulk.h"
#include "LBgeometry.h"


// COLLISION
template <typename DXQY>
inline void LbEqAll(const lbBase_t tau_inv, const lbBase_t *f, const lbBase_t rho, const lbBase_t* cu, const lbBase_t uu, lbBase_t* fEqAll)
{
  for (int q = 0; q < DXQY::nQ; ++q)
    fEqAll[q] = f[q] + tau_inv * (DXQY::w[q] * rho * (1 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*uu) ) - f[q]);
}

// COLLISION
template <typename DXQY>
inline void LbForceAll(const lbBase_t tau_factor, const lbBase_t* cu, const lbBase_t* cf, const lbBase_t uf, lbBase_t* forceAll)
{
    for (int q = 0; q < DXQY::nQ; ++q)
        forceAll[q] = DXQY::w[q]*tau_factor * (DXQY::c2Inv*cf[q] + DXQY::c4Inv * ( cf[q] * cu[q] - DXQY::c2 * uf));
}


// NOT TO BE PART OF LB CORE
bool insideDomain(int xNo, int yNo, int nX, int nY)
{
    return ( (xNo >= 0 ) && (xNo < nX) && (yNo >= 0) && (yNo < nY)   );
}


// NOT TO BE PART OF LB CORE => IMPORT GEOMETRY
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
    std::cout<<"No. bounce back nodes: "<<nBounceBackNodes<<std::endl;
    std::cout<<"No. periodic nodes: "<<nPeriodicNodes<<std::endl;
}


/* Find number of bulk nodes */
template <typename DXQY>
void numberOfBulkNodes(int& nBulkNodes, GeoField& geo, GridRegular<DXQY>& grid, int nX, int nY)
{
    nBulkNodes = 0;

    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            if (geo(0, grid.element(x, y)) == 0.0) {
	          nBulkNodes += 1;
                
            }
        }
    }

    std::cout<<"No. bulk nodes: "<<nBulkNodes<<std::endl;
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

template <typename DXQY>
void setupBulk(Bulk& bulk, GridRegular<DXQY>& grid,  GeoField& geo, int nX, int nY)
{
    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            if (geo(0, grid.element(x, y)) == 0.0) { // Bulk node
                bulk.addBulkNode(grid.element(x, y));
            }
        }
    }
}


inline void updateMacroscopicFields(const lbBase_t* dist, lbBase_t& rhoNode, lbBase_t& rho, lbBase_t* velNode, lbBase_t* vel, lbBase_t* force){
  D2Q9::qSum(dist, rhoNode);
  D2Q9::qSumC(dist, velNode);
  for (int d = 0; d < D2Q9::nD; ++d) {
    velNode[d] = (velNode[d] + 0.5 * force[d]) /rhoNode;
    vel[d] = velNode[d];
  }
  
  rho = rhoNode;
}

inline void collision(const lbBase_t tau_inv, const lbBase_t* dist, const lbBase_t* velNode, const lbBase_t rhoNode,
		      const lbBase_t* force,const lbBase_t factor_force, lbBase_t* fEql, lbBase_t* forcel){
  
  lbBase_t uu, uF;
  uu = D2Q9::dot(velNode, velNode);
  uF = D2Q9::dot(velNode, force);
  
  
  
  // * Collision and propagation:
  lbBase_t  cul[D2Q9::nQ];
  D2Q9::cDotAll(velNode, cul);
  
  //lbBase_t fEql[D2Q9::nQ];
  LbEqAll<D2Q9>(tau_inv, dist, rhoNode, cul, uu, fEql);
  
  lbBase_t  cfl[D2Q9::nQ];
  D2Q9::cDotAll(force, cfl);
  
  //lbBase_t forcel[D2Q9::nQ];
  LbForceAll<D2Q9>(factor_force, cul, cfl, uF, forcel);

  
  
}

int main()
{
    std::cout << "Begin test";
    std::cout << std::endl;

    // Input data
    int nIterations;
    int nX, nY;
    int nPeriodicNodes, nBounceBackNodes, nBulkNodes;
    lbBase_t tau, tau_inv;

    lbBase_t force[2] = {1.0e-8, 0.0};
    lbBase_t factor_force;
    nIterations = 1;// 10000;
    nX = 250; nY = 101;

    tau = 0.7;
    tau_inv = 1.0 / tau;
    factor_force = (1 - 0.5 / tau);


    GridRegular<D2Q9> grid(nX, nY);
    GeoField geo(1, grid.nElements());
    makeGeometry(geo, grid, nX, nY);
    numberOfBoundaryNodes(nBounceBackNodes, nPeriodicNodes, geo, grid, nX, nY);
    numberOfBulkNodes(nBulkNodes, geo, grid, nX, nY);

    LbField f(1, D2Q9::nQ, grid.nElements()); // Bør f-field vite at det er nQ. Dette kan nok gjøre at ting blir gjort raskere
    LbField fTmp(1, D2Q9::nQ, grid.nElements()); // Do not need to be a LbField (as f), since its just a memory holder, that is immediately swaped after propagation.
    ScalarField rho(1, grid.nElements()); // Dette kan fikses for 'single rho fields' slik at man slipper å skrive rho(0, pos)
    VectorField vel(1, D2Q9::nD, grid.nElements()); // Bør vel-field vite at det er 2D. Dette kan nok gjøre at ting blir gjort raskere

    HalfWayBounceBack halfWayBounceBack(nBounceBackNodes, 4);
    Periodic periodic(nPeriodicNodes, 4);
    Bulk bulk(nBulkNodes);

    // SETUP BOUNDARY
    setupBoundary(halfWayBounceBack, periodic, grid, geo, nX, nY);
    // SETUP BULK
    setupBulk(bulk, grid, geo, nX, nY);

    
    
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
      //for (int y = 1; y < nY; y++ ) {
      //for (int x = 0; x < nX; x++) {
      //const int nodeNo = grid.element(x, y);
        for (int bulkNo = 0; bulkNo < bulk.nElements(); bulkNo++ ) {
            const int nodeNo = bulk.nodeNo(bulkNo);
            lbBase_t velNode[D2Q9::nD];
            lbBase_t rhoNode;
            // * Macrosopics
            // * * rho and vel

            //function that gives Macroscopic variables
            updateMacroscopicFields(&f(0,0, nodeNo), rhoNode, rho(0, nodeNo), velNode, &vel(0, 0, nodeNo), force);
            //D2Q9::qSum(&f(0,0, nodeNo), rhoNode);
            //D2Q9::qSumC(&f(0,0, nodeNo), velNode);



            //collision part start
            lbBase_t OmegaBGK_plus_f[D2Q9::nQ];
            lbBase_t deltaOmegaF[D2Q9::nQ];
            collision(tau_inv, &f(0, 0, nodeNo), velNode, rhoNode, force, factor_force, OmegaBGK_plus_f, deltaOmegaF);
            //collision part end

            //propagation part start

            for (int q = 0; q < D2Q9::nQ; q++) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = OmegaBGK_plus_f[q] + deltaOmegaF[q];//fTmp(0, q, grid.neighbor(q, nodeNo)) = fEql[q] + forcel[q];
            }
            //propagation part end
            // End collision and propagation */
            //}

        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);

        // BOUNDARY CONDITION Half way bounce back
        periodic.applyBoundaryCondition(f, grid);
        halfWayBounceBack.applyBoundaryCondition(f, grid);

    } // End iterations (LOOP TYPE 1)
    
    std::cout << std::setprecision(3) << vel(0, 0, grid.element(10, 10))  << std::endl;


    int** geoTest;
    makeGeometryX(5, 6, geoTest);
    printGeoScr(5, 6, geoTest);
    analyseGeometry<D2Q9>(5, 6, geoTest);
    printGeoScr(5, 6, geoTest);

    int** nodeLabels;
    makeNodeLabel(5, 6, nodeLabels);

    int nBulk = 0;
    nBulk = setBulkLabel(5, 6, geoTest, nodeLabels);
    int nNodes = 0;
    nNodes = setNonBulkLabel(nBulk, 5, 6, geoTest, nodeLabels);
    std::cout << "node numbers = " << nBulk << "   Node total = " << nNodes << std::endl;
    printGeoScr(5, 6, nodeLabels);

    Grid gridTmp(nNodes, D2Q9::nQ);
    setupGrid<D2Q9>(5, 6, nodeLabels, gridTmp);

    for (int n = 1; n <= nNodes; ++n) {
        std::cout << std::setw(3) << n <<": ";
        for (int q = 0; q < D2Q9::nQ; ++q)
        {
            std::cout << std::setw(3) << gridTmp.neighbor(q, n);
        }
        std::cout << std::endl;
    }
    deleteNodeLabel(5, 6, nodeLabels);
    deleteGeometryX(5, 6, geoTest);

    return 0;
}

