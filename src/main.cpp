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

#include <sstream>
#include <iostream>
#include <fstream>
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

void printSpeed(int iterNo, int nX, int nY, VectorField &vel, int** &geo, int** &lable)
{
    std::stringstream tmp_string;

    tmp_string << "vel_" << std::setw(8) << std::setfill('0') << iterNo << ".dat";
    std::string filename = tmp_string.str();
    std::cout << filename << std::endl;

    std::ofstream datafile;
    datafile.open(filename);

    for (int y = 0; y < nY; ++y)
    {
        for (int x = 0; x < nX; ++x)
        {
            if (geo[y][x] < 3)
                datafile << vel(0, 0, lable[y][x]) << "(" << geo[y][x] << ")" << " ";
        }
        datafile << std::endl;
    }

    datafile.close();
}

void printGeo(int nX, int nY, int** &geo)
{
    for (int y = 0; y < nY; ++y)
    {
        for (int x = 0; x < nX; ++x)
        {
            std::cout << geo[y][x] << " ";
        }
        std::cout  << std::endl;
    }
}

void printLabel(int nX, int nY, int** &label)
{
    for (int y = 0; y < nY; ++y)
    {
        for (int x = 0; x < nX; ++x)
        {
            std::cout << std::setw(3) <<label[y][x] << " ";
        }
        std::cout  << std::endl;
    }
}

/*--------------------END OF FUNCTION DEFINITIONS--------------------*/


int main()
{
    std::cout << "Begin test";
    std::cout << std::endl;


    // INPUT DATA
    int nIterations;
    int nX, nY;
    lbBase_t tau, tau_inv;

    lbBase_t force[2] = {1.0e-8, 0.0};
    lbBase_t factor_force;
    nIterations = 10000;
    nX = 250; nY = 101;

    tau = 0.7;
    tau_inv = 1.0 / tau;
    factor_force = (1 - 0.5 / tau);

    // SETUP GEOMETRY
    int ** geo;
    newGeometry(nX, nY, geo); // LBgeometry
    inputGeometry(nX, nY, geo); // LBgeometry
    analyseGeometry<D2Q9>(nX, nY, geo); // LBgeometry

    int ** labels;
    newNodeLabel(nX, nY, labels); // LBgeometry

    int nBulkNodes = 0;
    nBulkNodes = setBulkLabel(nX, nY, geo, labels); // LBgeometry
    int nNodes = 0;
    nNodes = setNonBulkLabel(nBulkNodes, nX, nY, geo, labels); // LBgeometry
    std::cout << "NUMBER OF NODES = " << nNodes << std::endl;

    // ADD 0 AS DEFUALT NODE
    nNodes += 1;

    // USE NODENO 0 AS DUMMY NODE.
    // Add one extra storage for the default dummy node
    nNodes += 1;

    // SETUP GRID
    std::cout << "NUMBER OF NODES = " << nNodes << std::endl;
    Grid<D2Q9> grid(nNodes); // object declaration: neigList_ stores node numbers of neighbors;  pos_ stores Cartesian coordinates of given node
    setupGrid(nX, nY, labels, grid); // LBgeometry, maybe move to LBgrid?

    // SETUP BULK
    std::cout << "NUMBER OF Bulk NODES = " << nBulkNodes << std::endl;
    Bulk bulk(nBulkNodes); // object declaration: nBulkNodes_ stores number of Bulk nodes; bulkNode_ stores node numbers of all bulk nodes
    setupBulk(nX, nY, geo, labels, bulk); // LBgeometry, maybe move to LBbulk?

    // SETUP BOUNDARY
    HalfWayBounceBack<D2Q9> boundary( nBoundaryNodes(1, nX, nY, geo) ); // object declaration: 
    setupBoundary(1, nX, nY, geo, labels, grid, boundary); // LBgeometry

    // SETUP FIELDS
    LbField f(1, D2Q9::nQ, nNodes); // Bør f-field vite at det er nQ. Dette kan nok gjøre at ting blir gjort raskere
    LbField fTmp(1, D2Q9::nQ, nNodes); // Do not need to be a LbField (as f), since its just a memory holder, that is immediately swaped after propagation.
    ScalarField rho(1, nNodes); // Dette kan fikses for 'single rho fields' slik at man slipper å skrive rho(0, pos)
    VectorField vel(1, D2Q9::nD, nNodes); // Bør vel-field vite at det er 2D. Dette kan nok gjøre at ting blir gjort raskere

    // INIT FILEDS
    lbBase_t velTmp[D2Q9::nD] = {0.0, 0.0};
    for (int n = 0; n < bulk.nElements(); n++ ) {
        const int nodeNo = bulk.nodeNo(n);
        lbBase_t cu[D2Q9::nQ], uu;
        D2Q9::cDotAll(velTmp, cu);
        uu = D2Q9::dot(velTmp, velTmp);
        for (int q = 0; q < D2Q9::nQ; q++) {
	  f(0, q, nodeNo) = D2Q9::w[q] * (1.0 + D2Q9::c2Inv*cu[q] + D2Q9::c4Inv0_5*(cu[q]*cu[q] - D2Q9::c2*uu)); // f_eq
        }
    }

    // -----------------MAIN LOOP------------------
    for (int i = 0; i < nIterations; i++) {
        for (int bulkNo = 0; bulkNo < bulk.nElements(); bulkNo++ ) {
            // Find current node number
            const int nodeNo = bulk.nodeNo(bulkNo);

            // UPDATE MACROSCOPIC VARAIBLES
            lbBase_t velNode[D2Q9::nD];
            lbBase_t rhoNode;
            updateMacroscopicFields(&f(0,0, nodeNo), rhoNode, rho(0, nodeNo), velNode, &vel(0, 0, nodeNo), force);


            // COLLISION
            lbBase_t OmegaBGK_plus_f[D2Q9::nQ];
            lbBase_t deltaOmegaF[D2Q9::nQ];
            collision(tau_inv, &f(0, 0, nodeNo), velNode, rhoNode, force, factor_force, OmegaBGK_plus_f, deltaOmegaF);

            // PROPAGATION
            for (int q = 0; q < D2Q9::nQ; q++) {  // Collision should provide the right hand side must be
                fTmp(0, q,  grid.neighbor(q, nodeNo)) = OmegaBGK_plus_f[q] + deltaOmegaF[q];//fTmp(0, q, grid.neighbor(q, nodeNo)) = fEql[q] + forcel[q];
            }
        } // End nodes

        // Swap data_ from fTmp to f;
        f.swapData(fTmp);

        // BOUNDARY CONDITIONS
        boundary.apply(0, f, grid);

    } // End iterations (LOOP TYPE 1)
    // -----------------END MAIN LOOP------------------
    
    std::cout << std::setprecision(3) << vel(0, 0, labels[10][10])  << std::endl;


   // CLEANUP
    deleteNodeLabel(nX, nY, labels);
    deleteGeometry(nX, nY, geo);

    return 0;
}

