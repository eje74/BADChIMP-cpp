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
//#include "latticeboltzmann.h"


template <typename BASETYPE>
class Field
{
public:
    Field(int elementSize, int nElements)
    {
        elementSize_ = elementSize;
        nElements_ = nElements;
        data_ = new BASETYPE [elementSize_ * nElements_];
        std::cout << "Constructed a Field object" << std::endl;
    }

    ~Field()
    {
        delete [] data_;
        std::cout << "Destructed a Field object" << std::endl;
    }

    BASETYPE& operator () (int posInElement, int elementNo)
    {
        return data_[elementNo * elementSize_ + posInElement];  // Here we have choosen a collision based data structure
    }

    BASETYPE* operator () (int elementNo) // Returns a pointer to the first position in elementNo
    {
        return &data_[elementNo * elementSize_];
    }
    void swapData(Field& field)
    {
        BASETYPE* dataTmp;
        dataTmp = field.data_;
        field.data_ = data_;
        data_ = dataTmp;
    }

private:
    int elementSize_;
    int nElements_;
    BASETYPE* data_;
};

template <typename BASETYPE>
class LbField: public Field<BASETYPE> // public Field<double>
{
public:
    LbField(int elementSize, int nElements, Lattice* lattice) : Field<BASETYPE>(elementSize, nElements)
    {
        lattice_ = lattice;
        std::cout << "Constructed a LbField object" << std::endl;
    }

    void zerothMoment(int elementNo, double* returnScalar)
    {
        (*returnScalar) = 0.0;
        for (int q = 0; q < (*lattice_).nQ(); q++) {
            (*returnScalar) += (*this)(q, elementNo);
        }
    }

    void firstMoment(int elementNo, double* returnVector)
    {
        for (int d = 0; d < (*lattice_).nD(); d++)
            returnVector[d] = (*lattice_).innerProductQMajor(d, (*this)(elementNo));
    }

private:
    Lattice *lattice_;
};


// Different linke types
// beta  - betaHat  : one link crosses boundary. beta: unknow, and betaHat known.
// gamma - gammaHat : no links cross boundary. both known
// delta - deltaHat : two links cross boundary. both unknown
// example ||1,6(beta)     |3(gamma)     |4(delta)     ||
//         || 5,2(betaHat) | 7(gammaHat) | 0 (deltaHat)||
// linkList_ = [1, 6, 3, 4]
// nTypes_   = [2, 1, 1]
//
// Boundary. listOfAllWallBoundaryNodes = [2, 6, 8, 100]
// Boundary. linkList [1,4,5,6,   2,4,7,1,   2,4,7,1,   2,4,7,1,]
//
class Boundary // Ta høyde for beveglige grenser
{
public:
    Boundary(int nBoundaryNodes, int nLinkPairs)
    {
        nBoundaryNodes_ = nBoundaryNodes;
        nLinkPairs_ = nLinkPairs;
        restDirection_ = 2 * nLinkPairs_ + 1;
        nAddedNodes_ = 0;
        boundaryNode_ = new int [nBoundaryNodes_];
        linkList_ = new int [nBoundaryNodes_ * nLinkPairs_];
        betaListBegin_ = new int [nBoundaryNodes_];
        betaListEnd_ = new int [nBoundaryNodes_];
        gammaListBegin_ = new int [nBoundaryNodes_];
        gammaListEnd_ = new int [nBoundaryNodes_];
        deltaListBegin_ = new int [nBoundaryNodes_];
        deltaListEnd_ = new int [nBoundaryNodes_];
    }

    ~Boundary()

    {
        delete [] boundaryNode_;
        delete [] linkList_;
        delete [] betaListBegin_;
        delete [] betaListEnd_;
        delete [] gammaListBegin_;
        delete [] gammaListEnd_;
        delete [] deltaListBegin_;
        delete [] deltaListEnd_;
    }

    void addBoundaryNode( int nodeNo,
                          int nBetaDirs,
                          int* betaDirList,
                          int nGammaLinks,
                          int* gammaLinkList,
                          int nDeltaLinks,
                          int* deltaLinkList )
    {
        if (nAddedNodes_ < nBoundaryNodes_) {
            boundaryNode_[nAddedNodes_] = nodeNo;
        }
        else {
            std::cout << "ERROR: Added to many boundary nodes!" << std::endl;
            exit(1);
        }

        // Add the link information
        int linkListCounter = 0;
        // -- Beta dir
        betaListBegin_[nAddedNodes_] = 0;
        betaListEnd_[nAddedNodes_] = nBetaDirs;
        for (int q = 0; q < nBetaDirs; q++) {
            linkList_[linkListCounter + nLinkPairs_ * nAddedNodes_] = betaDirList[q];
            linkListCounter += 1;
        }
        // -- Gamma dir
        gammaListBegin_[nAddedNodes_] = betaListEnd_[nAddedNodes_];
        gammaListEnd_[nAddedNodes_] = nGammaLinks +  gammaListBegin_[nAddedNodes_];
        for (int q = 0; q < nGammaLinks; q++) {
            linkList_[linkListCounter + nLinkPairs_ * nAddedNodes_] = gammaLinkList[q];
            linkListCounter += 1;
        }
        // -- Delta dir
        deltaListBegin_[nAddedNodes_] = gammaListEnd_[nAddedNodes_];
        deltaListEnd_[nAddedNodes_] = nDeltaLinks + deltaListBegin_[nAddedNodes_];
        for (int q = 0; q < nDeltaLinks; q++) {
            linkList_[linkListCounter + nLinkPairs_ * nAddedNodes_] = deltaLinkList[q];
            linkListCounter += 1;
        }

        if (linkListCounter != nLinkPairs_) {
            std::cout << "ERROR: number of link pair added were " << linkListCounter << std::endl;
            std::cout << "Must be the same as nLinkPairs = " << nLinkPairs_ << std::endl;
        }

        nAddedNodes_ += 1;
    }
    // void removeBoundaryNode()


protected:
    int nBoundaryNodes_;
    int nLinkPairs_;
    int restDirection_;
    int nAddedNodes_; // Counter for number of added boundary nodes
    int* boundaryNode_;
    int* linkList_;
    int* betaListBegin_;
    int* betaListEnd_;
    int* gammaListBegin_;
    int* gammaListEnd_;
    int* deltaListBegin_;
    int* deltaListEnd_;

//    int* nTypes_;
};



class HalfWayBounceBack : public Boundary
{
public:
    HalfWayBounceBack(int nBoundaryNodes, int nLinkPairs): Boundary(nBoundaryNodes, nLinkPairs)
    {
    }

    template<typename BASETYPE>
    void applyBoundaryCondition(LbField<BASETYPE>& field, GridRegular& grid, Lattice& lattice)
    {        
        // Here we will go through all unknown directions. That is, beta and delta links.
        for (int n = 0; n < nBoundaryNodes_; n++) {
            for (int a = betaListBegin_[n]; a < betaListEnd_[n]; a++) {
                int direction = linkList_[a + nLinkPairs_ * n];
                int reverseDirection = lattice.reverseDirection(direction);
                int node = boundaryNode_[n];
                field(direction, node) = field(reverseDirection, grid.neighbor(reverseDirection, node) );
            }
            for (int a = deltaListBegin_[n]; a < deltaListEnd_[n]; a++) { // Remeber to use bounce back for both link pair directions
                int direction = linkList_[a + nLinkPairs_ * n];
                int reverseDirection = lattice.reverseDirection(direction);
                int node = boundaryNode_[n];
                field(direction, node) = field(reverseDirection, grid.neighbor(reverseDirection, node) );
                field(reverseDirection, node) = field(direction, grid.neighbor(direction, node) );
            }
        }
    }
};



/*
 *  For LbFields:
 *   Need to copy outgoing values FROM ghost nodes TO accompaning ghost nodes for incomming values
 *
*/
class Periodic: public Boundary
{
public:
    Periodic(int nBoundaryNodes, int nLinkPairs): Boundary(nBoundaryNodes, nLinkPairs)
    {
    }

/*    void saveNodeDistributionValues(LbField& field, GridRegular& grid, Lattice& lattice)
    {
    }

    void loadNodeDistributionValues(LbField& field, GridRegular& grid, Lattice& lattice)
    {

    } */

    int getnBoundaryNodes()
    {
        return nBoundaryNodes_;
    }

    template<typename BASETYPE>
    void applyBoundaryCondition(LbField<BASETYPE>& field, GridRegular& grid, Lattice& lattice)
    {
        // Need to pull all unknow values from ghost nodes. Those are the beta and delta links
        for (int n = 0; n < nBoundaryNodes_; ++n ) {
            for (int a = betaListBegin_[n]; a < betaListEnd_[n]; ++a) {
                int direction = linkList_[a + nLinkPairs_ * n];
                int reverseDirection = lattice.reverseDirection(direction);
                int node = boundaryNode_[n];

                field(direction, node) = field(direction, grid.periodicNeighbor(reverseDirection, node));
            }
            for (int a = deltaListBegin_[n]; a < deltaListEnd_[n]; ++a) {
                int direction = linkList_[a + nLinkPairs_ * n];
                int reverseDirection = lattice.reverseDirection(direction);
                int node = boundaryNode_[n];

                field(direction, node) = field(direction, grid.periodicNeighbor(reverseDirection, node));
                field(reverseDirection, node) = field(reverseDirection, grid.periodicNeighbor(direction, node));
            }
        }
    }
};



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


double LbEquilibirum(int qDirection, double rho, double* vel, Lattice &lattice)
{
    double cu, uu;

    uu = lattice.innerProduct(vel, vel);  // Only needed to be calulated once for each element, could be set at an input variable
    cu = lattice.innerProductDMajor(qDirection, vel);

    return lattice.w(qDirection) * rho * ( 1.0 + lattice.c2Inv_ * cu + 0.5*lattice.c4Inv_ * (cu*cu - lattice.c2_*uu) );
}


bool insideDomain(int xNo, int yNo, int nX, int nY)
{
    return ( (xNo >= 0 ) && (xNo < nX) && (yNo >= 0) && (yNo < nY)   );
}

template<typename BASETYPE>
void makeGeometry(Field<BASETYPE>& geo, GridRegular& grid, int nX, int nY)
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
bool containsLinkType(int linkType, int nodeNo, Field<BASETYPE>& geo, GridRegular& grid, Lattice& lattice)
{
    for (int q = 0; q < lattice.nQ() - 1; ++q) {
        if (geo(0, grid.neighbor(q, nodeNo)) == linkType)
                return true;
    }
    return false;
}


template<typename BASETYPE>
/* Find number of periodic nodes and bounce back */
void numberOfBoundaryNodes(int& nBounceBackNodes, int& nPeriodicNodes, Field<BASETYPE>& geo, GridRegular& grid, Lattice& lattice, int nX, int nY)
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
void addBoundaryNode(int node, int geoType, Field<BASETYPE>& geo, Boundary& boundary, GridRegular& grid, Lattice& lattice)
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
void setupBoundary(HalfWayBounceBack& bounceBackBoundary, Periodic& periodicBoundary, GridRegular& grid,  Field<BASETYPE>& geo, Lattice& lattice, int nX, int nY)
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
    double tau;

    double force[2] = {1.0e-8, 0.0};
    nIterations = 1000;
    nX = 250; nY = 101;
    nQ = 9;
    nD = 2;
    tau = 1.0;


    // Simulation
    Lattice lattice(9, 2);
    setd2q9(lattice);    

    GridRegular grid(nX, nY, &lattice);
    Field<int> geo(1, grid.nElements());
    makeGeometry(geo, grid, nX, nY);
    numberOfBoundaryNodes(nBounceBackNodes, nPeriodicNodes, geo, grid, lattice, nX, nY);

    LbField<double> f(nQ, grid.nElements(), &lattice); // Bør f-field vite at det er nQ. Dette kan nok gjøre at ting blir gjort raskere
    Field<double> fTmp(nQ, grid.nElements()); // Do not need to be a LbField (as f), since its just a memory holder, that is immediately swaped after propagation.
    Field<double> rho(1, grid.nElements()); // Dette kan fikses for 'single rho fields' slik at man slipper å skrive rho(0, pos)
    Field<double> vel(nD, grid.nElements()); // Bør vel-field vite at det er 2D. Dette kan nok gjøre at ting blir gjort raskere
    HalfWayBounceBack halfWayBounceBack(nBounceBackNodes, 4);
    Periodic periodic(nPeriodicNodes, 4);

    // SETUP BOUNDARY
    setupBoundary(halfWayBounceBack, periodic, grid, geo, lattice, nX, nY);

    // INIT:
    double velTmp[2] = {0.0, 0.0};
    for (int y = 0; y < nY; y++ ) {
        for (int x = 0; x < nX; x++) {
            for (int q = 0; q < nQ; q++) {
                f(q, grid.element(x, y)) = LbEquilibirum(q, 1.0, velTmp, lattice);
            }
        }
    }
    // - LOOP TYPE 1:
    for (int i = 0; i < nIterations; i++) {
   //     std::cout << "@ iteration " << i << std::endl;
        for (int y = 1; y < nY; y++ ) {
            for (int x = 0; x < nX; x++) {
                int nodeNo = grid.element(x, y);
                // * Macrosopics
                // * * rho and vel
                f.zerothMoment(nodeNo, rho(nodeNo));
                f.firstMoment(nodeNo, vel(nodeNo));
                for (int d = 0; d < lattice.nD(); ++d) {
                    vel(d, nodeNo) += 0.5 * force[d];
                    vel(d, nodeNo) /= rho(0, nodeNo);
                }

                // * Collision and propagation:
                for (int q = 0; q < nQ; q++) {  // Collision should provide the right hand side must be
                                                // collision(q, f(q, n), rho(n), u(n) \\ should be an array) ?
                    fTmp(q, grid.neighbor(q, nodeNo)) = f(q, nodeNo) - (1.0 / tau) * ( f(q, nodeNo) - LbEquilibirum(q, rho(0, nodeNo), vel(nodeNo), lattice) );
                    fTmp(q, grid.neighbor(q, nodeNo)) += lattice.w(q) * (1 - 0.5 / tau) * (lattice.c2Inv_ * lattice.innerProductDMajor(q, force)
                                                                                              + lattice.c4Inv_ *
                                                                                              (
                                                                                                  lattice.innerProductDMajor(q, force) * lattice.innerProductDMajor(q, vel(nodeNo)) - lattice.c2_ * lattice.innerProduct(vel(nodeNo), force))
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

    std::cout << std::setprecision(3) << vel(0, grid.element(10, 10))  << std::endl;

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

