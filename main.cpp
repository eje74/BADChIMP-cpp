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

/*! LATTICE class
 *  Contains all information about the basis vectors, and operators for vectors and tensors
 *  like scalarproducs.
 */
class Lattice // Can we make static class types to increase speed. And add it with templates
{
public:
    Lattice(int nDirections, int nDimensions)  // Lattice constructor
    {
        nDirections_ = nDirections;                         /* Number of spatial dimensions */
        nDimensions_ = nDimensions;                         /* Number of basis velocities */
        w_ = new double [nDirections_];                     /* List of lattice weights */
        cDMajor_ = new int [nDirections_ * nDimensions_];   /*  Basis velocities on the form c[alpha][dim]
                                                             *  So that c = [c(alpha=0)_x|c(alpha=0)_y, ...]
                                                             */
        cQMajor_ = new int [nDirections_ * nDimensions_];   /* Basis velocities on the form c[dim][alpha].
                                                             * So that c = [c(alpha=0)_x|c(alpha=1)_x, ...]
                                                             */
        reverseStep_ = nDirections_ / 2;                    /* Number of bounce back pairs */
        nNonZeroDirections_ = 2 * reverseStep_;             /* Number of non-rest particles velocities */

        std::cout << "Constructed a Lattice object" << std::endl;
    }

    ~Lattice()  // Lattice destructor
    {
        delete [] w_;
        delete [] cDMajor_;
        delete [] cQMajor_;
        std::cout << "Destructed a Lattice object" << std::endl;
    }

    void setW(int qDirection, double value)  // Setter function for lattice weights
    {
        w_[qDirection] = value;
    }

    double getW(int qDirection)  // Getter function for lattice weights
    {
        return w_[qDirection];
    }

    int getC(int qDirection, int dimension)  // Getter function for basis vector components
    {
        return cDMajor_[nDimensions_*qDirection + dimension];
    }

    void setC(int qDirection, int dimension, int value)  // Setter function for basis vector components
    {
        cDMajor_[qDirection * nDimensions_ + dimension] = value;
        cQMajor_[dimension * nDirections_ + qDirection] = value;
    }

    int reverseDirection(int qDirection)  // Returns the reverse direction
    {
        /* Here we have assumed that the reverse directions are given as:
         *     alpha_reverse = alpha + <number of non zero directons>/2 ,
         * and that the rest particle velocity is that last element in the
         * basis vector.
         */
        return (qDirection + reverseStep_) % nNonZeroDirections_;
    }

    int getNQ()  // Getter for the number of lattice directions
    {
        return nDirections_;
    }

    int getND()  // Getter for the number of spatial dimensions
    {
        return nDimensions_;
    }

    double innerProductDMajor(int qDirection, double* vec) // Used for inner products with Cartesian vecotrs (c_{\alpha i} u_i)
    {
        double ret = 0;
        for (int d = 0; d < nDimensions_; d++) {
            ret += vec[d]*cDMajor_[qDirection*nDimensions_ + d];
        }
        return ret;
    }

    double innerProductQMajor(int dimension, double* vec) // Used to calculate first moments (\sum_\alpha c_{\alpha i} f_\alpha)
    {
        double ret = 0;
        for (int q = 0; q < nDirections_; q++) {
            ret += vec[q]*cQMajor_[dimension*nDirections_ + q];
        }
        return ret;
    }

    double innerProduct(double *leftVec, double *rightVec) // Standard Cartesion scalar product
    {
        double ret = 0;
        for (int d = 0; d < nDimensions_; d++) {
            ret += leftVec[d]*rightVec[d];
        }
        return ret;
    }

    // Different powers of the sound speed.
    double c2Inv_ = 3.0;
    double c2_ = 1.0 / c2Inv_;
    double c4Inv_ = 9.0;
    double c4_ = 1.0 / c4Inv_;

private:
    int nDirections_;  // Size of velocity basis
    int nDimensions_;  // Number of spacial dimensions
    double *w_; // Lattice weights
    int *cDMajor_;  // Velocity basis on the form c_\alpha,i = c[nD*alpha + i]
    int *cQMajor_;  // Velocity basis on the form c_\alpha,i = c[nQ*i + alpha]
    int reverseStep_;  // Used to reverse direction (assuming that \hat{\alpha} = \alpha + reverseStep_
    int nNonZeroDirections_;  // = nQ_ if no reste particle, = (nQ_ - 1) if rest particle
};


// Hva skal Grid inneholder?
// - Nabonoder i hver gridretning
// - Oversikt over bulknoder
// - Oversikt over 'boundary nodes'

class GridRegular  // Use ghost nodes. Should be a child a master Grid class in the finished code
{
public:
    GridRegular(int nX, int nY, Lattice *lattice)  // Input can also be one or many files ...
    {
        nX_ = nX + 2;
        nY_ = nY + 2;
        lattice_ = lattice;
        nElements_ = nX_ * nY_;
        std::cout << "Constructed a Grid object" << std::endl;
    }

    ~GridRegular()
    {
        std::cout << "Destructed a Grid object" << std::endl;
    }    

    int getnElements()
    {
        return nElements_;
    }


    void position(int &xNo, int &yNo, int elementNo) // Returns the position of the element
    {
        xNo = (elementNo % nX_) - 1;
        yNo = (elementNo / nX_) - 1;
    }

    int element(int xNo, int yNo) // The user should not need to care about the ghost nodes
    {
        return (xNo + 1) + (yNo + 1) * nX_;
    }

    int neighbor(int qDirection, int elementNo) // This should be generic to all grids
    {
        return elementNo + (*lattice_).getC(qDirection, 0) + nX_ * (*lattice_).getC(qDirection, 1);
    }

    int periodicNeighbor(int qDirection, int elementNo) /* elementNo = (xNo + 1) + nX_ * (yNo + 1) */
    {
        int xNo, yNo;
        position(xNo, yNo, elementNo);
        int xNeigNo = xNo + (*lattice_).getC(qDirection, 0);
        int yNeigNo = yNo + (*lattice_).getC(qDirection, 1);
        xNeigNo %= nX_ - 2;
        xNeigNo = (xNeigNo < 0) ? (xNeigNo + nX_ - 2) : (xNeigNo);
        yNeigNo %= nY_ - 2;
        yNeigNo = (yNeigNo < 0) ? (yNeigNo + nY_ - 2) : (yNeigNo);
        return neighbor((*lattice_).reverseDirection(qDirection), element(xNeigNo, yNeigNo));
    }

private:
    int nX_, nY_;
    int nElements_;
    Lattice *lattice_;
};


// template <typename double>
class Field
{
public:
    Field(int elementSize, int nElements)
    {
        elementSize_ = elementSize;
        nElements_ = nElements;
        data_ = new double [elementSize_ * nElements_];
        std::cout << "Constructed a Field object" << std::endl;
    }

    ~Field()
    {
        delete [] data_;
        std::cout << "Destructed a Field object" << std::endl;
    }

    double& operator () (int posInElement, int elementNo)
    {
        return data_[elementNo * elementSize_ + posInElement];  // Here we have choosen a collision based data structure
    }

    double* operator () (int elementNo) // Returns a pointer to the first position in elementNo
    {
        return &data_[elementNo * elementSize_];
    }
    void swapData(Field& field)
    {
        double* dataTmp;
        dataTmp = field.data_;
        field.data_ = data_;
        data_ = dataTmp;
    }

private:
    int elementSize_;
    int nElements_;
    double * data_;
};


class LbField: public Field // public Field<double>
{
public:
    LbField(int elementSize, int nElements, Lattice* lattice) : Field(elementSize, nElements)
    {
        lattice_ = lattice;
        std::cout << "Constructed a LbField object" << std::endl;
    }

    void zerothMoment(int elementNo, double* returnScalar)
    {
        (*returnScalar) = 0.0;
        for (int q = 0; q < (*lattice_).getNQ(); q++) {
            (*returnScalar) += (*this)(q, elementNo);
        }
    }

    void firstMoment(int elementNo, double* returnVector)
    {
        for (int d = 0; d < (*lattice_).getND(); d++)
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

    void applyBoundaryCondition(LbField& field, GridRegular& grid, Lattice& lattice)
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

    void applyBoundaryCondition(LbField& field, GridRegular& grid, Lattice& lattice)
    {
        // Need to pull all unknow values from ghost nodes. Those are the beta and delta links
        for (int n = 0; n < nBoundaryNodes_; ++n ) {
            for (int a = betaListBegin_[n]; a < betaListEnd_[n]; ++a) {
                int direction = linkList_[a + nLinkPairs_ * n];
                int reverseDirection = lattice.reverseDirection(direction);
                int node = boundaryNode_[n];

/*                grid.position(xNo, yNo, node);
                std::cout << "(" << xNo << "," << yNo << ")  ";
                grid.position(xNo, yNo, pushNode);
                std::cout << "(" << xNo << "," << yNo << ")  ";
                std::cout << " beta = " << direction << std::endl; */

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

    return lattice.getW(qDirection) * rho * ( 1.0 + lattice.c2Inv_ * cu + 0.5*lattice.c4Inv_ * (cu*cu - lattice.c2_*uu) );
}


bool insideDomain(int xNo, int yNo, int nX, int nY)
{
    return ( (xNo >= 0 ) && (xNo < nX) && (yNo >= 0) && (yNo < nY)   );
}


void makeGeometry(Field& geo, GridRegular& grid, int nX, int nY)
{

    // Periodic
    for (int y = -1; y < nY+1; ++y) {
        for (int x = -1; x < nX+1; ++x)
            geo(0, grid.element(x, y)) = 2;
    }
    // Fluid
    for (int y = 0; y < nY; ++y)
        for (int x = 0; x < nX; ++x)
            geo(0, grid.element(x, y)) = 0.0;

    // Solid
    for (int x = 0; x < nX; ++x) {
        geo(0, grid.element(x, 0)) = 1.0;
    }

    // Fill ghost nodes with solid
    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            if ( (y==0) || (x==0) || (y == nY-1) || (x = nX-1) ) {
                if (geo(0, grid.element(x, y)) == 1.0) {
                    for (int q = 0; q < 8; ++q) {
                        int neigNode = grid.neighbor(q, grid.element(x, y));
                        int xNeig, yNeig;
                        grid.position(xNeig, yNeig, neigNode);
                        if (!insideDomain(xNeig, yNeig, nX, nY)) {
                            geo(0, grid.periodicNeighbor(q, grid.element(x, y))) = 1.0;
                        }
                    }
                }
            }
        }
    }
}

bool containsLinkType(double linkType, int nodeNo, Field& geo, GridRegular& grid, Lattice& lattice)
{
    for (int q = 0; q < lattice.getNQ() - 1; ++q) {
        if (geo(0, grid.neighbor(q, nodeNo)) == linkType)
                return true;
    }
    return false;
}


/* Find number of periodic nodes and bounce back */
void numberOfBoundaryNodes(int& nBounceBackNodes, int& nPeriodicNodes, Field& geo, GridRegular& grid, Lattice& lattice, int nX, int nY)
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


void addBoundaryNode(int node, double geoType, Field& geo, Boundary& boundary, GridRegular& grid, Lattice& lattice)
{
    int nBeta=0, nGamma=0, nDelta=0;
    int betaList[4], gammaList[4], deltaList[4];
    bool linkPresent, reverseLinkPresent;

    if (geo(0, node) != 0.0)
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


void setupBoundary(HalfWayBounceBack& bounceBackBoundary, Periodic& periodicBoundary, GridRegular& grid,  Field& geo, Lattice& lattice, int nX, int nY)
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
    nIterations = 100;
    nX = 250; nY = 101;
    nQ = 9;
    nD = 2;
    tau = 1.0;


    // Simulation
    Lattice lattice(9, 2);
    setd2q9(lattice);    

    GridRegular grid(nX, nY, &lattice);
    Field geo(1, grid.getnElements());
    makeGeometry(geo, grid, nX, nY);
    numberOfBoundaryNodes(nBounceBackNodes, nPeriodicNodes, geo, grid, lattice, nX, nY);

    LbField f(nQ, grid.getnElements(), &lattice); // Bør f-field vite at det er nQ. Dette kan nok gjøre at ting blir gjort raskere
    Field fTmp(nQ, grid.getnElements()); // Do not need to be a LbField (as f), since its just a memory holder, that is immediately swaped after propagation.
    Field rho(1, grid.getnElements()); // Dette kan fikses for 'single rho fields' slik at man slipper å skrive rho(0, pos)
    Field vel(nD, grid.getnElements()); // Bør vel-field vite at det er 2D. Dette kan nok gjøre at ting blir gjort raskere
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
                for (int d = 0; d < lattice.getND(); ++d) {
                    vel(d, nodeNo) += 0.5 * force[d];
                    vel(d, nodeNo) /= rho(0, nodeNo);
                }

                // * Collision and propagation:
                for (int q = 0; q < nQ; q++) {  // Collision should provide the right hand side must be
                                                // collision(q, f(q, n), rho(n), u(n) \\ should be an array) ?
                    fTmp(q, grid.neighbor(q, nodeNo)) = f(q, nodeNo) - (1.0 / tau) * ( f(q, nodeNo) - LbEquilibirum(q, rho(0, nodeNo), vel(nodeNo), lattice) );
                    fTmp(q, grid.neighbor(q, nodeNo)) += lattice.getW(q) * (1 - 0.5 / tau) * (lattice.c2Inv_ * lattice.innerProductDMajor(q, force)
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

