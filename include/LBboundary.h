#ifndef LBBOUNDARY_H
#define LBBOUNDARY_H


#include <iostream>
#include "LBglobal.h"
#include "LBfield.h"
#include "LBgrid.h"
#include "LBd2q9.h"
//#include "LBlattice.h"

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
    Boundary(int nBoundaryNodes, int nLinkPairs);
    ~Boundary();

    void addBoundaryNode( int nodeNo,
                          int nBetaDirs,
                          int* betaDirList,
                          int nGammaLinks,
                          int* gammaLinkList,
                          int nDeltaLinks,
                          int* deltaLinkList );

protected:
    int nBoundaryNodes_; // Number of boundary nodes
    int nLinkPairs_; // Number of linke pairs, a lattice constant
    int restDirection_; // Rest velocities direction.
    int nAddedNodes_; // Counter for number of added boundary nodes
    int* boundaryNode_;
    int* linkList_;
    int* betaListBegin_;
    int* betaListEnd_;
    int* gammaListBegin_;
    int* gammaListEnd_;
    int* deltaListBegin_;
    int* deltaListEnd_;
};



class HalfWayBounceBack : public Boundary
{
public:
    HalfWayBounceBack(int nBoundaryNodes, int nLinkPairs);

    template <typename DXQY>
    void applyBoundaryCondition(LbField& field, GridRegular<DXQY>& grid);

};


template <typename DXQY>
inline void HalfWayBounceBack::applyBoundaryCondition(LbField& field, GridRegular<DXQY>& grid)
{
    // Here we will go through all unknown directions. That is, beta and delta links.
    for (int n = 0; n < nBoundaryNodes_; n++) {
        for (int a = betaListBegin_[n]; a < betaListEnd_[n]; a++) {
            int direction = linkList_[a + nLinkPairs_ * n];
            int reverseDirection = DXQY::reverseDirection(direction);
            int node = boundaryNode_[n];
            field(0, direction, node) = field(0, reverseDirection, grid.neighbor(reverseDirection, node) );
        }
        for (int a = deltaListBegin_[n]; a < deltaListEnd_[n]; a++) { // Remeber to use bounce back for both link pair directions
            int direction = linkList_[a + nLinkPairs_ * n];
            int reverseDirection = DXQY::reverseDirection(direction);
            int node = boundaryNode_[n];
            field(0, direction, node) = field(0, reverseDirection, grid.neighbor(reverseDirection, node) );
            field(0, reverseDirection, node) = field(0, direction, grid.neighbor(direction, node) );
        }
    }
}


/*
 *  For LbFields:
 *   Need to copy outgoing values FROM ghost nodes TO accompaning ghost nodes for incomming values
 *
*/
class Periodic: public Boundary
{
public:
    Periodic(int nBoundaryNodes, int nLinkPairs);

/*    void saveNodeDistributionValues(LbField& field, GridRegular& grid, Lattice& lattice)
    {
    }

    void loadNodeDistributionValues(LbField& field, GridRegular& grid, Lattice& lattice)
    {

    } */

    int getnBoundaryNodes();

    template <typename DXQY>
    void applyBoundaryCondition(LbField& field, GridRegular<DXQY>& grid);
};


inline int Periodic::getnBoundaryNodes()
{
    return nBoundaryNodes_;
}

template <typename DXQY>
inline void Periodic::applyBoundaryCondition(LbField& field, GridRegular<DXQY>& grid)
{
    // Need to pull all unknow values from ghost nodes. Those are the beta and delta links
    for (int n = 0; n < nBoundaryNodes_; ++n ) {
        for (int a = betaListBegin_[n]; a < betaListEnd_[n]; ++a) {
            int direction = linkList_[a + nLinkPairs_ * n];
            int reverseDirection = DXQY::reverseDirection(direction);
            int node = boundaryNode_[n];

            field(0, direction, node) = field(0, direction, grid.periodicNeighbor(reverseDirection, node));
        }
        for (int a = deltaListBegin_[n]; a < deltaListEnd_[n]; ++a) {
            int direction = linkList_[a + nLinkPairs_ * n];
            int reverseDirection = DXQY::reverseDirection(direction);
            int node = boundaryNode_[n];

            field(0, direction, node) = field(0, direction, grid.periodicNeighbor(reverseDirection, node));
            field(0, reverseDirection, node) = field(0, reverseDirection, grid.periodicNeighbor(direction, node));
        }
    }
}



#endif // LBBOUNDARY_H
