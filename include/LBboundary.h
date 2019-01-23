#ifndef LBBOUNDARY_H
#define LBBOUNDARY_H


#include <iostream>
#include "LBglobal.h"
#include "LBfield.h"
#include "LBgrid.h"
#include "LBd2q9.h"
#include <assert.h>
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

/************************************************************
 * class BOUNDARY: super class used in boundary condtions.
 *
 * Classification of boundary nodes:
 *  * A node pair is defined as lattice direction and its reverse.
 *     That is, the speeds along a given non zero speed link direction
 *  * There are three types of node pairs beta, gamma and delta, which
 *     are charcterized according to if the values are known after
 *     streaming. This is, are values streamed from fluid (known) or
 *     solid (unknown) nodes.
 *  * Note: just one direction for each link is recorded. So to get all
 *     direction use the revDir of the recorded directions as well.
 *
 *  - beta          : Unknown
 *  - beta_revers   : Known
 *
 *  - gamma         : Known
 *  - gamma_reverse : Known
 *
 *  - delta         : Unknown
 *  - delta_reverse : Unknown
 *
 *  Example:
 *    // Going through all boundary nodes and write all directions
 *    for (int n = 0; n < nBoundaryNodes_; n++) {
 *       int node = nodeNo(n); // Get node number
 *       // Print all beta values
 *       std::cout << "beta =";
 *       for (int q = 0; q < nBeta_; ++q) {
 *          std::cout << " " << beta(q, n) << "(Unknown)";
 *          std::cout << " " << revDir( beta(q, n) ) << "(Known)";
 *       }
 *       std::cout << std::endl;
 *       // Print all gamma values
 *       std::cout << "gamma =";
 *       for (int q = 0; q < nGamma_; ++q) {
 *          std::cout << " " << gamma(q, n) << "(Known)";
 *          std::cout << " " << revDir( gamma(q, n) ) << "(Known)";
 *       }
 *       std::cout << std::endl;
 *       // Print all delta values
 *       std::cout << "delta =";
 *       for (int q = 0; q < nDelta_; ++q) {
 *          std::cout << " " << delta(q, n) << "(Unknown)";
 *          std::cout << " " << revDir( delta(q, n) ) << "(Unknown)";
 *       }
 *       std::cout << std::endl;
 ************************************************************/
template <typename DXQY>
class Boundary // Ta høyde for beveglige grenser
{
public:
    Boundary(int nBoundaryNodes);
    ~Boundary();

    void addNode( int nodeNo,
                  int nBetaLinks,
                  int* betaDirList,
                  int nGammaLinks,
                  int* gammaLinkList,
                  int nDeltaLinks,
                  int* deltaLinkList );
    int beta(const int dirNo, const int bndNo) const;
    int gamma(const int dirNo, const int bndNo) const;
    int delta(const int dirNo, const int bndNo) const;
    int dirRev(const int dir) const;
    int nodeNo(const int bndNo) const;

protected:
    int nBoundaryNodes_; // Number of boundary nodes
    int nAddedNodes_; // Counter for number of added boundary nodes
    int nAddedLinks_; // Connter for the number of added links
    int* boundaryNode_;
    int* linkList_;
    int* nBeta_;
    int* nGamma_;
    int* nDelta_;
};

template <typename DXQY>
Boundary<DXQY>::Boundary(int nBoundaryNodes)
{
    nBoundaryNodes_ = nBoundaryNodes;
    nAddedNodes_ = 0;
    nAddedLinks_ = 0;
    boundaryNode_ = new int [nBoundaryNodes_];
    linkList_ = new int [nBoundaryNodes_ * DXQY::nDirPairs_];
    nBeta_ = new int [nBoundaryNodes_];
    nGamma_ = new int [nBoundaryNodes_];
    nDelta_ = new int [nBoundaryNodes_];
}

template <typename DXQY>
Boundary<DXQY>::~Boundary()
{
    std::cout << "BOUNDARY destructor begin " << std::endl;
    delete [] boundaryNode_;
    delete [] linkList_;
    delete [] nBeta_;
    delete [] nGamma_;
    delete [] nDelta_;
    std::cout << "BOUNDARY destructor end " << std::endl;
}

template <typename DXQY>
void Boundary<DXQY>::addNode( int nodeNo,
                              int nBetaLinks,
                              int* betaDirList,
                              int nGammaLinks,
                              int* gammaLinkList,
                              int nDeltaLinks,
                              int* deltaLinkList )
{
    assert(nAddedNodes_ < nBoundaryNodes_);
    boundaryNode_[nAddedNodes_] = nodeNo;

    assert( (nBetaLinks + nGammaLinks + nDeltaLinks) == DXQY::nDirPairs_);
    nBeta_[nAddedNodes_] = nBetaLinks;
    for (int n = 0; n < nBetaLinks; ++n) {
        linkList_[nAddedLinks_] = betaDirList[n];
        nAddedLinks_ += 1;
    }
    nGamma_[nAddedNodes_] = nGammaLinks;
    for (int n = 0; n < nGammaLinks; ++n) {
        linkList_[nAddedLinks_] = gammaLinkList[n];
        nAddedLinks_ += 1;
    }
    nDelta_[nAddedNodes_] = nDeltaLinks;
    for (int n = 0; n < nDeltaLinks; ++n) {
        linkList_[nAddedLinks_] = deltaLinkList[n];
        nAddedLinks_ += 1;
    }

    nAddedNodes_ += 1;
}

template <typename DXQY>
inline int Boundary<DXQY>::beta(const int dirNo, const int bndNo) const
{
    return linkList_[bndNo * DXQY::nDirPairs_ + dirNo];
}


template <typename DXQY>
inline int Boundary<DXQY>::gamma(const int dirNo, const int bndNo) const
{
    return linkList_[bndNo * DXQY::nDirPairs_ + nBeta_[bndNo] + dirNo];
}

template <typename DXQY>
inline int Boundary<DXQY>::delta(const int dirNo, const int bndNo) const
{
    return linkList_[bndNo * DXQY::nDirPairs_ + nBeta_[bndNo] + nGamma_[bndNo]+ dirNo];
}

template <typename DXQY>
inline int Boundary<DXQY>::dirRev(const int dir) const
{
    return DXQY::reverseDirection(dir);
}

template <typename DXQY>
inline int Boundary<DXQY>::nodeNo(const int bndNo) const
{
    return boundaryNode_[bndNo];
}

template <typename DXQY>
class HalfWayBounceBack : public Boundary<DXQY>
{
public:
    HalfWayBounceBack(int nBoundaryNodes) : Boundary<DXQY>(nBoundaryNodes) {}
   // ~HalfWayBounceBack() {this->~Boundary<DXQY>();}
    void apply(const int fieldNo, LbField &f, const Grid<DXQY> &grid) const;
};

template <typename DXQY>
inline void HalfWayBounceBack<DXQY>::apply(const int fieldNo, LbField &f, const Grid<DXQY> &grid) const
{
    for (int n = 0; n < this->nBoundaryNodes_; ++n) {
        int node = this->nodeNo(n);

        // Bounce back for the beta directions (beta unknow)
        for (int q = 0; q < this->nBeta_[n]; ++q) {
            int beta = this->beta(q, n);
            int beta_rev = this->dirRev(beta);
            f(fieldNo, beta, node) = f(fieldNo, beta_rev, grid.neighbor(beta_rev, node));
        }

        // Bounce back for the delta directions (delta and delta.rev unknown)
        for (int q  = 0; q < this->nDelta_[n]; ++q) {
            int delta = this->delta(q, n);
            int delta_rev = this->dirRev(delta);

            f(fieldNo, delta, node) = f(fieldNo, delta_rev, grid.neighbor(delta_rev, node));
            f(fieldNo, delta_rev, node) = f(fieldNo, delta, grid.neighbor(delta, node));
        }
    }
}


#endif // LBBOUNDARY_H
