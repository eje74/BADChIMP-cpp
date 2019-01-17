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
        int node = this->boundaryNode_[n];

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
