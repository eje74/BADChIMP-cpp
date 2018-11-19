#include "LBboundary.h"


Boundary::Boundary(int nBoundaryNodes, int nLinkPairs)
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

Boundary::~Boundary()
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

void Boundary::addBoundaryNode( int nodeNo,
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



HalfWayBounceBack::HalfWayBounceBack(int nBoundaryNodes, int nLinkPairs): Boundary(nBoundaryNodes, nLinkPairs)
{
}

Periodic::Periodic(int nBoundaryNodes, int nLinkPairs): Boundary(nBoundaryNodes, nLinkPairs)
{
}



