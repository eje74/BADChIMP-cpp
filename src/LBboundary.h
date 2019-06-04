#ifndef LBBOUNDARY_H
#define LBBOUNDARY_H


#include <cassert> // Contains assert functions (see eq. addNode function)
#include <vector>
#include "LBglobal.h"
#include "LBlatticetypes.h"

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
class Boundary
{
public:
    Boundary(int nBoundaryNodes);
    //~Boundary();

    void addNode( int nodeNo,
                  int nBetaLinks, int* betaDirList,
                  int nGammaLinks, int* gammaLinkList,
                  int nDeltaLinks, int* deltaLinkList );

    inline std::vector<int> beta(const int bndNo) const;
    inline std::vector<int> gamma(const int bndNo) const;
    inline std::vector<int> delta(const int bndNo) const;
    inline int betaOld(const int dirNo, const int bndNo) const;
    inline int gammaOld(const int dirNo, const int bndNo) const;
    inline int deltaOld(const int dirNo, const int bndNo) const;
    inline int dirRev(const int dir) const;
    inline int nodeNo(const int bndNo) const;

    // Getter for number of boundary nodes (Should we call it size ? )
    inline int size() const {return nBoundaryNodes_;}
    inline int nBeta(const int bndNo) const {return nBeta_[bndNo];}
    inline int nGamma(const int bndNo) const {return nGamma_[bndNo];}
    inline int nDelta(const int bndNo) const {return nDelta_[bndNo];}

protected:
    int nBoundaryNodes_;  // Number of boundary nodes
    //int* boundaryNode_;  // List containing the node number (tag) of the boundary nodes
    //int* linkList_;  // List that contains the information about which directions are in the beta-, gamma-, delta-categories
    //int* nBeta_;  // List of the number of beta links for each boundary node
    //int* nGamma_;  // List of the number of gamma links for each boundary node
    //int* nDelta_;  // List of the number of delta links for each boundary node
    int nAddedNodes_; // Counter for number of added boundary nodes
    int nAddedLinks_; // Counter for the number of added links
    std::vector<int> nodes_;
    std::vector<int> boundaryNode_; // List containing the node number (tag) of the boundary nodes
    std::vector<int> linkList_; // List that contains the information about which directions are in the beta-, gamma-, delta-categories
    std::vector<int> nBeta_; // List of the number of beta links for each boundary node
    std::vector<int> nGamma_;  // List of the number of gamma links for each boundary node
    std::vector<int> nDelta_;  // List of the number of delta links for each boundary node
    
};

template <typename DXQY>
Boundary<DXQY>::Boundary(int nBoundaryNodes) : nBoundaryNodes_(nBoundaryNodes),
  nodes_(static_cast<std::size_t>(nBoundaryNodes)),
  boundaryNode_(static_cast<std::size_t>(nBoundaryNodes)),
  linkList_(nBoundaryNodes_ * DXQY::nDirPairs_),
  nBeta_(static_cast<std::size_t>(nBoundaryNodes)),
  nGamma_(static_cast<std::size_t>(nBoundaryNodes)),
  nDelta_(static_cast<std::size_t>(nBoundaryNodes))
/* Constructor for a Boundary object. Allocates memory for
 * all lists used by the object.
 *
 * nBoundaryNodes : number of boundary nodes
 */
{
  //nBoundaryNodes_ = nBoundaryNodes;
    nAddedNodes_ = 0;
    nAddedLinks_ = 0;
    //boundaryNode_ = new int [nBoundaryNodes_];
    //linkList_ = new int [nBoundaryNodes_ * DXQY::nDirPairs_];
    //nBeta_ = new int [nBoundaryNodes_];
    //nGamma_ = new int [nBoundaryNodes_];
    //nDelta_ = new int [nBoundaryNodes_];
}

/*
template <typename DXQY>
Boundary<DXQY>::~Boundary()
// Destructor for a Boundary object. Deletes all list created
// by the constructor.
/ /
{
    //removed delete [] boundaryNode_;
    //removed delete [] linklist_;
    delete [] nBeta_;
    delete [] nGamma_;
    delete [] nDelta_;
}
*/

template <typename DXQY>
void Boundary<DXQY>::addNode( int nodeNo,
                              int nBetaLinks,
                              int* betaLinkList,
                              int nGammaLinks,
                              int* gammaLinkList,
                              int nDeltaLinks,
                              int* deltaLinkList )
/* addNode adds all lattice information about a boundary node, and updates the
 * boundary-class object. A lattice pair is represented by one of its lattice
 * direction. For gamma and delta links either one will do. But for beta links
 * it is assumed that it is the unknown direction that is given as input.
 *
 * nodeNo : the node number (tag) to the boundary node.
 * nBetaLinks : the number of link pairs that are of beta type.
 * betaLinkList : list of the unknown directions the beta type link pair
 * nGammaLinks : the number of link pairs that are of gamma type.
 * gammaLinkList : list of one of the directions in the gamma type link pair
 * nDeltaLinks : the number of link pairs that are of delta type.
 * deltaLinkList : list of one of the directions in the delta type link pair
 *
 * The structure of the linkList_ is first to list all beta links for a given boundary node,
 * then all the gamma links and finally all delta links. This is then repeated for all
 * boundary nodes.
 */
{
    assert(nAddedNodes_ < nBoundaryNodes_); // Checks that it is room for a new node
    boundaryNode_[nAddedNodes_] = nodeNo;
    nodes_[nAddedNodes_] = nodeNo;

    assert( (nBetaLinks + nGammaLinks + nDeltaLinks) == DXQY::nDirPairs_); // Checks that it is room for the link information
    // Add beta links
    nBeta_[nAddedNodes_] = nBetaLinks;
    for (int n = 0; n < nBetaLinks; ++n) {
        linkList_[nAddedLinks_] = betaLinkList[n];
        nAddedLinks_ += 1;
    }
    // Add gamma links
    nGamma_[nAddedNodes_] = nGammaLinks;
    for (int n = 0; n < nGammaLinks; ++n) {
        linkList_[nAddedLinks_] = gammaLinkList[n];
        nAddedLinks_ += 1;
    }
    // Add delta links
    nDelta_[nAddedNodes_] = nDeltaLinks;
    for (int n = 0; n < nDeltaLinks; ++n) {
        linkList_[nAddedLinks_] = deltaLinkList[n];
        nAddedLinks_ += 1;
    }

    nAddedNodes_ += 1;  // Increase the local counter
}


template <typename DXQY>
inline std::vector<int> Boundary<DXQY>::beta(const int bndNo) const
/* beta returns a vector of the unknown direction for the beta links for the the
 *  given boundary nodes. That is, the directions pointing into the fluid.
 *
 * bndNo : boundary number: legal values are from 0 to nBoundaryNodes_
 */
{
    auto ptrBegin = linkList_.data() + bndNo * DXQY::nDirPairs_;
    return std::vector<int>(ptrBegin, ptrBegin + nBeta_[bndNo]);
}

template <typename DXQY>
inline std::vector<int> Boundary<DXQY>::gamma(const int bndNo) const
/* gamma returns a vector to stored gamma links. Gamma links contains
 *  no unknown directions.
 * NB remember to include bounce back directions in loops, as just one
 *  direction in each bounce back direction pair is returned.
 *
 * bndNo : boundary number: legal values are from 0 to nBoundaryNodes_
 */
{
    auto ptrBegin = linkList_.data() + bndNo * DXQY::nDirPairs_ + nBeta_[bndNo];
    return std::vector<int>(ptrBegin, ptrBegin + nGamma_[bndNo]);
}

template <typename DXQY>
inline std::vector<int> Boundary<DXQY>::delta(const int bndNo) const
/* delta returns a vector to direction in a delta link. Delta links contains
 *  only unknown directions.
 *  NB remember to include bounce back directions in loops, as just one
 *  direction in each bounce back direction pair is returned.
 *
 * bndNo : boundary number: legal values are from 0 to nBoundaryNodes_
 */
{
    auto ptrBegin = linkList_.data() + bndNo * DXQY::nDirPairs_ + nBeta_[bndNo] + nGamma_[bndNo];
    return std::vector<int>(ptrBegin, ptrBegin + nDelta_[bndNo]);
}


template <typename DXQY>
inline int Boundary<DXQY>::betaOld(const int dirNo, const int bndNo) const
/* beta returns the unknown direction in a beta link. That is the direction
 *   pointing into the fluid.
 *
 * dirNo : link number: legal values are from 0 to nBeta_
 * bndNo : boundary number: legal values are from 0 to nBoundaryNodes_
 */
{
    return linkList_[bndNo * DXQY::nDirPairs_ + dirNo];
}


template <typename DXQY>
inline int Boundary<DXQY>::gammaOld(const int dirNo, const int bndNo) const
/* gamma returns the stored direction in a gamma link. Gamma links contains
 *  no unknown directions.
 *
 * dirNo : link number: legal values are from 0 to nGamma_
 * bndNo : boundary number: legal values are from 0 to nBoundaryNodes_
 */
{
    return linkList_[bndNo * DXQY::nDirPairs_ + nBeta_[bndNo] + dirNo];
}

template <typename DXQY>
inline int Boundary<DXQY>::deltaOld(const int dirNo, const int bndNo) const
/* delta returns the stored direction in a delta link. Delta links contains
 *  only unknown directions.
 *
 * dirNo : link number: legal values are from 0 to nDelta_
 * bndNo : boundary number: legal values are from 0 to nBoundaryNodes_
 */
{
    return linkList_[bndNo * DXQY::nDirPairs_ + nBeta_[bndNo] + nGamma_[bndNo]+ dirNo];
}

template <typename DXQY>
inline int Boundary<DXQY>::dirRev(const int dir) const
/* dirRev returns the reverse direction of non-zero velocity dir.
 *
 * dir : lattice direction
 */
{
    return DXQY::reverseDirection(dir);
}

template <typename DXQY>
inline int Boundary<DXQY>::nodeNo(const int bndNo) const
/* nodeNo returns the node number (tag) of the given boundary node.
 *
 * bndNo : boundary node number
 */
{
    return boundaryNode_[bndNo];
}


#endif // LBBOUNDARY_H
