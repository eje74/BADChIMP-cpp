#ifndef LBBOUNDARYBASIC_H
#define LBBOUNDARYBASIC_H


#include <cassert> // Contains assert functions (see eq. addNode function)
#include <vector>
#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBgrid.h"

/************************************************************
 * class BOUNDARY: super class used in boundary condtions.
 *
 * Note on node classification:
 *  * a solid node is a term, here used, as a node that do not stream
 *    a value. Hence, if a neighbor node is a solid the distribution
 *    propageted from that node will be treated as unknown.
 *  * A solid node will not alway be a solid part of the geometry.
 *    For instance will the ghost node neighbors of inlet and outlet
 *    bundaries be treated as solids.
 *  * in the Node class. isSolid will be assumed to mean a node that
 *    do not propegate a distribution. This will be used in the Boundary
 *    constructor.
 *
 * Classification of boundary nodes:
 * * A direction will classified as  known or unknown
 *    - known : means that a distribution is propegated from a
 *              neighboring node
 *    - unknown : means that no distribution is available for
 *                 propagation in that direction
 ************************************************************/

template <typename DXQY>
class BoundaryBasic
{
public:
    BoundaryBasic(const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid);

    // Node Numbers
    inline int size() const {
        return nBoundaryNodes_;
    }
    inline int nodeNo(const int bndNo) const {
        return nodeNumbers_[bndNo];
    }

    // Know directions
    inline int nKnownDirs(const int bndNo) const {
        return nodeNumKnownDirs_[bndNo];
    }
    std::vector<int> knownDirs(const int bndNo) const;

    // Unknown directions
    inline int nUnknownDirs(const int bndNo) const {
        return DXQY::nQ - nodeNumKnownDirs_[bndNo];
    }
    std::vector<int> unknownDirs(const int bndNo) const;

    // Reverse directions
    inline int dirRev(const int dir) const;

    template<typename T>
    void addNode(const int bndNo, const int nodeNo, const int nKnownDir, const T& knownDirList);

protected:
    int nBoundaryNodes_;  // Number of boundary nodes
    std::vector<int> nodeNumKnownDirs_;
    std::vector<int> nodeDirections_;
    std::vector<int> nodeNumbers_; // List containing the node number (tag) of the boundary nodes

};



template <typename DXQY>
BoundaryBasic<DXQY>::BoundaryBasic(const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid) :
    nBoundaryNodes_(static_cast<int>(bndNodes.size())),
    nodeNumKnownDirs_(static_cast<std::size_t>(nBoundaryNodes_)),
    nodeNumbers_(static_cast<std::size_t>(nBoundaryNodes_)),
    nodeDirections_(nBoundaryNodes_ * DXQY::nQ)
{
    // Counter for boundary number
    int bndNo = 0;  
    
    // Loop over nodes in the boundary
    for (auto nodeNo: bndNodes) 
    {  
        // Number of known directions
        int numKnownDir = 0;
        // List of know directions
        std::array<int, DXQY::nQ> knowDirs;
        
        // loop over all neighbors
        for (int q = 0; q < DXQY::nQ; ++q)
        {
            // Neighbors node number
            int neigNo = grid.neighbor(q, nodeNo);
            
            // Add direction 
            if (nodes.isFluid(neigNo) )
            {
                knowDirs[numKnownDir++] = dirRev(q);
            }
        }
        
        // Add boundary node information
        addNode(bndNo, nodeNo, numKnownDir, knowDirs);
        
        // Update bondary number
        bndNo += 1;
    } // For all boundary nodes
}



template <typename DXQY>
std::vector<int> BoundaryBasic<DXQY>::knownDirs(const int bndNo) const
{
    // Position of first element in the list of know directions
    auto ptrBegin = nodeDirections_.begin() + DXQY::nQ * bndNo;

    return std::vector<int>(ptrBegin, ptrBegin + nKnownDirs(bndNo));
}



template <typename DXQY>
std::vector<int> BoundaryBasic<DXQY>::unknownDirs(const int bndNo) const
{
    // Position of first element in the list of unknow directions
    auto ptrBegin = nodeDirections_.begin() + DXQY::nQ * bndNo + nKnownDirs(bndNo);

    return std::vector<int>(ptrBegin, ptrBegin + nUnknownDirs(bndNo));
}



template <typename DXQY>
inline int BoundaryBasic<DXQY>::dirRev(const int dir) const
/* dirRev returns the reverse direction of non-zero velocity dir.
 *
 * dir : lattice direction
 */
{
    return DXQY::reverseDirection(dir);
}



template <typename DXQY>
template <typename T>
void BoundaryBasic<DXQY>::addNode(const int bndNo, const int nodeNo, const int nKnownDir, const T& knownDirList)
{
    // List of directions added
    std::vector<bool> addedDir(DXQY::nQ, false);

    // Set node number
    nodeNumbers_[bndNo] = nodeNo;

    // Set number of known directions
    nodeNumKnownDirs_[bndNo] = nKnownDir;

    // Add know directions
    int ind = 0;
    for (int q = 0; q < nKnownDir; ++q) {
        int dir = knownDirList[q];
        nodeDirections_[DXQY::nQ*bndNo + ind] = dir;
        addedDir[dir] = true;
        ind += 1;
    }

    // Add unknow directions
    for (int q = 0; q < DXQY::nQ; ++q) {
        if (!addedDir[q]) {
            nodeDirections_[DXQY::nQ*bndNo + ind] = q;
            ind += 1;
        }
    }
}


#endif // LBBOUNDARYBASIC_H
