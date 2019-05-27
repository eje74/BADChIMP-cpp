#ifndef LBNODES_H
#define LBNODES_H

#include <vector>
#include <numeric>
#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "Input.h"
#include "LBgrid.h"


/*
 * Nodes: -contains additional informatino about a node, eg. rank, node type etc.
 *        -class Grid will still contain the cartesion indices of the nodes and
 *         the list of node-numbers of the nodes in a node's neighborhood.
 *
 * Work in progress: moving node-information from class Grid to class Nodes
*/

template<typename DXQY>
class Nodes
{
public:
    Nodes(const int nNodes): nNodes_(nNodes), nodeType_(static_cast<std::size_t>(nNodes)) {}
    inline int size() const {return nNodes_;}
    void setup(MpiFile<DXQY> &mfs, MpiFile<DXQY> &rfs);
    inline int getType(const int nodeNo) const {return nodeType_[nodeNo];}
    inline int getRank(const int nodeNo) const {return nodeType_[nodeNo]-1;}
    static Nodes<DXQY> makeObject(MpiFile<DXQY> &mfs, MpiFile<DXQY> &rfs, const Grid<DXQY> &grid);

private:
    int nNodes_;
    std::vector<int> nodeType_; // Zero for solid and rank = value-1
};


template<typename DXQY>
void Nodes<DXQY>::setup(MpiFile<DXQY> &mfs, MpiFile<DXQY> &rfs)
{
    for (int pos=0; pos < static_cast<int>(mfs.size()); ++pos)
    {
        auto nodeNo = mfs.template getVal<int>(); // Gets node number. The template name is needed
        auto nodeType =  rfs.template getVal<int>(); // get the node type
        if ( mfs.insideDomain(pos) ) { // Update the grid object
            if (nodeNo > 0) { // Only do changes if it is a non-default node
                // Add node type
                nodeType_[nodeNo] = nodeType;
            }
        }
    }
    mfs.reset();
    rfs.reset();
}

template<typename DXQY>
Nodes<DXQY> Nodes<DXQY>::makeObject(MpiFile<DXQY> &mfs, MpiFile<DXQY> &rfs, const Grid<DXQY> &grid)
/* Makes a grid object using the information in the file created by our
 * python program, for each mpi-processor.
 *
 * mfs : local node number file
 * rfs : rank number file
 *
 * We assume that the file contains:
 *  1) Dimesions of the system (including the rim)
 *  2) The global Cartesian coordiantes of the local origo (first node
 *       in the list of nodes)
 *  3) Rim-width/thickness in number of nodes.
 *  4) List of all nodes including rim-nodes for this processor.
 */
{
    Nodes<DXQY> newNodes(grid.size());
    newNodes.setup(mfs, rfs);
    return newNodes;
}

#endif // LBNODES_H
