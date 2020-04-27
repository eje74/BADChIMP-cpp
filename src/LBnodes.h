#ifndef LBNODES_H
#define LBNODES_H

#include <vector>
#include <numeric>
#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "Input.h"
#include "LBgrid.h"

/*******************************************************
 * Node LABEL functions. Each node in the computational
 * domains should be given its unique label-tag (int number)
 *
 *  Example:
 *  geo-matrix:     label-matrix:
 *    4 4 4            0  0  0
 *    3 3 3           10 11 12
 *    1 1 1            1  2  3
 *    0 0 0            4  5  6
 *    1 1 1            7  8  9
 *    3 3 3           13 14 15
 *
 *  Here we have made the choices that
 *  - 0 is a default bulk solid tag. We will allocate
 *    memory for it, but we assume that it contains
 *    only dummy values
 *  - We also want the bulk node values to be numbered
 *    consequently from 1 (number one). In this example
 *    both fluid bulk and fluid boundary are treated as
 *    bulk nodes in the LB iterations for macroscopic
 *    varaible evaluation, collision and propagation.
 *
 *******************************************************/

/*
 * Nodes: -contains additional informatino about a node, eg. rank, node type etc.
 *        -class Grid will still contain the cartesion indices of the nodes and
 *         the list of node-numbers of the nodes in a node's neighborhood.
 *
 * Work in progress: (FINISHED) moving node-information from class Grid to class Nodes
*/

class newNodes;
template<typename DXQY>
class Nodes
{
public:
    Nodes(const int nNodes, const int myRank):
        nNodes_(nNodes),
        myRank_(myRank),
        nodeRank_(static_cast<std::size_t>(nNodes), -1),
        nodeType_(static_cast<std::size_t>(nNodes), -1)
    {}
    inline int size() const {return nNodes_;}
    void setup(MpiFile<DXQY> &mfs, MpiFile<DXQY> &rfs, const Grid<DXQY> &grid);
    inline int getType(const int nodeNo) const {return nodeType_[nodeNo];}
    inline int getRank(const int nodeNo) const {return nodeRank_[nodeNo];}
    inline bool isMyRank(const int nodeNo) const {return nodeRank_[nodeNo] == myRank_;}
    inline bool isSolid(const int nodeNo) const {return (nodeType_[nodeNo] < 2);} /* Solid: -1, 0 or 1*/
    inline bool isBulkSolid(const int nodeNo) const {return (nodeType_[nodeNo] == 0);}
    inline bool isSolidBoundary(const int nodeNo) const {return (nodeType_[nodeNo] == 1);}
    inline bool isFluid(const int nodeNo) const {return (nodeType_[nodeNo] > 1);} /* Fluid: 2, 3 or 4 */
    inline bool isBulkFluid(const int nodeNo) const {return (nodeType_[nodeNo] == 3);}
    inline bool isFluidBoundary(const int nodeNo) const {return (nodeType_[nodeNo] == 2);}
    inline bool isMpiBoundary(const int nodeNo) const {return (nodeType_[nodeNo] == 4);}

    static Nodes<DXQY> makeObject(MpiFile<DXQY> &mfs, MpiFile<DXQY> &rfs, const int myRank, const Grid<DXQY> &grid);

private:
    int nNodes_;
    int myRank_;
    std::vector<int> nodeRank_; // The actual node type
    std::vector<short int> nodeType_; // The actual node type

};


template<typename DXQY>
void Nodes<DXQY>::setup(MpiFile<DXQY> &mfs, MpiFile<DXQY> &rfs, const Grid<DXQY> &grid)
/*
 * mfs : local node number file
 * rfs : rank number file
 *
 * Node type
 * -1 : default
 *  0 : solid (bulk solid)
 *  1 : solid boundary
 *  2 : fluid boundary
 *  3 : fluid (bulk fluid)
 *  4 : fluid on another process
 */
{
    mfs.reset();
    rfs.reset();


    // Assign each node as either fluid (3) or solid (0)
    // Set the rank of each node
    for (int pos=0; pos < static_cast<int>(mfs.size()); ++pos)
    {
        auto nodeNo = mfs.template getVal<int>(); // Gets node number. The template name is needed
        auto nodeType =  rfs.template getVal<int>(); // get the node type

        // If rank is already set to myRank do not change it
        if (nodeNo > 0) {
            // Set the rank of the node (-1 if no rank is set, eg. for the default node)
            if (nodeRank_[nodeNo] != myRank_)
                nodeRank_[nodeNo] = nodeType - 1;
            // Set node to solid (0) or fluid (3)
            if (nodeType == 0) nodeType_[nodeNo] = 0;
            else nodeType_[nodeNo] = 3;
        }
    }

    mfs.reset();
    rfs.reset();

    // Refine the nodeType, depending on neighborhood nodes.
    for (int nodeNo = 1; nodeNo < grid.size(); ++nodeNo) {
        if (isSolid(nodeNo)) {
            bool hasFluidNeig = false;
            for (auto neigNode: grid.neighbor(nodeNo)) {
                if (isFluid(neigNode))  hasFluidNeig = true;
            }
            if (hasFluidNeig) nodeType_[nodeNo] = 1;
            else nodeType_[nodeNo] = 0;
        }
        else if (isFluid(nodeNo)) {
            if (getRank(nodeNo) != myRank_) {
                nodeType_[nodeNo] = 4;
            } else {
                bool hasSolidNeig = false;
                for (auto neigNode: grid.neighbor(nodeNo)) {
                    if (isSolid(neigNode))  hasSolidNeig = true;
                }
                if (hasSolidNeig)  nodeType_[nodeNo] = 2;
                else nodeType_[nodeNo] = 3;
            }
        } else {
            std::cout << "WARNING: NodeNo = " << nodeNo << " is neither fluid or solid" << std::endl;
        }
    }
}

template<typename DXQY>
Nodes<DXQY> Nodes<DXQY>::makeObject(MpiFile<DXQY> &mfs, MpiFile<DXQY> &rfs, const int myRank, const Grid<DXQY> &grid)
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
    Nodes<DXQY> newNodes(grid.size(), myRank);

    mfs.reset();
    rfs.reset();

    newNodes.setup(mfs, rfs, grid);

    mfs.reset();
    rfs.reset();

    return newNodes;
}

#endif // LBNODES_H
