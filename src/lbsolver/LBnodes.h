#ifndef LBNODES_H
#define LBNODES_H

#include <vector>
#include <numeric>
#include <unordered_map>   // Get nodeNo from pos
#include "LBglobal.h"
#include "LBlatticetypes.h"
//#include "../io/Input.h"
#include "LBvtk.h"
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
 * Nodes: -contains additional informatino about a node, eg. rank, node-type
 *        -class Grid will still contain the cartesion indices of the nodes and
 *         the list of node-numbers of the nodes in a node's neighborhood.
 *
 * Work in progress: (FINISHED) moving node-information from class Grid to class Nodes
*/


/*******************************************************
 * NODES type codes
 *  SOLIDS             | FLUIDS
 *  -1: default ghost  |  2: boundary fluid
 *   0: bulk solid     |  3: bulk fluid
 *   1: boundary solid |  [4: MPI boundary]. No longer in use. Use isMpiBoundary()  
 * 
 *******************************************************/
template<typename DXQY>
class Nodes
{
public:
    Nodes(const int nNodes, const int myRank):
        nNodes_(nNodes),
        myRank_(myRank),
        nodeRank_(static_cast<std::size_t>(nNodes), myRank_),
        nodeType_(static_cast<std::size_t>(nNodes), -1)
    {  
    }
    
    Nodes(LBvtk<DXQY> &vtk, const Grid<DXQY> &grid);
    
    inline int size() const {return nNodes_;}

    void setupNodeType(const Grid<DXQY> &grid);
    inline int getType(const int nodeNo) const {return nodeType_[nodeNo];}
    inline int getRank(const int nodeNo) const {return nodeRank_[nodeNo];}
    inline int getTag(const int nodeNo) const {return nodeTag_[nodeNo];}
    inline bool isDefault(const int nodeNo) const {return nodeType_[nodeNo] == -1;}
    inline bool isMyRank(const int nodeNo) const {return nodeRank_[nodeNo] == myRank_;}
    inline bool isSolid(const int nodeNo) const {return (nodeType_[nodeNo] < 2);} /* Solid: -1, 0 or 1*/
    inline bool isBulkSolid(const int nodeNo) const {return (nodeType_[nodeNo] == 0);}
    inline bool isSolidBoundary(const int nodeNo) const {return (nodeType_[nodeNo] == 1);}
    inline bool isFluid(const int nodeNo) const {return (nodeType_[nodeNo] > 1);} /* Fluid: 2, 3 or 4 */
    inline bool isBulkFluid(const int nodeNo) const {return (nodeType_[nodeNo] == 3);}
    inline bool isFluidBoundary(const int nodeNo) const {return (nodeType_[nodeNo] == 2);}
    inline bool isMpiBoundary(const int nodeNo) const {return (nodeRank_[nodeNo] != myRank_) && (!isDefault(nodeNo));}
    inline bool isPressureBoundary(const int nodeNo) const {return (nodeType_[nodeNo] == 5);}

    void addNodeRank(const int nodeRank, const int nodeNo) {nodeRank_[nodeNo] = nodeRank;}
    void addNodeType(const int nodeType, const int nodeNo) {nodeType_[nodeNo] = nodeType;}
    void addNodeTag(const int nodeTag, const int nodeNo) {nodeTag_[nodeNo] = nodeTag;}
    inline void setTag(const int nodeTag, const int nodeNo) {nodeTag_[nodeNo] = nodeTag;}

    // JLV
    // Return a vector where 1 is solid nodes and 0 is not solid nodes
    std::vector<int> geo(const Grid<DXQY>& grid, const LBvtk<DXQY>& vtklb) const {
        std::vector<int> vec(grid.size(), -1);
        for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
            vec[nodeNo] = isSolid(nodeNo) ? 1 : 0;        
        }
        return vec;
    }
    // JLV


private:
    int nNodes_;
    int myRank_;
    std::vector<int> nodeRank_; // Node rank
    std::vector<short int> nodeType_; // The actual node type
    std::vector<short int> nodeTag_; // Tag ie. a pressure boundary
};


template<typename DXQY>
Nodes<DXQY>::Nodes(LBvtk<DXQY> &vtk,const Grid<DXQY> &grid):
        nNodes_(grid.size()),
        myRank_(vtk.getRank()),
        nodeRank_(static_cast<std::size_t>(nNodes_), myRank_),
        nodeType_(static_cast<std::size_t>(nNodes_), -1),
        nodeTag_(static_cast<std::size_t>(nNodes_), -1)
{
    // Set rank for processors shared nodes
    addNodeRank(-1, 0); // Set default to -1
    for (int n = 0; n < vtk.getNumNeigProc(); ++n) {
        int nodeRank = vtk.getNeigRank(n);
        std::vector<int> neigNodes = vtk.getNeigNodesNo(n);
        for (auto nodeNo: neigNodes) {
            addNodeRank(nodeRank, nodeNo);
        }
    }
 
    // Setup node types
    vtk.toAttribute("nodetype");
    for (int n = vtk.beginNodeNo(); n < vtk.endNodeNo(); ++n) {
        int nodeType = vtk.template getScalarAttribute<int>();
        addNodeType(nodeType, n);
    }
    setupNodeType(grid);
    
}


template <typename DXQY>
void Nodes<DXQY>::setupNodeType(const Grid<DXQY> &grid)
/*
 * Setsup the node types based on grid and nodeRank, so these has
 *  to be up to date.
 * 
 * Node type numerical code
 * -1 : default
 *  0 : solid (bulk solid)
 *  1 : solid boundary
 * 
 *  2 : fluid boundary
 *  3 : fluid (bulk fluid)
 *  4 : fluid on another process (This is just a temporary label. Will change to either 2 or 3 during mpi setup)
 */
{
    // Set the default
    //addNodeType(-1, 0);
    nodeType_[0] = -1;
    // Test for rank of default node
    if (getRank(0) != -1) {
        std::cout << "WARNING: Node number 0 has rank " << getRank(0) << ". Should this be -1?" << std::endl;
    }
    // INITIATE TYPES
    // -- Set all solid nodes to 0 and all fluid to 3
    for (int nodeNo = 1; nodeNo < size(); ++nodeNo)
    {
        if (nodeType_[nodeNo] == 0) nodeType_[nodeNo] = 0;
        else nodeType_[nodeNo] = 3;    
    }
    
    for (int nodeNo = 1; nodeNo < size(); ++nodeNo)
    {
        if ( isSolid(nodeNo) )
        {
            // Is the solid node a boundary (needs fluid neighbors)
            // or bulk (only solid neighbors) ?
            bool hasFluidNeig = false;
            for (auto neigNode: grid.neighbor(nodeNo))
            {
                if (isFluid(neigNode))  hasFluidNeig = true;
            }
            nodeType_[nodeNo] =  hasFluidNeig ? 1 : 0;
        }
        else if ( isFluid(nodeNo) )
        {
            // Is the fluid a node on a neighboring process?
            if (myRank_ != getRank(nodeNo))
            {
                nodeType_[nodeNo] = 4;
            } 
            else
            {
                // Is it a fluid boundary node (needs solid neighbors)
                // or is it a bulk fluid (only fluid neighbors)
                bool hasSolidNeig = false;
                for (auto neigNode: grid.neighbor(nodeNo)) 
                {
                    if (isSolid(neigNode))  hasSolidNeig = true;
                }
                nodeType_[nodeNo]  = hasSolidNeig ? 2 : 3;
            }
        }
        else
        {
            std::cout << "WARNING: NodeNo = " << nodeNo << " is neither fluid or solid" << std::endl;
        }
    }
}


#endif // LBNODES_H
