#ifndef LBGRID_H
#define LBGRID_H

#include <vector>
#include <numeric>
#include <array>
#include <unordered_map>
#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "../io/Input.h"
#include "LBvtk.h"



// Class for handeling nodeNo functionality. That is, obtaining the node number from the position
template <typename DXQY>   // ND: number of dimensions
class NodeNumber
{
public:
    NodeNumber(LBvtk<DXQY> &vtk);
    template <typename T>
    void insertPair(const T &pos, const int nodeNo) {
        posNodeNoPair_.insert( {multi2flat(pos), nodeNo});
    }
    template <typename T>
    int operator[](const T &pos) {return this->posNodeNoPair_[this->multi2flat(pos)];}
private:
    template <typename T>
    inline int multi2flat(const T &pos);    std::array<int, DXQY::nD> posMul_;
    std::unordered_map<int, int> posNodeNoPair_;
};


template <typename DXQY>
NodeNumber<DXQY>::NodeNumber(LBvtk<DXQY> &vtk)
{
    // Setup position multipliers
    posMul_[0] = 1;
    for (int i = 1; i < DXQY::nD; ++i) {
        posMul_[i] = posMul_[i-1]*vtk.getGlobaDimensions(i-1);
    }

    posNodeNoPair_.reserve(vtk.endNodeNo());
}


template<typename DXQY>
template<typename T>
inline int NodeNumber<DXQY>::multi2flat(const T &pos)
{
    int ret = pos[0];
    for (int i = 1; i < DXQY::nD; ++i) {
        ret += pos[i]*posMul_[i];
    }
    return ret;
}



/*********************************************************
 * class GRID:  Contains node indices (i,j,k) and the
 *   node number of neighbors.
 *
 * Functions:
 *  - Grid::neighbor(lattice direction (q), given node number (n))
 *     returns the node number of the neighbor at position
 *     r_n + c_q
 *  - Grid::pos(given node number)
 *     returns the cartesian coordinates of the node as
 *     an integer tuple (that is a constant integer pointer)
 *  - Grid::type(given node number)
 *     Gives the type of node. Used to initilize boundary
 *     objects. Especially
 *
 *********************************************************/
template <typename DXQY>
class Grid
{
public:
//    Grid(int nNodes);  // Constructor
    Grid(LBvtk<DXQY> &vtk);
    inline int neighbor(const int qNo, const int nodeNo) const;  // See general comment
    inline std::vector<int> neighbor(const int nodeNo) const;
    inline const std::vector<int> pos(const int nodeNo) const;  // See general comment
    inline int& pos(const int nodeNo, const int index);
    template <typename T>
    inline int nodeNo(const T &pos) {return nodeNumbers_[pos];}
    inline const int& pos(const int nodeNo, const int index) const;
    void addNeigNode(const int qNo, const int nodeNo, const int nodeNeigNo);  // Adds link
    void addNeighbors(const std::vector<int> &neigNodes, const int nodeNo);
    void addNodePos(const std::vector<int>& ind, const int nodeNo); // adds node position in n-dim
    inline int size() const {return nNodes_;}
    std::vector<std::vector<int>> getNodePos(const std::vector<int> &nodeNoList) const;
    std::vector<std::vector<int>> getNodePos(const int beginNodeNo, const int endNodeNo) const;
     
    
private:
    int nNodes_;   // Total number of nodes
    std::vector<int> neigList_;  // List of neighbors [neigNo(dir=0),neigNo(dir=1),neigNo(dir=2)...]
    std::vector<int> pos_;
    NodeNumber<DXQY> nodeNumbers_;
    // Get 
};



template <typename DXQY>
Grid<DXQY>::Grid(LBvtk<DXQY> &vtk) :nNodes_(vtk.endNodeNo()), neigList_(nNodes_ * DXQY::nQ, 0), pos_(nNodes_ * DXQY::nD, -1), nodeNumbers_(vtk)
{
    // Set grid's positions
    vtk.toPos();
    for (int n=vtk.beginNodeNo(); n < vtk.endNodeNo(); ++n) {
        // Needed to add template as the compiler may think that < could refere to larger than ...?
        std::vector<int> pos = vtk.template getPos<int>(); 
        addNodePos(pos, n);
        nodeNumbers_.insertPair(pos, n);
    }


    // Set node neighborhoods
    vtk.toNeighbors();
    for (int n = vtk.beginNodeNo(); n < vtk.endNodeNo(); ++n) {
        std::vector<int> neigNodes = vtk.template getNeighbors<int>();
        addNeighbors(neigNodes, n);
    }   

    
}


template <typename DXQY>
void Grid<DXQY>::addNeigNode(const int qNo, const int nodeNo, const int nodeNeigNo)
/* Adds a neighbor link to Grid neigList_.
 *
 * qNo        : lattice direction from nodeNo to nodeNeigNo
 * nodeNo     : current node
 * nodeNeigNo : neighbor node in direction c_qNo.
 */
{
    neigList_[nodeNo * DXQY::nQ + qNo] = nodeNeigNo;
}


template<typename DXQY>
void Grid<DXQY>::addNeighbors(const std::vector<int> &neigNodes, const int nodeNo)
/* Add a list of neighbor nodes
 * 
 * neigNodes: list of neighbors. Number of entries need to match DXQY::nQ
 * nodeNo : current Node
 */
{
   for (int q = 0; q < DXQY::nQ; ++q) {
       addNeigNode(q, nodeNo, neigNodes[q]);
   } 
}


template <typename DXQY>
void Grid<DXQY>::addNodePos(const std::vector<int>& ind,  const int nodeNo)
{
    for (unsigned d = 0; d < ind.size(); ++d)
        pos(nodeNo, d) = ind[d];

}


template <typename DXQY>
inline int Grid<DXQY>::neighbor(const int qNo, const int nodeNo) const
/* Returns the node number of nodeNo's neighbor in the qNo direction.
 *
 * qNo    : lattice direction from nodeNo to its neighbor in qNo direction
 * nodeNo : currnet node number
 *
 * return : node number of the neighbor in the qNo direction
 */
{
    return neigList_[nodeNo * DXQY::nQ + qNo];
}


template <typename DXQY>
inline std::vector<int> Grid<DXQY>::neighbor(const int nodeNo) const
/* Reurns iterator to the neighborhood nodeNo
 *
 * nodeNo : node number
 *
 * return : array of neighbor node numbers
 */
{
    auto pos_begin = neigList_.data() + nodeNo * DXQY::nQ;
    return std::vector<int>(pos_begin, pos_begin + DXQY::nQ);
}


template <typename DXQY>
inline const std::vector<int> Grid<DXQY>::pos(const int nodeNo) const
/* Returns a pointer the the Cartesian positon array of nodeNo.
 * Example:
 * int* indices = grid.pos(current_node_number):
 * int x = indices[0]; // first Cartesian coordinate
 * int y = indices[1]; // second Cartesian coordinate
 * int z = indices[2]; // if 3d then third Cartesian coordinate
 *
 * nodeNo : currnet node number
 *
 * return : Pointer to array of Catesian positions
*/
{
    auto pos_begin = pos_.data() + DXQY::nD*nodeNo;
    return std::vector<int>(pos_begin, pos_begin + DXQY::nD);
}


template <typename DXQY>
inline int& Grid<DXQY>::pos(const int nodeNo, const int index)
/* Returns a reference to nodeNo's cartesian coordinate given by index.

  Example: y-coordinate of node 25 is given by
                Grid.pos(25, 1)

 * int x = indices[0]; // first Cartesian coordinate
 * int y = indices[1]; // second Cartesian coordinate
 * int z = indices[2]; // if 3d then third Cartesian coordinate
 *
 * nodeNo : currnet node number
 * index : Cartesian dimension (eq 0, 1, 2)
 * return : a coordiante refernce
*/
{
    return pos_[DXQY::nD*nodeNo + index];
}


template <typename DXQY>
inline const int& Grid<DXQY>::pos(const int nodeNo, const int index) const
{
    return pos_[DXQY::nD*nodeNo + index];
}


template <typename DXQY>
std::vector<std::vector<int>> Grid<DXQY>::getNodePos(const std::vector<int> &nodeNoList) const
/* Returns a vector of positions of the nodes given in nodeNoList. Each position has  
 * tre elements where the last elemen is one if the system only has two sptatial dimensions.
 * This function is used to supply the node_pos input used by the 
 * Output objects.
 */ 
{
    std::vector<std::vector<int>> nodePos;
    nodePos.reserve(nodeNoList.size());
    for (const auto &nodeNo: nodeNoList) {
        std::vector<int> tmp(3, 0);
        for (int i=0; i < DXQY::nD; ++i)
            tmp[i] = pos(nodeNo, i);
        nodePos.push_back(tmp);
    }
    return nodePos;
}

template <typename DXQY>
std::vector<std::vector<int>> Grid<DXQY>::getNodePos(const int beginNodeNo, const int endNodeNo) const
{
    std::vector<std::vector<int>> nodePos;
    
    // Check of node numer limits
    if ( beginNodeNo <= 0 ) {
        std::cout << "Error in Grid.getNodePos: Smallest node value must be larger than 0. Now it is " << beginNodeNo << "." << std::endl;
        exit(1);
    }
    if ( endNodeNo > size()) {
        std::cout << "Error in Grid.getNodePos: Largest node value must be less than or equal Grid's size." << std::endl;
        exit(1);
    }
    if ( endNodeNo < beginNodeNo) {
        std::cout << "Error in Grid.getNodePos: End node number is smaller than begin node number." << std::endl;
        exit(1);
    }
    
    nodePos.reserve(endNodeNo-beginNodeNo);
    for (int nodeNo = beginNodeNo; nodeNo < endNodeNo; ++nodeNo)
    {
        std::vector<int> tmp(3, 0);
        for (int i=0; i < DXQY::nD; ++i)
            tmp[i] = pos(nodeNo, i);
        nodePos.push_back(tmp);
    }
    return nodePos;
}



#endif // LBGRID_H
