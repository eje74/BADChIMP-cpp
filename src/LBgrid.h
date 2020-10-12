#ifndef LBGRID_H
#define LBGRID_H

#include <vector>
#include <numeric>
#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "Input.h"
#include "LBvtk.h"
//#include "Geo.h"

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
};


//template <typename DXQY>
//Grid<DXQY>::Grid(const int nNodes) :nNodes_(nNodes), neigList_(nNodes_ * DXQY::nQ, 0), pos_(nNodes_ * DXQY::nD, -1)
  /* Constructor of a Grid object. Allocates memory
   *  for the neighbor list (neigList_) and the positions (pos_)
   *  pos_ is initated to -1 so that inital check for 'pos inside domain' will return false.
   * Usage:
   *   Grid<D2Q9> grid(number_of_nodes);
   */
//{
//}


template <typename DXQY>
Grid<DXQY>::Grid(LBvtk<DXQY> &vtk) :nNodes_(vtk.endNodeNo()), neigList_(nNodes_ * DXQY::nQ, 0), pos_(nNodes_ * DXQY::nD, -1)
{
    // Set grid's positions
    vtk.toPos();
    for (int n=vtk.beginNodeNo(); n < vtk.endNodeNo(); ++n) {
        // Needed to add template as the compiler may think that < could refere to larger than ...?
        std::vector<int> pos = vtk.template getPos<int>(); 
        addNodePos(pos, n);
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







// NBNBNBNBNB****** OLD CODE BELOW **************

// Her lager vi sikkert en ny
// Hva skal Grid inneholder?
// - Nabonoder i hver gridretning
// - Oversikt over bulknoder
// - Oversikt over 'boundary nodes'

template <typename DXQY>
class GridRegular  // Use ghost nodes. Should be a child a master Grid class in the finished code
{
public:
    // Grid constructor and destructor
    GridRegular<DXQY>(const int nX, const int nY);
    ~GridRegular();

    // Getter for number of elements (including ghost nodes)
    int nElements() const;

    // Returns the position of the element (without ghost nodes)
    void position(int &xNo, int &yNo, const int elementNo) const;

    // Returns the element number from a position (without ghost nodes)
    int element(const int xNo, const int yNo) const; // The user should not need to care about the ghost nodes

    // Returns the element number the neighbor of _elementNo_ in direction _qDirection_
    int neighbor(const int qDirection, const int elementNo) const; // This should be generic to all grids

    // Returns the periodic node of _elementNo_ in direction _qDirection_
    // NB! This function is probably just temporary. Should be fixed in preliminary work.
    int periodicNeighbor(const int qDirection, const int elementNo) const;

private:
    const int nX_;
    const int nY_;
    const int nElements_;
    int neighborStride[DXQY::nQ];
};

// Constructor and destructor
template <typename DXQY>
GridRegular<DXQY>::GridRegular(const int nX, const int nY): nX_(nX + 2), nY_(nY + 2), nElements_(nX_ * nY_)  // Input can also be one or many files ...
{
    //lattice_ = lattice;
    // Setup neighbor list
    for (int q = 0; q < DXQY::nQ; ++q)
        neighborStride[q] = DXQY::c(q, 0) + nX_ * DXQY::c(q, 1);
}

template <typename DXQY>
GridRegular<DXQY>::~GridRegular()
{}


// Getter for number of elements (including ghost nodes)
template <typename DXQY>
inline int GridRegular<DXQY>::nElements() const
{
    return nElements_;
}

// Returns the position of the element (without ghost nodes)
template <typename DXQY>
inline void GridRegular<DXQY>::position(int &xNo, int &yNo, const int elementNo) const // Returns the position of the element
{
    xNo = (elementNo % nX_) - 1;
    yNo = (elementNo / nX_) - 1;
}

// Returns the element number from a position (without ghost nodes)
template <typename DXQY>
inline int GridRegular<DXQY>::element(const int xNo, const int yNo) const // The user should not need to care about the ghost nodes
{
    return (xNo + 1) + (yNo + 1) * nX_;
}

// Returns the element number the neighbor of _elementNo_ in direction _qDirection_
template <typename DXQY>
inline int GridRegular<DXQY>::neighbor(const int qDirection, const int elementNo) const// This should be generic to all grids
{
    return elementNo + neighborStride[qDirection];
}

// Returns the periodic node of _elementNo_ in direction _qDirection_
// NB! This function is probably just temporary. Should be fixed in preliminary work.
template <typename DXQY>
inline int GridRegular<DXQY>::periodicNeighbor(const int qDirection, const int elementNo) const
{
    int xNo, yNo;
    position(xNo, yNo, elementNo);
    int xNeigNo = xNo + DXQY::c(qDirection, 0);
    int yNeigNo = yNo + DXQY::c(qDirection, 1);
    /* elementNo = (xNo + 1) + nX_ * (yNo + 1) */
    xNeigNo %= nX_ - 2;
    xNeigNo = (xNeigNo < 0) ? (xNeigNo + nX_ - 2) : (xNeigNo);
    yNeigNo %= nY_ - 2;
    yNeigNo = (yNeigNo < 0) ? (yNeigNo + nY_ - 2) : (yNeigNo);
    return neighbor(DXQY::reverseDirection(qDirection), element(xNeigNo, yNeigNo));
}

#endif // LBGRID_H
