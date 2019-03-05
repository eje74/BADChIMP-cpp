#ifndef LBGRID_H
#define LBGRID_H

#include <vector>
#include "LBglobal.h"
//#include "LBd2q9.h"
#include "LBlatticetypes.h"
#include "Geo.h"

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
 *
 *********************************************************/
template <typename DXQY>
class Grid
{
public:
    Grid(const int nNodes);  // Constructor
    ~Grid();  // Destructor
    void setup(const Geo &geo);
    int neighbor(const int qNo, const int nodeNo) const;  // See general comment
    int* neighbor(const int nodeNo) const;  // Check if this is in use. Possibly redundant
    int* pos(const int nodeNo) const;  // See general comment
    void addNeigNode(const int qNo, const int nodeNo, const int nodeNeigNo);  // Adds link
    void addNodePos(const int x, const int y, const int nodeNo);  // adds node position in 2d
    void addNodePos(const int x, const int y, const int z, const int nodeNo); // adds node position in 3d
    void addNodePos(const std::vector<int>& ind, const int nodeNo); // adds node position in n-dim

private:
//public:
    int nNodes_;   // Total number of nodes
    int* neigList_;  // List of neighbors [neigNo(dir=0),neigNo(dir=1),neigNo(dir=2)...]
    int* pos_;  // list of cartesian coordinates [x_1,y_1,z_1,x_2,y_2, ...]
    std::vector<int> neigh_list_;
    std::vector<int> xyz_;

public:
    //void print_pos() const {std::cout << "POS: " << xyz_ << std::endl;};
    inline int get_pos(const int i) const {return pos_[i];};
    inline int num_nodes() const {return nNodes_;};
};


template <typename DXQY>
Grid<DXQY>::Grid(const int nNodes) :nNodes_(nNodes), neigh_list_(nNodes_ * DXQY::nQ), xyz_(nNodes_ * DXQY::nD)
  /* Constructor of a Grid object. Allocates memory
   *  for the neighbor list (neigList_) and the positions (pos_)
   * Usage:
   *   Grid<D2Q9> grid(number_of_nodes);
   */
{
    neigList_ = new int [nNodes_ * DXQY::nQ];
    pos_ = new int [nNodes_ * DXQY::nD];
}


template <typename DXQY>
Grid<DXQY>::~Grid()
/* Grid descructor
*/
{
    delete [] neigList_;
    delete [] pos_;
}

template<typename DXQY>
void Grid<DXQY>::setup(const Geo &geo) {
  const std::vector<int>& labels = geo.get_labels();
  const std::vector<int>& n = geo.local_.n_;
  for (size_t pos=0; pos<labels.size(); ++pos) {
    if (labels[pos] > 0) {
      int nodeNo = labels[pos];
      std::vector<int> xyz = geo.get_index(pos);
      //std::cout << xyz << std::endl;
      addNodePos(xyz, nodeNo); //LBgrid
      for (int q = 0; q < DXQY::nQ; ++q) {
        std::vector<int> nn = (xyz + DXQY::c(q) + n) % n;  // periodicity and non-negative nn
        addNeigNode(q, nodeNo, labels[geo.get_pos(nn,n)]); //LBgrid
      }
    }
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


template <typename DXQY>
void Grid<DXQY>::addNodePos(const int x, const int y, const int nodeNo)
/* Adds the Cartesian indices to Grid's position array
 *
 * x      : 1st index
 * y      : 2nd index
 * nodeNo : current node
 */
{
    pos_[nodeNo * DXQY::nD] = x;
    pos_[nodeNo * DXQY::nD + 1] = y;
}

template <typename DXQY>
void Grid<DXQY>::addNodePos(const int x, const int y, const int z, const int nodeNo)
/* Adds the Cartesian indices to Grid's position array
 *
 * x      : 1st index
 * y      : 2nd index
 * z      : 3rd index
 * nodeNo : current node
 */
{
    pos_[nodeNo * DXQY::nD] = x;
    pos_[nodeNo * DXQY::nD + 1] = y;
    pos_[nodeNo * DXQY::nD + 2] = z;
}

template <typename DXQY>
void Grid<DXQY>::addNodePos(const std::vector<int>& ind, const int nodeNo) {
  auto start = xyz_.begin() + nodeNo*DXQY::nD;
  int i = 0;
  for (auto it=start; it!=start+ind.size(); ++it) {
    *it = ind[i++];
  }
  int a = nodeNo * DXQY::nD;
  for (size_t i=a; i<a+ind.size(); ++i)
    pos_[i] = xyz_[i];
//    xyz_[nodeNo * DXQY::nD] = ind[0];
//    xyz_[nodeNo * DXQY::nD + 1] = ind[1];
//    xyz_[nodeNo * DXQY::nD + 2] = ind[2];
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
inline int* Grid<DXQY>::neighbor(const int nodeNo) const
/* Returns a pointer to nodeNo's neighbor list.
 * Example:
 *  int* list = grid.neighbor(current_node_number);
 *  int neighbor_node = list[qNo]; // node number of the current node's neibhor in qNo direction.
 *
 * nodeNo : currnet node number
 *
 * return : Pointer to neigbhor list.
 */
{
    return &neigList_[nodeNo * DXQY::nQ];
}


template <typename DXQY>
inline int* Grid<DXQY>::pos(const int nodeNo) const
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
    return &pos_[DXQY::nD*nodeNo];
}

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
//    Lattice *lattice_;
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
