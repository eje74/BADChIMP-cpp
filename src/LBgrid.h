#ifndef LBGRID_H
#define LBGRID_H

#include <vector>
#include <numeric>
#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "Input.h"
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
 *  - Grid::type(given node number)
 *     Gives the type of node. Used to initilize boundary
 *     objects. Especially
 *
 *********************************************************/
template <typename DXQY>
class Grid
{
public:
    Grid(const int nNodes);  // Constructor

    void setup(MpiFile<DXQY> &mfs);

    inline int neighbor(const int qNo, const int nodeNo) const;  // See general comment
    inline std::vector<int> neighbor(const int nodeNo) const;
    inline const std::vector<int> pos(const int nodeNo) const;  // See general comment
    inline int& pos(const int nodeNo, const int index);
    inline const int& pos(const int nodeNo, const int index) const;
    void addNeigNode(const int qNo, const int nodeNo, const int nodeNeigNo);  // Adds link
    void addNodePos(const std::vector<int>& ind, const int nodeNo); // adds node position in n-dim
//    void addNodeType(const int type, const int nodeNo);
 //   inline int getType(const int nodeNo) const {return nodeType_[nodeNo];}
 //   inline int getRank(const int nodeNo) const {return nodeType_[nodeNo]-1;}
    inline int size() const {return nNodes_;}

    static Grid<DXQY> makeObject(MpiFile<DXQY> &mfs);

private:
    int nNodes_;   // Total number of nodes
    std::vector<int> neigList_;  // List of neighbors [neigNo(dir=0),neigNo(dir=1),neigNo(dir=2)...]
    std::vector<int> pos_;
};


template <typename DXQY>
Grid<DXQY>::Grid(const int nNodes) :nNodes_(nNodes), neigList_(nNodes_ * DXQY::nQ), pos_(nNodes_ * DXQY::nD)
  /* Constructor of a Grid object. Allocates memory
   *  for the neighbor list (neigList_) and the positions (pos_)
   * Usage:
   *   Grid<D2Q9> grid(number_of_nodes);
   */
{}


template <typename DXQY>
Grid<DXQY> Grid<DXQY>::makeObject(MpiFile<DXQY> &mfs)
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

    // Finds the largest node label, which is equal to the number
    // of nodes excluding the default node.
    int numNodes = 0;
    for (std::size_t n=0; n < mfs.size(); ++n) {
        auto nodeNo = mfs.template getVal<int>();
        numNodes = nodeNo > numNodes ? nodeNo : numNodes;
    }

    // Allocate memory for all nodes in this processor,
    //  include the default node (node label = 0)
    Grid<DXQY> ret(numNodes+1);

    // Use grid's setup to initiate the neighbor lists and
    // node positions
    mfs.reset();
    ret.setup(mfs);

    // Must reset files so that they can be read by the next function
    mfs.reset();
    return ret;
}


/*****************************************  GRID
 * Setup grid using the lokal processor-file
 *
 * We assume that the file is ordered as a
 * structured Cartesian grid . We assume that
 * the size of the file is given as
 * nX, nY, ... That is, the size of the
 * axis in the Cartesian grid
 * (including the ghost node rim)
 * We assume that the size of the rim
 * is given. If it is not given we assume
 * that it is 1.
 *
 * We assume that the position, 'pos' of the node
 * with indecies nx, ny, nz,... is given
 * by :
 *     'pos' = nx + nX*[ny + nY*[nz + ...]]
 *
 * We can keep a list of only the read values that we need
 * to create the grid neighborhood.
 *
 * Algorithm:
 * -------------
 * Calculate the strides for each directions:
 * c_x + nX*[c_y + nY*[c_z]].
 * We will only have access to the ones that
 * has a negative stride. We assume that we
 * do not need to store the current position.
 * Hence we can take the size of the vector
 * holding the read node labels to be:
 * -min(strides) - 1.
 *
 * When position is read we can then  find
 * all its neighbors with negative strides,
 * and assuming that all lattice vectors has
 * a bounce back direction we can the also
 * update that direction for the read node label.
 *
 * NB! We will not update any values to the zero label
 * as this can delete grid information.
 *
 * When a node is read and all information is
 * updated we will add this node the the first
 * of the list and reomve the first read (if
 * the list if full). I could not find
 * a std container that did this so, we
 * need to make one.
******************************************/

class lbFifo {
    /* lbFifo (first in, first out) is a type of queue container used when
     * reading grid information from a file.
     * So that we can access prevously read values.
     *
     * It works in the following way:
     * lbFifo[-1] is the previously read value
     * lbFifo[-data_size] is the oldest value that is remember
     *
     */
public:
    lbFifo(int data_size): data_(static_cast<size_t>(data_size), 0), size_(data_size)
    {
        zero_pos_ = 0;
    }

    int & operator [] (int pos) {return data_[static_cast<size_t>((pos + zero_pos_ + size_) % size_)];}
    void push(int value) {data_[static_cast<size_t>(zero_pos_)] = value; zero_pos_ = (1 + zero_pos_) % size_; }
private:
    std::vector<int> data_;
    int size_;
    int zero_pos_;
};


template<typename DXQY>
void Grid<DXQY>::setup(MpiFile<DXQY> &mfs)
/* reads the input file and setup the grid object.
 *
 * mfs : local node number file
 * rfs : rank number file
 */
{
    int dir_stride[DXQY::nDirPairs_];  // Number of entries between current node and neighbornodes
    int dir_q[DXQY::nDirPairs_];  // q-direction of read neighbornodes
    int dim_stride[DXQY::nD];  // stride in the spatial indices (flattened matrix)

    // Make the dimension stride varaiables
    // [1, nx, nx*ny, ...]
    dim_stride[0] = 1;
    for (size_t d = 0; d < DXQY::nD - 1; ++d)
        dim_stride[d+1] = dim_stride[d]*mfs.dim(d);

    // Setup the neighborstrides
    int max_stride = 0;
    int counter = 0;
    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
        int stride = DXQY::cDot(q, dim_stride);
        if (stride < 0) { // Add entry to stride. Only negative strides so that 1) its read into the buffer
                          //                                                    2) we do not double count a link.
            dir_stride[counter] = stride;
            dir_q[counter] = q;
            max_stride = (stride < max_stride) ? stride : max_stride;
            counter += 1;
        }
    }

    // ** Read the rank node label file **
    lbFifo node_buffer(-max_stride);

    for (int pos=0; pos < static_cast<int>(mfs.size()); ++pos)
    {
        auto nodeNo = mfs.template getVal<int>();

        if ( mfs.insideDomain(pos) ) { // Update the grid object
            if (nodeNo > 0) { // Only do changes if it is a non-default node
                // Add the cartesian position of the node
                std::vector<int> localCartInd(DXQY::nD, 0);
                mfs.getPos(pos, localCartInd);
                addNodePos(localCartInd, nodeNo);

                // Update a link in the nodes neighborhood
                addNeigNode(DXQY::nQNonZero_, nodeNo, nodeNo); // Add the rest particle
                for (int q = 0; q < DXQY::nDirPairs_; ++q) { // Lookup all previously read neighbors
                    int neigNodeNo = node_buffer[dir_stride[q]];
                    if (neigNodeNo > 0) {
                        addNeigNode(dir_q[q], nodeNo, neigNodeNo); // Add the link information
                        addNeigNode(DXQY::reverseDirection( dir_q[q] ), neigNodeNo, nodeNo); // Add the reverse link information
                    }
                }
            }
        }
        node_buffer.push(nodeNo);
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
