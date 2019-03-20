#ifndef LBGRID_H
#define LBGRID_H

#include <vector>
#include <numeric>
#include "LBglobal.h"
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
    void setup(std::ifstream &ifs, std::vector<int> &dim, std::vector<int> &origo,  int rim_size);
    int neighbor(const int qNo, const int nodeNo) const;  // See general comment
    int* neighbor(const int nodeNo) const;  // Check if this is in use. Possibly redundant
    int* pos(const int nodeNo) const;  // See general comment
    int& pos(const int nodeNo, const int index);
    void addNeigNode(const int qNo, const int nodeNo, const int nodeNeigNo);  // Adds link
    void addNodePos(const int x, const int y, const int nodeNo);  // adds node position in 2d
    void addNodePos(const int x, const int y, const int z, const int nodeNo); // adds node position in 3d
    void addNodePos(const std::vector<int>& ind, const std::vector<int>& origo, const int nodeNo); // adds node position in n-dim

    static Grid<DXQY> makeObject(std::string fileName);

private:
//public:
    int nNodes_;   // Total number of nodes
    int* neigList_;  // List of neighbors [neigNo(dir=0),neigNo(dir=1),neigNo(dir=2)...]
    int* pos_;  // list of cartesian coordinates [x_1,y_1,z_1,x_2,y_2, ...]
    std::vector<int> neigh_list_;
    std::vector<int> xyz_;

public:
    //void print_pos() const {std::cout << "POS: " << xyz_ << std::endl;};
    inline int get_pos(const int i) const {return pos_[i];}
    inline int num_nodes() const {return nNodes_;}
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

template <typename DXQY>
Grid<DXQY> Grid<DXQY>::makeObject(std::string fileName)
/* Makes a grid object using the information in the file created by our
 * python program, for each mpi-processor.
 *
 * We assume that the file contains:
 *  1) Dimesions of the system (including the rim)
 *  2) The global Cartesian coordiantes of the local origo (first node
 *       in the list of nodes)
 *  3) Rim-width/thickness in number of nodes.
 *  4) List of all nodes including rim-nodes for this processor.
 */
{
    std::ifstream ifs;
    ifs.open(fileName);

    // READ preamble
    //  -- dimension
    std::vector<int> dim(DXQY::nD, 0);
    std::string entry_type;
    ifs >> entry_type;
    for (auto &d: dim) ifs >> d;
    std::cout << entry_type << " = " << dim << std::endl;

    //  -- origo
    std::vector<int> origo(DXQY::nD, 0);
    ifs >> entry_type;
    for (auto &d: origo) ifs >> d;
    std::cout << entry_type << " = " << origo << std::endl;

    //  -- rim
    int rim_width = 0;
    ifs >> entry_type >> rim_width;
    std::cout << entry_type << " = " << rim_width << std::endl;

    std::getline(ifs, entry_type);  // Read the remainer of the rim_width line
    std::getline(ifs, entry_type);  // Read the <lable int> line

    // Finds the largest node label, which is equal to the number
    // of nodes excluding the default node.
    int numNodes = 0;
    int nodeNo;
    auto nEntries = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<int>());
    for (int n=0; n < nEntries; ++n) {
        ifs >> nodeNo;
        numNodes = nodeNo > numNodes ? nodeNo : numNodes;
    }
    std::cout << numNodes << std::endl;
    ifs.close();

    // Allocate memory for all nodes in this processor,
    //  include the default node (node label = 0)
    Grid<DXQY> ret(numNodes+1);

    // Use grid's setup to initiate the neighbor lists and
    // node positions
    ifs.open(fileName);
    // Read the 3 lines of preamble
    for (int l = 0; l < 4; ++l)  std::getline(ifs, entry_type);
    ret.setup(ifs, dim, origo, rim_width);
    ifs.close();

    return ret;
}

/***************************************** SETUP GRID
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
    lbFifo(int data_size): data_(static_cast<size_t>(data_size)), size_(data_size)
    {
        zero_pos_ = 0;
        for (std::size_t i = 0; i < data_.size(); ++i)
            data_[i] = 0;
    }

    int & operator [] (int pos) {return data_[static_cast<size_t>((pos + zero_pos_ + size_) % size_)];}
    void push(int value) {data_[static_cast<size_t>(zero_pos_)] = value; zero_pos_ = (1 + zero_pos_) % size_; }
private:
    std::vector<int> data_;
    int size_;
    int zero_pos_;
};


static bool insideDomain(int pos, std::vector<int> dim, int rim_size)
/* insideDomain return false if a a node at pos is part of the rim,
 *  and true if it is not part of the rim.
 *
 * pos : current position
 * dim : Cartesian dimension including the rim
 * rim_size : size of the rim (in number of nodes)
 *
 */
{
    int ni = pos  % dim[0];

    if ( (dim[0] - rim_size) <= ni ) return false;
    if ( rim_size > ni) return false;

    for (size_t d = 0; d < dim.size() - 1; ++d) {
        pos = pos / dim[d];
        ni = pos % dim[d + 1];
        if ( (dim[d+1] - rim_size) <= ni ) return false;
        if ( rim_size > ni) return false;
    }

    return true;
}

static void getPos(int pos, std::vector<int> dim, int rim_size, std::vector<int> &cartesianPos)
/* Find the cartesian indecies to position pos*/
{
    int ni = pos  % dim[0];
    cartesianPos[0] = ni - rim_size;

    for (size_t d = 0; d < dim.size() - 1; ++d) {
        pos = pos / dim[d];
        ni = pos % dim[d + 1];
        cartesianPos[d+1] = ni - rim_size;
    }


}

template<typename DXQY>
void Grid<DXQY>::setup(std::ifstream &ifs, std::vector<int> &dim, std::vector<int> &origo, int rim_size)
/* reads the input file and setup the grid object.
 *
 * ifs : in file stream. Assuming that all premable is read,
 *       so that it is only 'node label map' left to read
 * dim : Cartesian dimension including the rim
 * rim_size : size of the rim (in number of nodes)
 */
{
    int dir_stride[DXQY::nDirPairs_];  // Number of entries between current node and neighbornodes
    int dir_q[DXQY::nDirPairs_];  // q-direction of read neighbornodes
    int dim_stride[DXQY::nD];  // stride in the spatial indices (flattened matrix)

    // Make the dimension stride varaiables
    // [1, nx, nx*ny, ...]
    dim_stride[0] = 1;
    for (size_t d = 0; d < DXQY::nD - 1; ++d)
        dim_stride[d+1] = dim_stride[d]*dim[d];

    // Setup the neighborstrides
    int max_stride = 0;
    int counter = 0;
    for (int q = 0; q < DXQY::nQNonZero_; ++q) {
        int stride = DXQY::cDot(q, dim_stride);
        if (stride < 0) { // Add entry to stride
            dir_stride[counter] = stride;
            dir_q[counter] = q;
            max_stride = (stride < max_stride) ? stride : max_stride;
            counter += 1;
        }
    }

    // ** Read the rank node label file **
    lbFifo node_buffer(-max_stride);
    int nodeNo;
    auto nEntries = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<int>());
    for (int pos=0; pos < nEntries; ++pos)
    {
        ifs >> nodeNo;
        if ( insideDomain(pos, dim, rim_size) ) { // Update the grid object
            if (nodeNo > 0) { // Only do changes if it is a non-default node
                // Add the cartesian position of the node
                std::vector<int> local_index(DXQY::nD, 0);
                getPos(pos, dim, rim_size, local_index);
                addNodePos(local_index, origo, nodeNo);

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
        // pos += 1;
    }
}


/* template<typename DXQY>
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
 */


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
void Grid<DXQY>::addNodePos(const std::vector<int>& ind, const std::vector<int>& origo,  const int nodeNo)
{
    for (std::size_t d = 0; d < ind.size(); ++d)
        pos(nodeNo, d) = ind[d] + origo[0];

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
