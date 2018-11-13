#ifndef LBGRID_H
#define LBGRID_H

#include "LBlattice.h"

// Hva skal Grid inneholder?
// - Nabonoder i hver gridretning
// - Oversikt over bulknoder
// - Oversikt over 'boundary nodes'

class GridRegular  // Use ghost nodes. Should be a child a master Grid class in the finished code
{
public:
    // Grid constructor and destructor
    GridRegular(const int nX, const int nY, Lattice *lattice);
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
    Lattice *lattice_;
};

// Getter for number of elements (including ghost nodes)
inline int GridRegular::nElements() const
{
    return nElements_;
}

// Returns the position of the element (without ghost nodes)
inline void GridRegular::position(int &xNo, int &yNo, const int elementNo) const // Returns the position of the element
{
    xNo = (elementNo % nX_) - 1;
    yNo = (elementNo / nX_) - 1;
}

// Returns the element number from a position (without ghost nodes)
inline int GridRegular::element(const int xNo, const int yNo) const // The user should not need to care about the ghost nodes
{
    return (xNo + 1) + (yNo + 1) * nX_;
}

// Returns the element number the neighbor of _elementNo_ in direction _qDirection_
inline int GridRegular::neighbor(const int qDirection, const int elementNo) const// This should be generic to all grids
{
    return elementNo + (*lattice_).c(qDirection, 0) + nX_ * (*lattice_).c(qDirection, 1);
}

// Returns the periodic node of _elementNo_ in direction _qDirection_
// NB! This function is probably just temporary. Should be fixed in preliminary work.
inline int GridRegular::periodicNeighbor(const int qDirection, const int elementNo) const
{
    int xNo, yNo;
    position(xNo, yNo, elementNo);
    int xNeigNo = xNo + (*lattice_).c(qDirection, 0);
    int yNeigNo = yNo + (*lattice_).c(qDirection, 1);
    /* elementNo = (xNo + 1) + nX_ * (yNo + 1) */
    xNeigNo %= nX_ - 2;
    xNeigNo = (xNeigNo < 0) ? (xNeigNo + nX_ - 2) : (xNeigNo);
    yNeigNo %= nY_ - 2;
    yNeigNo = (yNeigNo < 0) ? (yNeigNo + nY_ - 2) : (yNeigNo);
    return neighbor((*lattice_).reverseDirection(qDirection), element(xNeigNo, yNeigNo));
}



#endif // LBGRID_H
