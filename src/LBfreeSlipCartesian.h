#ifndef LBFREESLIPCARTESIAN_H
#define LBFREESLIPCARTESIAN_H

#include "LBglobal.h"
#include "LBboundary.h"
#include "LBgrid.h"
#include "LBfield.h"


/*********************************************************
 * class FREESLIPCARTESIAN: class that performes the
 * free slip boundary condition for static walls directed 
 * along a Cartesian direction, and has
 * the Boundary class as a parent.
 *
 * Here we have assumed that boundary solids are included
 * in the data structure so that we can extract the
 * f-values that are propagated from the solid into the
 * wall.
 *
 * FreeSlipCartesian.apply(...) needs to be run straight 
 * after propagation.
 *
 *********************************************************/

template <typename DXQY>
class FreeSlipCartesian : public Boundary<DXQY>
{
public:
    FreeSlipCartesian(const std::vector<int> &normVec, const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid);     
    void apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const;
    
private:    
    const std::vector<int> n_vec;  // Normal vector
    int q_wall;
    std::vector<int> beta_ref;  // List of reflected beta values
};


template <typename DXQY>
FreeSlipCartesian<DXQY>::FreeSlipCartesian(const std::vector<int> &normVec, const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
: n_vec(normVec.begin(), normVec.end()),  Boundary<DXQY>(bndNodes, nodes, grid), beta_ref(DXQY::nQ)
{
    // Set the direction of the free slip wall
    q_wall = this->dirRev(DXQY::c2q(n_vec));

    // Setup the reflecte
    for (int q = 0; q < DXQY::nQ; ++q ) {
        std::vector<int> c_ref(DXQY::nQ, 0);
        for (int d = 0; d < DXQY::nD; ++d) {
            c_ref[d] = DXQY::c(q, d) - 2 * DXQY::cDot(q, n_vec) * n_vec[d];
        }
        beta_ref[q] = DXQY::c2q(c_ref);
    }
}


template <typename DXQY>
inline void FreeSlipCartesian<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const
/* apply : performs the free slip condition, at the bondary nodes.
 *
 * fieldNo : the lB-field number
 * f       : the field object
 * grid    : grid object
 *
 * Use 'this->' to access functions and variables in the parent class, Boundary<DXQY>.
 *
 */
{
    for (int n = 0; n < this->nBoundaryNodes_; ++n) {
        int node = this->nodeNo(n);
        int node_wall = grid.neighbor(q_wall, n);
        
        for (auto beta: this->beta(n))
        {
            f(fieldNo, beta, node) = f(fieldNo, beta_ref[beta], node_wall);
        }
    }
}


#endif // LBFREESLIPCARTESIAN_H
