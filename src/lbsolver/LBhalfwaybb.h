#ifndef LBHALFWAYBB_H
#define LBHALFWAYBB_H

#include "LBglobal.h"
//#include "LBboundary.h"
#include "LBgrid.h"
#include "LBfield.h"
#include "LBhalfwayhelperclass.h"

/*********************************************************
 * class HALFWAYBOUNCEBACK: class that performes the
 *  half way bounce back for static walls, and has
 *  the Boundary class as a parent.
 *
 * Here we have assumed that boundary solids are included
 * in the data structure so that we can extract the
 * f-values that are propagated from the solid into the
 * wall.
 *
 * HalfWayBounceBack.apply(...) needs to be run straight
 * after propagation.
 *
 *********************************************************/
template <typename DXQY>
class HalfWayBounceBack : public BoundaryHalwWayHelper<DXQY>
{
public:
    HalfWayBounceBack(const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid) : BoundaryHalwWayHelper<DXQY>(bndNodes, nodes, grid) {}
//    HalfWayBounceBack(Boundary<DXQY> base) : Boundary<DXQY>(base.size()) {}
    void apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const;
    void apply(LbField<DXQY> &f, const Grid<DXQY> &grid) const;
    
};



template <typename DXQY>
inline void HalfWayBounceBack<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const
/* apply : performs the half way bounce back, the bondary nodes.
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
	
        for (auto beta: this->beta(n)) {
            int beta_rev = this->dirRev(beta);
            f(fieldNo, beta, node) = f(fieldNo, beta_rev, grid.neighbor(beta_rev, node));
        }

        for (auto delta: this->delta(n)) {
            int delta_rev = this->dirRev(delta);
            f(fieldNo, delta, node) = f(fieldNo, delta_rev, grid.neighbor(delta_rev, node));
            f(fieldNo, delta_rev, node) = f(fieldNo, delta, grid.neighbor(delta, node));
        }
    }
}

template <typename DXQY>
inline void HalfWayBounceBack<DXQY>::apply(LbField<DXQY> &f, const Grid<DXQY> &grid) const
{
    for (int n=0; n < f.num_fields(); ++n) {
        apply(n, f, grid);
    }
}

#endif // LBHALFWAYBB_H
