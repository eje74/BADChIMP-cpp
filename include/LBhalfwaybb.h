#ifndef LBHALFWAYBB_H
#define LBHALFWAYBB_H

#include "LBglobal.h"
#include "LBboundary.h"
#include "LBgrid.h"
#include "LBfield.h"


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
class HalfWayBounceBack : public Boundary<DXQY>
{
public:
    HalfWayBounceBack(int nBoundaryNodes) : Boundary<DXQY>(nBoundaryNodes) {}
    void apply(const int fieldNo, LbField &f, const Grid<DXQY> &grid) const;
};

template <typename DXQY>
inline void HalfWayBounceBack<DXQY>::apply(const int fieldNo, LbField &f, const Grid<DXQY> &grid) const
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

        // Bounce back for the beta directions (beta unknow)
        for (int q = 0; q < this->nBeta_[n]; ++q) {
            int beta = this->beta(q, n);
            int beta_rev = this->dirRev(beta);
            f(fieldNo, beta, node) = f(fieldNo, beta_rev, grid.neighbor(beta_rev, node));
        }

        // Bounce back for the delta directions (delta and delta.rev unknown)
        for (int q  = 0; q < this->nDelta_[n]; ++q) {
            int delta = this->delta(q, n);
            int delta_rev = this->dirRev(delta);

            f(fieldNo, delta, node) = f(fieldNo, delta_rev, grid.neighbor(delta_rev, node));
            f(fieldNo, delta_rev, node) = f(fieldNo, delta, grid.neighbor(delta, node));
        }
    }
}


#endif // LBHALFWAYBB_H
