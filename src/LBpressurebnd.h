#ifndef LBPRSSUREBND_H
#define LBPRSSUREBND_H


#include "LBglobal.h"
#include "LBboundary.h"
#include "LBgrid.h"
#include "LBfield.h"

template <typename DXQY>
class PressureBnd : public Boundary<DXQY>
{
public:
    PressureBnd(const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid) : Boundary<DXQY>(bndNodes, nodes, grid) {}
    void apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid, const ScalarField &rho) const;
};

template<typename DXQY>
void PressureBnd<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid, const ScalarField &rho) const
/* apply : performs the pressure boundary conditions.
 *
 * fieldNo : the lB-field number
 * f       : the field object
 * grid    : grid object
 *
 * Use 'this->' to access functions and variables in the parent class, Boundary<DXQY>.
 *
 */
{
    for (int n=0; n < this->size(); ++n) {
        int nodeNo = this->nodeNo(n);
        for (auto beta: this->beta(n)) {
            f(fieldNo, beta, grid.neighbor(beta, nodeNo)) = DXQY::w[beta]*rho(fieldNo, nodeNo);
        }
        for (auto delta: this->delta(n)) {
            auto delta_rev = this->dirRev(delta);
            f(fieldNo, delta, grid.neighbor(delta, nodeNo)) = DXQY::w[delta]*rho(fieldNo, nodeNo);
            f(fieldNo, delta_rev, grid.neighbor(delta_rev, nodeNo)) = DXQY::w[delta_rev]*rho(fieldNo, nodeNo);
        }
    }
}

template<typename DXQY>
class InletOutlet : public Boundary<DXQY>
{
public:
    InletOutlet(const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid) : Boundary<DXQY>(bndNodes, nodes, grid) {}
    void apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid, const lbBase_t &rho, const std::vector<lbBase_t> vel) const;
};

template<typename DXQY>
void InletOutlet<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid, const lbBase_t &rho, const std::vector<lbBase_t> vel) const
{
    lbBase_t u_sq = DXQY::dot(vel, vel);
    std::valarray<lbBase_t> cu = DXQY::cDotAll(vel);
    for (int n=0; n < this->size(); ++n) {
        int nodeNo = this->nodeNo(n);
        for (auto beta: this->beta(n)) {
            f(fieldNo, beta, grid.neighbor(beta, nodeNo)) = rho * DXQY::w[beta]*( 1.0 + DXQY::c2Inv*cu[beta] + DXQY::c4Inv0_5*(cu[beta]*cu[beta] - DXQY::c2*u_sq) );
        }
        for (auto delta: this->delta(n)) {
            auto delta_rev = this->dirRev(delta);
            f(fieldNo, delta, grid.neighbor(delta, nodeNo)) = rho * DXQY::w[delta]*( 1.0 + DXQY::c2Inv*cu[delta] + DXQY::c4Inv0_5*(cu[delta]*cu[delta] - DXQY::c2*u_sq) );
            f(fieldNo, delta_rev, grid.neighbor(delta_rev, nodeNo)) = rho * DXQY::w[delta_rev]*( 1.0 + DXQY::c2Inv*cu[delta_rev] + DXQY::c4Inv0_5*(cu[delta_rev]*cu[delta_rev] - DXQY::c2*u_sq) );
        }
    }
}

#endif
