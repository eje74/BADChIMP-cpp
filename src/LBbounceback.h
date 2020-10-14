#ifndef LBBOUNCEBACK_H
#define LBBOUNCEBACK_H

#include "LBboundary.h"

/* Direction classification:
 * 
 *  - beta          : Unknown
 *  - beta_revers   : Known
 *
 *  - gamma         : Known
 *  - gamma_reverse : Known
 *
 *  - delta         : Unknown
 *  - delta_reverse : Unknown
 * 
 */ 

template<typename DXQY>
class SolidBounceBack: public Boundary<DXQY> 
{
public:
    SolidBounceBack(const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid) : Boundary<DXQY>(bndNodes, nodes, grid) {}
    void apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const;
};


template<typename DXQY>
void SolidBounceBack<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const
{
    for (int bndNo = 0; bndNo < this->size(); ++bndNo) {
        int nodeNo = this->nodeNo(bndNo);
        for (auto beta: this->beta(bndNo)) {
            int betaRev = this->dirRev(beta);
            f(fieldNo, beta, grid.neighbor(beta, nodeNo)) = f(fieldNo, betaRev, nodeNo);
        }
        for (auto gamma: this->gamma(bndNo)) {
            int gammaRev = this->dirRev(gamma);
            f(fieldNo, gamma, grid.neighbor(gamma, nodeNo)) = f(fieldNo, gammaRev, nodeNo);
            f(fieldNo, gammaRev, grid.neighbor(gammaRev, nodeNo)) = f(fieldNo, gamma, nodeNo);
        }
    }
}

#endif