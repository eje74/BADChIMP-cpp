#ifndef LBBOUNCEBACKBASIC_H
#define LBBOUNCEBACKBASIC_H

#include <vector>
#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBgrid.h"
#include "LBboundarybasic.h"


template <typename DXQY>
class BounceBackBasic
{
public:
    BounceBackBasic(const std::vector<int> & boundaryNodes,
                    const Nodes<DXQY>      & nodes, 
                    const Grid<DXQY>       & grid);
                    
    void apply(const int          fieldNo, 
               LbField<DXQY>    & f, 
               const Grid<DXQY> & grid) const;
                    
private:
    BoundaryBasic<DXQY> bnd_;
};



template <typename DXQY>
BounceBackBasic<DXQY>::BounceBackBasic(const std::vector<int> & boundaryNodes,
                                       const Nodes<DXQY>      & nodes, 
                                       const Grid<DXQY>       & grid
                                       ): 
                                       bnd_(boundaryNodes, nodes, grid)
{
}



template <typename DXQY>
void BounceBackBasic<DXQY>::apply(const int          fieldNo, 
                                  LbField<DXQY>    & f, 
                                  const Grid<DXQY> & grid
                                  ) const
{
    for (int bndNo = 0; bndNo < bnd_.size(); ++bndNo)
    {
        // Bondary node number
        auto nodeNo = bnd_.nodeNo(bndNo);

        for (auto q: bnd_.unknownDirs(bndNo))
        {
            auto qRev = bnd_.dirRev(q);
            auto wallNodeNo = grid.neighbor(qRev, nodeNo);
            f(fieldNo, q, nodeNo) = f(fieldNo, qRev, wallNodeNo);
        }
    }
}



#endif