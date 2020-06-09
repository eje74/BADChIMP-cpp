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
    FreeSlipCartesian(const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid) : Boundary<DXQY>(bndNodes, nodes, grid) {}
    void apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const;
    
};

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
  
}
#endif // LBFREESLIPCARTESIAN_H
