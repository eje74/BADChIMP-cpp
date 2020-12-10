#ifndef LBTURBINEBOUNDARY
#define LBTURBINEBOUNDARY

#include "../src/LBglobal.h"
#include "../src/LBboundary.h"
#include "../src/LBgrid.h"
#include "../src/LBnodes.h"
#include "../src/LBfield.h"


template<typename DXQY>
class InletBoundary : public Boundary<DXQY>
{
public:
    InletBoundary(const std::vector<lbBase_t> &normalVector, const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
        :Boundary<DXQY>(bndNodes, nodes, grid) 
    {


    }
};

#endif
