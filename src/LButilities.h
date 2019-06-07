#ifndef LBUTILITIES_H
#define LBUTILITIES_H

#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBgrid.h"

template <typename DXQY>
inline std::valarray<lbBase_t> grad(const ScalarField &sField, const int fieldNum, const int nodeNo, const Grid<DXQY> &grid)
{
    std::valarray<lbBase_t> scalarTmp(DXQY::nQ);
    for (int q = 0; q < DXQY::nQ; ++q) {
        int neigNode = grid.neighbor(q, nodeNo);
        scalarTmp[q] = sField(fieldNum, neigNode);
    }

    return DXQY::grad(scalarTmp);
}

template <typename DXQY, typename T>
inline lbBase_t vecNorm(const T &vec)
{
    return sqrt(DXQY::dot(vec, vec));
}

#endif // LBUTILITIES_H
