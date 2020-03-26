#ifndef LBUTILITIES_H
#define LBUTILITIES_H

#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBgrid.h"
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort

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

template <typename DXQY>
inline lbBase_t divGrad(const ScalarField &sField, const int fieldNum, const int nodeNo, const Grid<DXQY> &grid)
{
    std::valarray<lbBase_t> scalarTmp(DXQY::nQ);
    for (int q = 0; q < DXQY::nQ; ++q) {
        int neigNode = grid.neighbor(q, nodeNo);
        scalarTmp[q] = sField(fieldNum, neigNode);
    }

    return DXQY::divGrad(scalarTmp);
}

template <typename DXQY, typename T>
inline lbBase_t vecNorm(const T &vec)
{
    return sqrt(DXQY::dot(vec, vec));
}

template <typename T>
inline std::valarray<T> inputAsValarray(Block& input)
{    std::vector<T> tmpVec = input;
     return std::valarray<lbBase_t>(tmpVec.data(), tmpVec.size());
}

template <typename DXQY>
void setFieldFromFile(MpiFile<DXQY> &valueFile, MpiFile<DXQY> &localFile)
{

}

template <typename T>
std::vector<std::size_t> sort_indexes(const std::vector<T> &v)
/* sorting the indexes so that
 *
 *  for (auto i: idx)
 *     std::cout << v[i] << std::endl
 *
 * prints the sorted list.
 * Taken form stackowerflow: https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
 *
 *
 */
{

  // initialize original index locations
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

#endif // LBUTILITIES_H
