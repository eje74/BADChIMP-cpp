#ifndef LBGLOBAL_H
#define LBGLOBAL_H

#include <limits>
#include <vector>
#include <valarray>
#include <algorithm>
#include <numeric>

typedef double lbBase_t;
#define SQRT2 1.4142135623730950488
constexpr lbBase_t lbBaseEps = std::numeric_limits<lbBase_t>::epsilon();

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

#endif // LBGLOBAL_H
