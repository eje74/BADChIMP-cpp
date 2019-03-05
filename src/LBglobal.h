#ifndef LBGLOBAL_H
#define LBGLOBAL_H

#include <limits>

typedef double lbBase_t;
#define SQRT2 1.4142135623730950488
constexpr lbBase_t lbBaseEps = std::numeric_limits<lbBase_t>::epsilon();

#endif // LBGLOBAL_H
