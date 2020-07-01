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

template <typename DXQY, typename T>
inline std::valarray<lbBase_t> calcfeq(const lbBase_t& rho, const lbBase_t& u_sq, const T &cu)
/* calcfeq : calculates the equilibrium distribution
 *
 * rho      : density
 * u_sq     : square of the velocity
 * cu       : array of scalar product of all lattice vectors and the velocity.
 * omegaBGK : array of the BGK-collision term in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = rho * DXQY::w[q]*(1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq) );
    }
    return ret;
}

template <typename DXQY, typename T1, typename T2>
  inline std::valarray<lbBase_t> calcRegDist(const T1 &feq, const T2 &PiNeqLowTri)
/* calcOmegaBGK : sets the BGK-collision term in the lattice boltzmann equation
 *
 * feq          : precalculated array of local equilibriums distributions.
 * PiNeqLowTri    : Lower Triangular Matrix representation of 2nd moment fneq: fneq*c*c
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t c2TracePi = DXQY::c2*DXQY::traceLowTri(PiNeqLowTri);
    
    for (int q = 0; q < DXQY::nQ; ++q)
    {
      std::valarray<lbBase_t> ccLowTri(DXQY::nD*(DXQY::nD+1)/2);
      std::vector<int> cq = DXQY::c(q);
      int it=0;

      for(int i = 0; i < DXQY::nD; ++i){
	for(int j = 0; j <= i ; ++j){
	  ccLowTri[it]= cq[i]*cq[j];
	  it++;
	}
      }
      ret[q] = feq[q]
	+ DXQY::w[q]*DXQY::c4Inv0_5*(DXQY::contractionLowTri(PiNeqLowTri,ccLowTri)-c2TracePi);
      
      /*
      ret[q] = feq[q]
      + DXQY::w[q]*DXQY::c4Inv0_5*(DXQY::dot(DXQY::contractionLowTriVec(PiNeqLowTri,DXQY::c(q)),DXQY::c(q))-c2TracePi);
      */
      
      
    }
    return ret;
}

template <typename DXQY, typename T1, typename T2>
  inline std::valarray<lbBase_t> calcRegDist(const T1 &matLowTri, const lbBase_t& rho, const lbBase_t& u_sq, const T2 &cu)
/* calcOmegaBGK : sets the BGK-collision term in the lattice boltzmann equation
 *
 * matLowTri    : Lower Triangular Matrix representation of 2nd moment fneq: fneq*c*c
 * rho          : density
 * u_sq         : square of the velocity
 * cu           : array of scalar product of all lattice vectors and the velocity.
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    std::valarray<lbBase_t> feq = calcfeq<DXQY>(rho, u_sq, cu);
    for (int q = 0; q < DXQY::nQ; ++q)
    {
      ret[q] = feq[q]
      + DXQY::w[q]*DXQY::c4Inv0_5*(DXQY::dot(DXQY::contractionLowTriVec(matLowTri,DXQY::c(q)),DXQY::c(q))-DXQY::c2*DXQY::traceLowTri(matLowTri));
    }
    return ret;
}


#endif // LBUTILITIES_H
