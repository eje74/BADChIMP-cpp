#ifndef LBUTILITIES_H
#define LBUTILITIES_H

#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBgrid.h"
//#include "Input.h"
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
inline std::valarray<lbBase_t> gradHigher(const ScalarField &sField, const int fieldNum, const int nodeNo, const Grid<DXQY> &grid)
{
    std::valarray<lbBase_t> scalarTmp(DXQY::nQ);
    std::valarray<lbBase_t> scalar2Tmp(DXQY::nQ);
    for (int q = 0; q < DXQY::nQ; ++q) {
        int neigNode = grid.neighbor(q, nodeNo);
	int nextNeigNode = grid.neighbor(q, neigNode);
        scalarTmp[q] = sField(fieldNum, neigNode);
	scalar2Tmp[q] = sField(fieldNum, nextNeigNode);
    }

    lbBase_t alpha = 1;//1/3.;
    
    //return (DXQY::grad(scalarTmp) + 2.0*DXQY::grad(scalar2Tmp))*0.25;
    //return (8.0*DXQY::grad(scalarTmp)-DXQY::grad(scalar2Tmp))/6.;
    //return 0.5*DXQY::grad(scalar2Tmp);

    //return (1-alpha)*DXQY::grad(scalarTmp) + alpha*0.5*(4*DXQY::grad(scalarTmp)-DXQY::grad(scalar2Tmp));
    return (1-alpha)*DXQY::grad(scalarTmp) + alpha*(8*DXQY::grad(scalarTmp)-DXQY::grad(scalar2Tmp))/6;
    //return 0.66666666666667*DXQY::grad(scalar2Tmp);
    //return (8.0*DXQY::grad(scalarTmp)- 1.0*DXQY::grad(scalar2Tmp))/6.;
    //return DXQY::grad(scalarTmp);
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
inline std::valarray<T> inputAsValarray(const Block& input)
{    std::vector<T> tmpVec = input;
     return std::valarray<T>(tmpVec.data(), tmpVec.size());
}

template <typename DXQY, typename T>
inline std::valarray<lbBase_t> calcfeq(const lbBase_t& rho, const lbBase_t& u_sq, const T &cu)
/* calcfeq : calculates the equilibrium distribution
 *
 * rho      : density
 * u_sq     : square of the velocity
 * cu       : array of scalar product of all lattice vectors and the velocity.
 * ret      : array of the equilibrium distribution in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = rho * DXQY::w[q]*(1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq) );
    }
    return ret;
}

template <typename DXQY>
inline std::valarray<lbBase_t> calcfeq_TEST(const lbBase_t& Gamma0, const lbBase_t& GammaNonZero, const lbBase_t& rho, const lbBase_t& u_sq, const std::valarray<lbBase_t> &cu)
/* calcfeq_TEST : calculates the equilibrium distribution
 *
 * Gamma0       : static mass fraction
 * GammaNonZero : compensation for static mass fraction
 * rho          : density
 * u_sq         : square of the velocity
 * cu           : array of scalar product of all lattice vectors and the velocity.
 * ret          : array of the equilibrium distribution in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    
    for (int q = 0; q < DXQY::nQNonZero_; ++q)
    {
      ret[q] = rho*DXQY::w[q]*(GammaNonZero + (DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq)) ) ;
    }
    ret[DXQY::nQNonZero_] = rho * DXQY::w[DXQY::nQNonZero_]*(Gamma0-0.5*DXQY::c2Inv*u_sq);
    return ret;
}



template <typename DXQY, typename T>
inline std::valarray<lbBase_t> calcfeqHANS(const lbBase_t& rho, const lbBase_t& u_sq, const T &cu, const lbBase_t& momCoefHANS)
/* calcfeq : calculates the equilibrium distribution
 *
 * rho         : density
 * u_sq        : square of the velocity
 * cu          : array of scalar product of all lattice vectors and the velocity.
 * momCoefHANS : scalar proportionality coefficient, weigthing the non-linear momentum contribution
 * ret         : array of the equilibrium distribution in each lattice direction
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    for (int q = 0; q < DXQY::nQ; ++q)
    {
        ret[q] = rho * DXQY::w[q]*(1.0 + DXQY::c2Inv*cu[q] + momCoefHANS*DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq) );
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

template <typename DXQY, typename T1, typename T2, typename T3>
inline std::valarray<lbBase_t> calcRegDist(const T1 &feq, const lbBase_t& MNeq, const T2 &MiNeq, const T3 &PiNeqLowTri)
/* calcOmegaBGK : sets the BGK-collision term in the lattice boltzmann equation
 *
 * feq          : precalculated array of local equilibriums distributions.
 * qSrc         : Mass source
 * cF           : array of scalar product of all lattice vectors and the Force.
 * cu           : array of scalar product of all lattice vectors and the velocity.
 * PiNeqLowTri    : Lower Triangular Matrix representation of 2nd moment fneq: fneq*c*c
 */
{
    std::valarray<lbBase_t> ret(DXQY::nQ);
    lbBase_t c2TracePi = DXQY::c2*DXQY::traceLowTri(PiNeqLowTri);

    const auto cMiNeq = DXQY::cDotAll(MiNeq);
    
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
	+ DXQY::w[q]*(MNeq + DXQY::c2Inv*cMiNeq[q] + DXQY::c4Inv0_5*(DXQY::contractionLowTri(PiNeqLowTri,ccLowTri)-c2TracePi));
      
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
