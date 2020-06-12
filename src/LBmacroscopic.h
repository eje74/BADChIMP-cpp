#ifndef LBMACROSCOPIC_H
#define LBMACROSCOPIC_H

#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBfield.h"
#include "LBboundary.h"


// CALCULATION OF MACROSCOPIC VALUES
template <typename DXQY, typename T>  // Using typename T, makes the function work for both vectors and pointer
inline lbBase_t calcRho(const T &f)
/* calcRho : Calculates the macroscopic denisty at a node
 *
 * f   : pointer to the distribution at a node
 * rho : reference to the varabel where the density value at a node is stored
 */
{
     return DXQY::qSum(f);
}


template <typename DXQY, typename T1>
inline std::valarray<lbBase_t> calcVel(const T1 &f, const lbBase_t &rho)
/* calcVel : Calculates the macroscopic velocity at a node, using Guo's method.
 *
 * f     : pointer to the distribution at a node
 * rho   : reference to the varabel where the density value at a node is stored
 */
{
    return DXQY::qSumC(f) / rho;
}


template <typename DXQY, typename T1, typename T2>
inline std::valarray<lbBase_t> calcVel(const T1 &f, const lbBase_t &rho, const T2 &force)
/* calcVel : Calculates the macroscopic velocity at a node, using Guo's method.
 *
 * f     : pointer to the distribution at a node
 * rho   : reference to the varabel where the density value at a node is stored
 * force : pointer to the array where the force vector for a node is stored
 */
{
    return (DXQY::qSumC(f) + 0.5*force) / rho;
}



template <typename DXQY, typename T1, typename T2>
  inline std::valarray<lbBase_t> calcStrainRateTildeLowTri(const T1 &f, const lbBase_t &rho, const T2 &vel, const T2 &force, const lbBase_t &source)
/* calcVel : Calculates the LB strain rate, without 1/(2*rho*c2*tau), at a node, and returns 
 * it as an array with elements of a lower triangular matrix
 *
 * f     : pointer to the distribution at a node
 * rho   : reference to the varabel where the density value at a node is stored
 * force : pointer to the array where the force vector for a node is stored
 */
{
  std::valarray<lbBase_t> ret = DXQY::qSumCCLowTri(f);
  int it=0;
  for(int i = 0; i < DXQY::nD; ++i){
    for(int j = 0; j <= i ; ++j){
      ret[it]+= - (rho-0.5*source)*(DXQY::c2*DXQY::UnitMatrixLowTri[it] +vel[i]*vel[j]) + 0.5*(vel[i]*force[j]+vel[j]*force[i]);
      //ret[it]+= - (DXQY::qSumC(f))*(DXQY::c2*DXQY::UnitMatrixLowTri[it] +vel[i]*vel[j]) + 0.5*(vel[i]*force[j]+vel[j]*force[i]);
      it++;
    }
  }
  return -ret;
}

template <typename DXQY, typename T1, typename T2>
inline std::valarray<lbBase_t> calcShearRateTildeLowTri(const T1 &f, const lbBase_t &rho, const T2 &vel, const T2 &force, const lbBase_t &source)
/* calcVel : Calculates the LB shear rate, without 1/(2*rho*c2*tau), at a node, and returns 
 * it as an array with elements of a lower triangular matrix
 *
 * f     : pointer to the distribution at a node
 * rho   : reference to the varabel where the density value at a node is stored
 * force : pointer to the array where the force vector for a node is stored
 * source   : reference to the varabel where the mass source value at a node is stored
 */
{
  std::valarray<lbBase_t> ret = calcStrainRateTildeLowTri<DXQY>(f, rho, vel, force, source);
  lbBase_t ETrace = DXQY::traceLowTri(ret);
  
  for(int i = 0; i < DXQY::nD*(DXQY::nD+1)/2; ++i){
    ret -=  ETrace*DXQY::UnitMatrixLowTri[i]/DXQY::nD;
  }
  
  return ret;
}

template <typename DXQY, typename T1, typename T2>
inline std::valarray<lbBase_t> calcStrainRateTilde(const T1 &f, const lbBase_t &rho, const T2 &vel, const T2 &force, const lbBase_t &source)
/* calcVel : Calculates the LB strain rate without, 1/(2*rho*c2*tau), at a node.
 *
 * f     : pointer to the distribution at a node
 * rho   : reference to the varabel where the density value at a node is stored
 * force : pointer to the array where the force vector for a node is stored
 * source   : reference to the varabel where the mass source value at a node is stored
 */
{
  std::valarray<lbBase_t> ELowTri = calcStrainRateTildeLowTri<DXQY>(f, rho, vel, force, source);    
  std::valarray<lbBase_t> ret(DXQY::nD*DXQY::nD);
  
  int it=0;
  for(int i = 0; i < DXQY::nD; ++i){
    for(int j = 0; j <= i ; ++j){
      ret[i*DXQY::nD + j] = ELowTri[it];
      it++;
    }
  }  
  for(int i = 0; i < DXQY::nD; ++i)
    for(int j = i+1; j < DXQY::nD; ++j)
      ret[i*DXQY::nD + j] = ret[j*DXQY::nD + i];
  
  return ret;
}

template <typename DXQY, typename T1, typename T2>
inline std::valarray<lbBase_t> calcShearRateTilde(const T1 &f, const lbBase_t &rho, const T2 &vel, const T2 &force, const lbBase_t &source)
/* calcShearRateTilde : Calculates the LB shear rate, without 1/(2*rho*c2*tau), at a node.
 *
 * f     : pointer to the distribution at a node
 * rho   : reference to the varabel where the density value at a node is stored
 * force : pointer to the array where the force vector for a node is stored
 * source   : reference to the varabel where the mass source value at a node is stored
 */
{
  std::valarray<lbBase_t> SLowTri = calcShearRateTildeLowTri<DXQY>(f, rho, vel, force, source);    
  std::valarray<lbBase_t> ret(DXQY::nD*DXQY::nD);
  
  int it=0;
  for(int i = 0; i < DXQY::nD; ++i){
    for(int j = 0; j <= i ; ++j){
      ret[i*DXQY::nD + j] = SLowTri[it];
      it++;
    }
  }  
  for(int i = 0; i < DXQY::nD; ++i)
    for(int j = i+1; j < DXQY::nD; ++j)
      ret[i*DXQY::nD + j] = ret[j*DXQY::nD + i];
  
  return ret;
}


// FILL FIELDS
inline void setFieldToConst(const lbBase_t rhoConst, const int &fieldNo, ScalarField &rho)
/* sets all density values for field 'fieldNo' to a given density
 *
 * rhoConst : value copied to all node
 * fieldNo  : the number of the field to fill
 * rho      : reference to the global density object
 */
{
    for (int n = 0; n < rho.size(); ++n)
        rho(fieldNo, n) = rhoConst;
}

template <typename DXQY>
inline void setFieldToConst(const lbBase_t* velConst, const int &fieldNo, VectorField<DXQY> &vel)
/* sets all velocity values for field 'fieldNo' to a given vector
 *
 * velConst : pointer to vector that is copied to all node
 * fieldNo  : the number of the field to fill
 * vel      : reference to the global velocity object
 */

{
    for (int n = 0; n < vel.getNumNodes(); ++n)
        for (int d = 0; d < DXQY::nD; ++d)
            vel(fieldNo, d, n) = velConst[d];
}

template <typename DXQY>
inline void setFieldToConst(const lbBase_t* fConst, const int &fieldNo, LbField<DXQY> &f)
/* sets all lb distribution  values for field 'fieldNo' to a given distribution
 *
 * fConst  : pointer to distribution that is copied to all node
 * fieldNo : the number of the field to fill
 * f       : reference to the global lb distribution object
 */
{
    for (int n = 0; n < f.getNumNodes(); ++n)
        for (int q = 0; q < DXQY::nQ; ++q)
            f(fieldNo, q, n) = fConst[q];
}


template <typename DXDY>
inline void setFieldToConst(const Boundary<DXDY> &bnd, const lbBase_t rhoConst, const int &fieldNo, ScalarField &rho)
/* sets all density values for field 'fieldNo' to a given density
 *
 * rhoConst : value copied to all node
 * fieldNo  : the number of the field to fill
 * rho      : reference to the global density object
 */
{
    for (int n = 0; n < bnd.size(); ++n)
        rho(fieldNo, bnd.nodeNo(n)) = rhoConst;
}



// SET GLOBAL VALUES
inline void storeGlobalValues(const int fieldNo, const int nodeNo, const lbBase_t rhoLocal, ScalarField &rho)
/* storeGlobalValues : store the local density in the global density object
 *
 * fieldNo  : the current field number
 * nodeNo   : the current node
 * rhoLocal : the local density
 * rho      : reference to the global density object
 */
{
    rho(fieldNo, nodeNo) = rhoLocal;
}

template <typename DXQY>
inline void storeGlobalValues(const int fieldNo, const int nodeNo, const lbBase_t rhoLocal, const lbBase_t * velLocal,  ScalarField &rho, VectorField<DXQY> &vel)
/* storeGlobalValues : store the local density and local velocity in the global density and global velocity objects
 *
 * fieldNo  : the current field number
 * nodeNo   : the current node
 * rhoLocal : the local density
 * velLocal : pointer to the local velocity vector
 * rho      : reference to the global density object
 * vel      : reference to the global velocity object
 */
{
    rho(fieldNo, nodeNo) = rhoLocal;
    for (int d = 0; d < DXQY::nD; ++d) {
        vel(fieldNo, d, nodeNo) = velLocal[d];
    }
}


#endif // LBMACROSCOPIC_H
