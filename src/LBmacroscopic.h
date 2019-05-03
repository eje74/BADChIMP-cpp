#ifndef LBMACROSCOPIC_H
#define LBMACROSCOPIC_H

#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBfield.h"
#include "LBboundary.h"


// CALCULATION OF MACROSCOPIC VALUES
template <typename DXQY>
inline void calcRho(const lbBase_t* f, lbBase_t &rho)
/* calcRho : Calculates the macroscopic denisty at a node
 *
 * f   : pointer to the distributioin at a node
 * rho : reference to the varabel where the density value at a node is stored
 */
{
     DXQY::qSum(f, rho);
}

template <typename DXQY>
inline void calcVel(const lbBase_t* f, const lbBase_t &rho, lbBase_t* vel, const lbBase_t* force)
/* calcVel : Calculates the macroscopic velocity at a node, using Guo's method.
 *
 * f     : pointer to the distributioin at a node
 * rho   : reference to the varabel where the density value at a node is stored
 * vel   : pointer to the array where the velocity vector for a node is stored
 * force : pointer to the array where the force vector for a node is stored
 */
{
    DXQY::qSumC(f, vel);
    for (int d = 0; d < DXQY::nD; ++d) {
      vel[d] = (vel[d] + 0.5 * force[d]) /rho;
    }
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
