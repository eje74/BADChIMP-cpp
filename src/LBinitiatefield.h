#ifndef LBINITIATEFIELD_H
#define LBINITIATEFIELD_H

#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBfield.h"
//#include "Field.h"

template <typename DXQY>
void initiateLbField(const int fieldNo, const std::vector<int> &bulk, const ScalarField &rho, const VectorField<DXQY> &vel, LbField<DXQY> &f)
/* initiateLbField : sets the lb distributions of the given field, given by fieldNo,
 *  to the equilibirum distribution with denisty and velocity given by the macroscopic
 *  fields rho and vel.
 *
 * fieldNo : The field number of the distribution that is to be initiated
 * bulk    : reference to a vector of bulk node numbers . Supplies the node numbers (tags)
 * rho     : reference to the density field object
 * vel     : reference to the velocity field object
 * f       : reference to the velocity lb distribution field object
 */
{
    for (auto nodeNo: bulk) {
        lbBase_t cu[DXQY::nQ], uu;
        DXQY::cDotAll(&vel(fieldNo, 0, nodeNo), cu);
        uu = DXQY::dot(&vel(fieldNo, 0, nodeNo), &vel(fieldNo, 0, nodeNo));
        for (int q = 0; q < DXQY::nQ; q++) { // Set the lb field to its equlibrium distribution
            f(fieldNo, q, nodeNo) = DXQY::w[q] * rho(fieldNo, nodeNo) * (1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*uu));
        }
    }
}


template <typename DXQY>
void initiateLbField(const int lbFieldNo, const int rhoFieldNo, const int velFieldNo,
                     const std::vector<int> &bulk, const ScalarField &rho, const VectorField<DXQY> &vel, LbField<DXQY> &f)
/* initiateLbField : sets the lb distributions of the given field, given by fieldNo,
 *  to the equilibirum distribution with denisty and velocity given by the macroscopic
 *  fields rho and vel. Here we can also choose which density and velocity fields to
 *  use.
 *
 * lbFieldNo : The field number of the distribution that is to be initiated
 * rhoFieldNo : The field number of the density field that is used
 * velFieldNo : The field number of the velocity field that is used
 * bulk    : reference to the bulk object. Supplies the node numbers (tags)
 * rho     : reference to the density field object
 * vel     : reference to the velocity field object
 * f       : reference to the velocity lb distribution field object
 */
{
    for (auto nodeNo: bulk) {
        std::valarray<lbBase_t> cu = DXQY::cDotAll(&vel(velFieldNo, 0, nodeNo));
        lbBase_t uu = DXQY::dot(&vel(velFieldNo, 0, nodeNo), &vel(velFieldNo, 0, nodeNo));
        for (int q = 0; q < DXQY::nQ; q++) { // Set the lb field to its equlibrium distribution
            f(lbFieldNo, q, nodeNo) = DXQY::w[q] * rho(rhoFieldNo, nodeNo) * (1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*uu));
        }
    }
}


#endif // LBINITIATEFIELD_H
