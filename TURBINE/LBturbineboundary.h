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
    InletBoundary(const double omega_x, const VectorField<DXQY> &pos, const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
        :Boundary<DXQY>(bndNodes, nodes, grid), vel_(1, bndNodes.size()) {
        // Set the direction of the normal
        std::vector<int> norm {1, 0, 0};
        qDirNormal_ =  DXQY::c2q(norm);
        // Set y and z velocity components. Assuming zero inertial velocity
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            int n = this->nodeNo(bndNo);
            vel_(0, 0, bndNo) =  0;
            vel_(0, 1, bndNo) =  omega_x*pos(0, 2, n);
            vel_(0, 2, bndNo) = -omega_x*pos(0, 1, n);
        }

    }

    void apply(const VectorField<DXQY> velFluid, const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid)
    /*
     *  - beta          : Unknown
     *  - beta_revers   : Known
     *
     *  - gamma         : Known
     *  - gamma_reverse : Known
     *
     *  - delta         : Unknown
     *  - delta_reverse : Unknown
     */
    {
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            // Set velocity
            std::valarray<lbBase_t> vel = vel_(0, bndNo);
            int nodeNo = this->nodeNo(bndNo);
            vel[0] = velFluid(fieldNo, 0, nodeNo);
            // Calculate values used feq
            lbBase_t u2 = DXQY::dot(vel, vel);
            // Need to define a lambda that evaluates the equilibirum function for a given direction
            // Need to calculate the forces
            // Need to define the force term
            for (auto beta: this->beta(bndNo)) {
                int alpha = this->dirRev(beta);
 //               f(fieldNo, qSlipDir_[alpha], nodeNoFluid) = f(fieldNo, alpha, nodeNo);
            }
            for (auto gamma: this->gamma(bndNo)) {
                int alpha = gamma;
//                f(fieldNo, qSlipDir_[alpha], nodeNoFluid) = f(fieldNo, alpha, nodeNo);
                alpha = this->dirRev(gamma);
//                f(fieldNo, qSlipDir_[alpha], nodeNoFluid) = f(fieldNo, alpha, nodeNo);
            }
        }
    }
 
        // cDot(const int qDir, const T* rightVec);
        // rho * DXQY::w[q]*(1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq) )

private:
    VectorField<DXQY> vel_;
    int qDirNormal_;
};

#endif
