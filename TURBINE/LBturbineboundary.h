#ifndef LBTURBINEBOUNDARY
#define LBTURBINEBOUNDARY

#include "../src/LBglobal.h"
#include "../src/LBboundary.h"
#include "../src/LBgrid.h"
#include "../src/LBnodes.h"
#include "../src/LBfield.h"
#include "LBturbineforce.h"



// ***************************************************************************************** INLET BOUNDARY
template<typename DXQY>
class InletBoundary : public Boundary<DXQY>
{
public:
    InletBoundary(const double omega_x, const VectorField<DXQY> &pos, const std::vector<int> &bndNodes,
                  const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
        :Boundary<DXQY>(bndNodes, nodes, grid), w_(omega_x), vel_(1, bndNodes.size()), force_(1, bndNodes.size()), dirNeig_(bndNodes.size()) {
        // Set the direction of the normal
        std::vector<int> norm {1, 0, 0};

        // Calculate the velocity and force due to the rotationg coordiante system
        // Set y and z velocity components. Assuming zero inertial velocity
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            int n = this->nodeNo(bndNo);
            std::valarray<lbBase_t> nodePos = pos(0, n);
            vel_(0, 0, bndNo) =  0;
            vel_(0, 1, bndNo) =  omega_x*pos(0, 2, n);
            vel_(0, 2, bndNo) = -omega_x*pos(0, 1, n);
            force_(0, bndNo) = rotatingForce(omega_x, 1.0, nodePos, vel_(0, bndNo));
        }
        // Find and set the neighbor direction
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            double val = -10;
            int qDir = -1;
            for (auto beta: this->beta(bndNo)) {
                int neigNo = grid.neighbor(beta, this->nodeNo(bndNo));
                if (!nodes.isFluidBoundary(neigNo)) {
                    double val_tmp = DXQY::dot(DXQY::c(beta), norm)/DXQY::cNorm[beta];
                    if (val < val_tmp) {
                        val = val_tmp;
                        qDir = beta;
                    }
                }
            }
            dirNeig_[bndNo] = qDir;
        }
    }
    void apply(const lbBase_t &tau, const std::valarray<lbBase_t> &bodyForce, const lbBase_t &rho, const VectorField<DXQY> &velFluid, const int &fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid)
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
            std::valarray<lbBase_t> velBnd = vel_(0, bndNo);
            int nodeNo = this->nodeNo(bndNo);
            velBnd[0] = 0.00; //velFluid(fieldNo, 0, grid.neighbor(dirNeig_[bndNo], nodeNo));
            // Calculate values used feq
            lbBase_t u2 = DXQY::dot(velBnd, velBnd);
            std::valarray<lbBase_t> cu = DXQY::cDotAll(velBnd);
            // Calculate the force
            std::valarray<lbBase_t> F = force_(0, bndNo);
            F += bodyForce;
            lbBase_t uF = DXQY::dot(velBnd, F);
            std::valarray<lbBase_t> cF = DXQY::cDotAll(F);
            // Collision terms
            auto feq = [&](int q) -> lbBase_t {return rho * DXQY::w[q]*(1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u2));};

            std::valarray<lbBase_t> OmegaF = calcDeltaOmegaF<DXQY>(tau, cu, uF, cF);

            // Need to define a lambda that evaluates the equilibirum function for a given direction
            // Need to calculate the forces
            // Need to define the force term
            for (auto beta: this->beta(bndNo)) {
                f(fieldNo, beta, grid.neighbor(beta, nodeNo)) = feq(beta) + OmegaF[beta];
            }

        }
    }

private:
    VectorField<DXQY> vel_;
    VectorField<DXQY> force_;
    std::vector<int> dirNeig_;
    lbBase_t w_;
};


// ***************************************************************************************** OUTLET BOUNDARY
template<typename DXQY>
class OutletBoundary : public Boundary<DXQY>
{
public:
    OutletBoundary(const double omega_x, const VectorField<DXQY> &pos, const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
        :Boundary<DXQY>(bndNodes, nodes, grid), pos_(1, bndNodes.size()), dirNeig_(bndNodes.size()) {
        // Set the direction of the normal
        std::vector<int> norm {-1, 0, 0};
        // rotation
        w_ = omega_x;
        // Set position for the boundary nodes
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            int n = this->nodeNo(bndNo);
            for (int d = 0; d < DXQY::nD; ++d)
                pos_(0, d, bndNo) = pos(0, d, n);
        }
        // Find and set the neighbor direction
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            double val = -10;
            int qDir = -1;
            for (auto beta: this->beta(bndNo)) {
                int neigNo = grid.neighbor(beta, this->nodeNo(bndNo));
                if (!nodes.isFluidBoundary(neigNo)) {
                    double val_tmp = DXQY::dot(DXQY::c(beta), norm)/DXQY::cNorm[beta];
                    if (val < val_tmp) {
                        val = val_tmp;
                        qDir = beta;
                    }
                }
            }
            dirNeig_[bndNo] = qDir;
        }
    }

    void apply(const lbBase_t &tau, const std::valarray<lbBase_t> &bodyForce, const lbBase_t &rho, const VectorField<DXQY> &velFluid,
               const int &fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid)
    /*
     *  - beta          : Unknown
     *  - beta_revers   : Known
     *
     *  - gamma         : Known
     *  - gamma_reverse : Known
     *
     *  - delta         : Unknown
     *  - delta_reverse : Unknown
     *
     *
     * For the inlet we will just copy the velocity from node in the upwind direction
     * or the node closes to the upwind direction
     */
    {
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            // Set velocity
            int nodeNo = this->nodeNo(bndNo);
            std::valarray<lbBase_t> velInit {0.0, 0.0, 0.0};
            std::valarray<lbBase_t> posNode = pos_(0, bndNo);

            std::valarray<lbBase_t> velBnd  =  velInert2Rot(w_, posNode , velInit);
            velBnd[0] = 0;
            velBnd[1] = w_*posNode[2];
            velBnd[2] = -w_*posNode[1];
                        
            // Calculate values used feq
            lbBase_t u2 = DXQY::dot(velBnd, velBnd);
            std::valarray<lbBase_t> cu = DXQY::cDotAll(velBnd);
            // Calculate the force

            std::valarray<lbBase_t> F = rotatingForce(w_, rho, posNode, velBnd);
            F += bodyForce;
            lbBase_t uF = DXQY::dot(velBnd, F);
            std::valarray<lbBase_t> cF = DXQY::cDotAll(F);
            // Collision terms
            auto feq = [&](int q) -> lbBase_t {return rho * DXQY::w[q]*(1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u2));};

            std::valarray<lbBase_t> OmegaF = calcDeltaOmegaF<DXQY>(tau, cu, uF, cF);

            // Need to define a lambda that evaluates the equilibirum function for a given direction
            // Need to calculate the forces
            // Need to define the force term
            for (auto beta: this->beta(bndNo)) {
                f(fieldNo, beta, grid.neighbor(beta, nodeNo)) = feq(beta) + OmegaF[beta];
            }
        }
    }

private:
    VectorField<DXQY> pos_;
    lbBase_t w_;
    std::vector<int> dirNeig_;
};




// ***************************************************************************************** FREEBOUNDARY BOUNDARY
template<typename DXQY>
class FreeBoundary : public Boundary<DXQY>
{
public:
    FreeBoundary(const double omega_x, const VectorField<DXQY> &pos, const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
        :w_(omega_x), Boundary<DXQY>(bndNodes, nodes, grid), pos_(1, bndNodes.size()), dirNeig_(bndNodes.size()), normal_(1, bndNodes.size()) {
        // Set the direction of the normal


        // Set position for the boundary nodes
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            int n = this->nodeNo(bndNo);
            std::valarray<lbBase_t> r(3);
            for (int d = 0; d < 3; ++d) {
                pos_(0, d, bndNo) = pos(0, d, n);
                r[d] = pos(0, d, n);
            }
            lbBase_t norm_inv = 1.0/sqrt(r[1]*r[1] + r[2]*r[2]);
            normal_(0, 0, bndNo) = 0;
            normal_(0, 1, bndNo) = -r[1]*norm_inv;
            normal_(0, 2, bndNo) = -r[2]*norm_inv;
        }


        // Find and set the neighbor direction
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            double val = -10;
            int qDir = -1;
            for (auto beta: this->beta(bndNo)) {
                int neigNo = grid.neighbor(beta, this->nodeNo(bndNo));
                if (!nodes.isFluidBoundary(neigNo)) {
                    double val_tmp = DXQY::dot(DXQY::c(beta), normal_(0, bndNo))/DXQY::cNorm[beta];
                    if (val < val_tmp) {
                        val = val_tmp;
                        qDir = beta;
                    }
                }
            }
            for (auto gamma: this->gamma(bndNo)) {
                int neigNo = grid.neighbor(gamma, this->nodeNo(bndNo));
                if (!nodes.isFluidBoundary(neigNo)) {
                    double val_tmp = DXQY::dot(DXQY::c(gamma), normal_(0, bndNo))/DXQY::cNorm[gamma];
                    if (val < val_tmp) {
                        val = val_tmp;
                        qDir = gamma;
                    }
                }
                int gammaRev = this->dirRev(gamma);
                neigNo = grid.neighbor(gammaRev, this->nodeNo(bndNo));
                if (!nodes.isFluidBoundary(neigNo)) {
                    double val_tmp = DXQY::dot(DXQY::c(gammaRev), normal_(0, bndNo))/DXQY::cNorm[gammaRev];
                    if (val < val_tmp) {
                        val = val_tmp;
                        qDir = gammaRev;
                    }
                }
            }
            dirNeig_[bndNo] = qDir;
        }
    }

    void apply(const lbBase_t &tau, const std::valarray<lbBase_t> &bodyForce, const lbBase_t &rho, const VectorField<DXQY> &velFluid,
               const int &fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid)
    /*
     *  - beta          : Unknown
     *  - beta_revers   : Known
     *
     *  - gamma         : Known
     *  - gamma_reverse : Known
     *
     *  - delta         : Unknown
     *  - delta_reverse : Unknown
     *
     *
     * For the inlet we will just copy the velocity from node in the upwind direction
     * or the node closes to the upwind direction
     */
    {
        for (int bndNo=0; bndNo < this->size(); ++bndNo) {
            // Set velocity
            if (dirNeig_[bndNo] != -1) {
                int nodeNo = this->nodeNo(bndNo);
                int neigNo = grid.neighbor(dirNeig_[bndNo], nodeNo);
                std::valarray<lbBase_t> velBnd = velFluid(0, neigNo);
                // velBnd = velRot2Inert(w_, pos_(0, bndNo), velBnd);
                velBnd -= DXQY::dot(normal_(0, bndNo), velBnd)*normal_(0, bndNo);
                // velBnd = velInert2Rot(w_, pos_(0, bndNo), velBnd);
                // Calculate values used feq
                lbBase_t u2 = DXQY::dot(velBnd, velBnd);
                std::valarray<lbBase_t> cu = DXQY::cDotAll(velBnd);
                // Calculate the force

                std::valarray<lbBase_t> F = rotatingForce(w_, rho, pos_(0, bndNo), velBnd);
                F += bodyForce;
                lbBase_t uF = DXQY::dot(velBnd, F);
                std::valarray<lbBase_t> cF = DXQY::cDotAll(F);
                // Collision terms
                auto feq = [&](int q) -> lbBase_t {return rho * DXQY::w[q]*(1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u2));};

                std::valarray<lbBase_t> OmegaF = calcDeltaOmegaF<DXQY>(tau, cu, uF, cF);

                // Need to define a lambda that evaluates the equilibirum function for a given direction
                // Need to calculate the forces
                // Need to define the force term
                for (auto beta: this->beta(bndNo)) {
                    f(fieldNo, beta, grid.neighbor(beta, nodeNo)) = feq(beta) + OmegaF[beta];
                }
                for (auto gamma: this->gamma(bndNo)) {
                    int alpha = gamma;
                    f(fieldNo, alpha, grid.neighbor(alpha, nodeNo)) = feq(alpha) + OmegaF[alpha];
                    alpha = this->dirRev(gamma);
                    f(fieldNo, alpha, grid.neighbor(alpha, nodeNo)) = feq(alpha) + OmegaF[alpha];
                }

            }

        }
    }

    inline int getDirNeig(const int &bndNo) {
        return dirNeig_[bndNo];
    }

private:
    VectorField<DXQY> pos_;
    VectorField<DXQY> normal_;
    std::vector<int> dirNeig_;
    lbBase_t w_;
};

#endif
