#ifndef LBFREESLIPSOLID_H
#define LBFREESLIPSOLID_H

#include "LBglobal.h"
#include "LBboundary.h"
#include "LBgrid.h"
#include "LBnodes.h"
#include "LBfield.h"

template<typename DXQY>
class SolidFreeSlip
{
public:
    SolidFreeSlip(const std::vector<int> &outwardNormal, const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid);
    void apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid);
    inline int nodeNo(const int bndNo) const {return boundaryNodes_[bndNo];}
    inline int size() const {return boundaryNodes_.size();}
    inline std::vector<int> alpha(const int bndNo) const;
    inline int nAlpha(const int bndNo) const {return nAlphasPerNode_[bndNo];}
private:
    int qFluid_; // Direction of the fluid node x_fluid = x + c_qFluid.
    std::vector<int> qSlipDir_;  // used as f_qSlipDir[q](x_fluid) = f_q(x)  
    std::vector<int> boundaryNodes_;  // List of boundary nodes
    std::vector<int> nAlphasPerNode_;
    std::vector<int> alphasPerNode_;
    std::vector<int> beginAlphaBlock_;
};


template<typename DXQY>
SolidFreeSlip<DXQY>::SolidFreeSlip(const std::vector<int> &outwardNormal, const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
{
    std::vector<int> vec(DXQY::nD);
    // Find the inward points normal
    for (int d=0; d < DXQY::nD; ++d) 
        vec[d] = -outwardNormal[d];
    // Grid direction for the inward pointing normal 
    qFluid_ = DXQY::c2q(vec);
    // Find the reflected direction. Set qSlipDir_
    qSlipDir_.reserve(DXQY::nQ); 
    int normal[DXQY::nD]; // Hack, cDot sould take std::vector as argument
    for (int d=0; d<DXQY::nD; ++d) normal[d] = outwardNormal[d];
    for (int q=0; q<DXQY::nQ; ++q) {
        for (int d=0; d<DXQY::nD; ++d) {
            vec[d] = DXQY::c(q, d) - 2*DXQY::cDot(q, normal)*outwardNormal[d];
        }
        qSlipDir_[q] = DXQY::c2q(vec);
        if (qSlipDir_[q] == -1) {
            std::cout << "ERROR in FreeSlipSolid constructor: Reflected direction not found" << std::endl;
            exit(1);
        }
    }
    
    // Setup boundary nodes with directions
    
    // FOR ALL BOUNDARY NODES
    // CHECK normal dot c_alpha > 0
    // CHECK x - c_alpha isFluid
    int cnt = 0;
    for (const auto &node : bndNodes) {
        bool dirAdded = false;
        int nDirAdded = 0;
        for (int q = 0; q < DXQY::nQ; ++q) {
            if ( DXQY::cDot(q, normal)>0 ) { // Crosses the boundary
                auto neig = grid.neighbor(DXQY::reverseDirection(q), node);
                if ( nodes.isFluid(neig) ) { // distribution is propagated from a fluid node
                    alphasPerNode_.push_back(q);
                    dirAdded = true;
                    nDirAdded += 1;
                }
            }
        }
        
        if (dirAdded) {
            boundaryNodes_.push_back(node);
            nAlphasPerNode_.push_back(nDirAdded);
            cnt += 1;        
        }
    }
    beginAlphaBlock_.reserve(nAlphasPerNode_.size());
    beginAlphaBlock_.push_back(0);
    for (int i = 1; i < nAlphasPerNode_.size(); ++i)
        beginAlphaBlock_.push_back(beginAlphaBlock_[i-1] + nAlphasPerNode_[i-1]);
}

template<typename DXQY>
inline std::vector<int> SolidFreeSlip<DXQY>::alpha(const int bndNo) const
{
    auto ptrBegin = alphasPerNode_.data() + beginAlphaBlock_[bndNo];
    return std::vector<int>(ptrBegin, ptrBegin + nAlphasPerNode_[bndNo]);    
}

template<typename DXQY>
void SolidFreeSlip<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY>& grid)
{
    for (int bndNo = 0; bndNo < size(); ++bndNo) {
        int node = nodeNo(bndNo);
        int nodeFluid = grid.neighbor(qFluid_, node);
        for (const auto &a: alpha(bndNo)) {
            f(fieldNo, qSlipDir_[a], nodeFluid) = f(fieldNo, a, node);
        }
    }
}


template<typename DXQY>
class SolidFreeSlipOld : public Boundary<DXQY>
// NOTE: Here we assume that the outwardNormal is one of the cartesian directions
{
public:
    SolidFreeSlipOld(const std::vector<int> &outwardNormal, const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid);
    void apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid);
private:
    int qFluid_; // Direction of the fluid node x_fluid = x + c_qFluid.
    std::vector<int> qSlipDir_;  // used as f_qSlipDir[q](x_fluid) = f_q(x)
};


template<typename DXQY>
SolidFreeSlipOld<DXQY>::SolidFreeSlipOld(const std::vector<int> &outwardNormal, const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
:Boundary<DXQY>(bndNodes, nodes, grid)
{
    std::vector<int> vec(DXQY::nD);
    // Find the inward points normal
    for (int d=0; d < DXQY::nD; ++d) 
        vec[d] = -outwardNormal[d];
    // Grid direction for the inward pointing normal 
    qFluid_ = DXQY::c2q(vec);
    // Find the reflected direction. Set qSlipDir_
    qSlipDir_.reserve(DXQY::nQ); 
    int nVec[DXQY::nD]; // Hack, cDot sould take std::vector as argument
    for (int d=0; d<DXQY::nD; ++d) nVec[d] = outwardNormal[d];
    for (int q=0; q<DXQY::nQ; ++q) {
        for (int d=0; d<DXQY::nD; ++d) {
            vec[d] = DXQY::c(q, d) - 2*DXQY::cDot(q, nVec)*outwardNormal[d];
        }
        qSlipDir_[q] = DXQY::c2q(vec);
        if (qSlipDir_[q] == -1) {
            std::cout << "ERROR in FreeSlipSolid constructor: Reflected direction not found" << std::endl;
            exit(1);
        }
    }
}


template<typename DXQY>
void SolidFreeSlipOld<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid)
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
        int nodeNo = this->nodeNo(bndNo);
        int nodeNoFluid = grid.neighbor(qFluid_, nodeNo);
        for (auto beta: this->beta(bndNo)) {
            int alpha = this->dirRev(beta);
            f(fieldNo, qSlipDir_[alpha], nodeNoFluid) = f(fieldNo, alpha, nodeNo);            
        }
        for (auto gamma: this->gamma(bndNo)) {
            int alpha = gamma;
            f(fieldNo, qSlipDir_[alpha], nodeNoFluid) = f(fieldNo, alpha, nodeNo);
            alpha = this->dirRev(gamma);
            f(fieldNo, qSlipDir_[alpha], nodeNoFluid) = f(fieldNo, alpha, nodeNo);            
        }
    }    
}

#endif