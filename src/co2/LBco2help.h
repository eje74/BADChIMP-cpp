#ifndef LBCO2HELP_H
#define LBCO2HELP_H

#include <string>
#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBvtk.h"
#include "../lbsolver/LBnodes.h"



template<typename DXQY>
class CGAtributes
{
public:
    template<typename T, typename U>
    CGAtributes(const int nFluidFields, const int nodeNo, const std::valarray<lbBase_t> cNormInv, const std::valarray<lbBase_t> &Gamma0, const T &rhoRelNode, const U &rhoRel, const Grid<DXQY> &grid);
    const int lowerTriangularSize_;
    const std::valarray<lbBase_t> GammaNonZero_;
    VectorField<DXQY> F_; //(1, (nFluidFields*(nFluidFields-1))/2);
    ScalarField FSquare_;//(1, F_.size());
    ScalarField FNorm_;//(1, F_.size());
    LbField<DXQY> cDotFRC_;//(1, F_.size());
    LbField<DXQY> cosPhi_;//(1, F_.size());
    VectorField<DXQY> gradNode_;//(1, nFluidFields);   
    //total effect of modified compressibility
    lbBase_t Gamma0TotNode_;
    lbBase_t GammaNonZeroTotNode_;
};

template<typename DXQY>
template<typename T, typename U>
CGAtributes<DXQY>::CGAtributes(const int nFluidFields, const int nodeNo, const std::valarray<lbBase_t> cNormInv, const std::valarray<lbBase_t> &Gamma0, const T &rhoRelNode, const U &rhoRel, const Grid<DXQY> &grid):
lowerTriangularSize_((nFluidFields*(nFluidFields-1))/2), GammaNonZero_((1-DXQY::w0*Gamma0)/(1-DXQY::w0)), F_(1, lowerTriangularSize_), FSquare_(1, lowerTriangularSize_),
FNorm_(1, lowerTriangularSize_), cDotFRC_(1, lowerTriangularSize_), cosPhi_(1, lowerTriangularSize_), gradNode_(1, nFluidFields)
{
    Gamma0TotNode_ = 0;
    GammaNonZeroTotNode_ = 0;
    int cnt = 0;
  
    for (int fieldNo_k=0; fieldNo_k<nFluidFields; ++fieldNo_k) {
        Gamma0TotNode_ += rhoRelNode(0, fieldNo_k)*Gamma0[fieldNo_k];
        GammaNonZeroTotNode_ += rhoRelNode(0, fieldNo_k)*GammaNonZero_[fieldNo_k];
        gradNode_.set(0, fieldNo_k) = grad<DXQY>(rhoRel, fieldNo_k, nodeNo, grid);
        for (int fieldNo_l = 0; fieldNo_l < fieldNo_k; ++fieldNo_l) {
            F_.set(0, cnt) = rhoRelNode(0, fieldNo_l)*gradNode_(0, fieldNo_k) - rhoRelNode(0, fieldNo_k)*gradNode_(0, fieldNo_l);
            cDotFRC_.set(0, cnt) = DXQY::cDotAll(F_(0,cnt));
            FSquare_(0, cnt) = DXQY::dot(F_(0, cnt), F_(0, cnt));
            FNorm_(0,cnt) = sqrt(FSquare_(0, cnt));
            if (abs(FNorm_(0, cnt)) < lbBaseEps)
                FNorm_(0, cnt) = lbBaseEps;
            cosPhi_.set(0, cnt) = cDotFRC_(0, cnt)*cNormInv/FNorm_(0,cnt);
            cnt++;                    
        }
    }    
};


template <typename DXQY>
void setScalarAttribute(ScalarField &field, const std::string &attributeName, LBvtk<DXQY> &vtklb) {
    for (int fieldNo=0; fieldNo < field.num_fields(); fieldNo++) {
        vtklb.toAttribute(attributeName + std::to_string(fieldNo));
        for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
            field(fieldNo, n) = vtklb.template getScalarAttribute<lbBase_t>(); //vtklb.getScalarAttribute<lbBase_t>();
        }
    }
}

template <typename DXQY>
void setScalarAttributeWall(ScalarField &field, const std::string &attributeName, LBvtk<DXQY> &vtklb, const Nodes<DXQY> &nodes) {
    for (int fieldNo=0; fieldNo < field.num_fields(); fieldNo++) {
        vtklb.toAttribute(attributeName + std::to_string(fieldNo));
        for (int n=vtklb.beginNodeNo(); n < vtklb.endNodeNo(); ++n) {
            const auto val = vtklb.template getScalarAttribute<lbBase_t>();
            if (nodes.isSolidBoundary(n)) {
                field(fieldNo, n) = val;
            }
        }
    }
}

template<typename T>
void normelizeScalarField(ScalarField &field, T nodeList)
{
    const int nFields = field.num_fields();
    for (auto nodeNo: nodeList) {
        lbBase_t tmp = 0;
        for (int fieldNo=0; fieldNo < nFields; ++fieldNo) {
            tmp += field(fieldNo, nodeNo);
        }
        for (int fieldNo=0; fieldNo < nFields; ++fieldNo) {
            field(fieldNo, nodeNo) /= tmp + lbBaseEps;
        }
    }
}

template<typename T, typename DXQY>
void calcDensityFields(ScalarField &rho, ScalarField &rhoRel, ScalarField &rhoTot, 
                       const T &nodeList, const LbField<DXQY> &f)
{
    const auto numFields = f.num_fields();
    for (auto nodeNo: nodeList) {
        rhoTot(0, nodeNo) = 0;
        for (int fieldNo=0; fieldNo < numFields; ++fieldNo) {
            const auto fNode = f(fieldNo, nodeNo);
            rho(fieldNo, nodeNo) = calcRho<DXQY>(fNode);
            rhoTot(0, nodeNo) += rho(fieldNo, nodeNo);
        }
        for (int fieldNo=0; fieldNo < numFields; ++fieldNo) {
            rhoRel(fieldNo, nodeNo) = rho(fieldNo, nodeNo)/rhoTot(0, nodeNo);
        }
    }
}



#endif