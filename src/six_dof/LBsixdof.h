#ifndef LBSIXDOF_H
#define LBSIXDOF_H

#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBvtk.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBboundary.h"

#include <array>

//=====================================================================================
//
//                          W E T  N O D E  B O U N D A R Y
//
//=====================================================================================
template<typename DXQY>
class WetNodeBoundary
{
public:
    WetNodeBoundary(LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
    void applyBoundaryCondition(const int fieldNo, LbField<DXQY> &f, const VectorField<DXQY> & force, const Grid<DXQY> &grid) const;
private:
    template<typename T>
    std::vector<T> readScalarValues(const std::string attributeName, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
    template<typename T>
    std::vector<std::valarray<T>> readVectorValues(const std::string attributeName, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);

    const std::valarray<lbBase_t> w_;
    int numBoundaries_;
    struct DiffusionBoundaryNodes
    {
        std::vector<std::valarray<lbBase_t>> normals;
        std::vector<lbBase_t> qs;
        std::vector<int> indicators;
        std::vector<int> neighbors;
        Boundary<DXQY> bnd;
    } wallBoundary_, pressureBoundary_, wallPressureBoundary_;
};


//                               WetNodeBoundary
//----------------------------------------------------------------------------------- WetNodeBoundary
template<typename DXQY>
WetNodeBoundary<DXQY>::WetNodeBoundary(LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
:w_(DXQY::w, DXQY::nQ)
//-----------------------------------------------------------------------------------
{
    std::vector<int> wallBoundaryNodes;
    std::vector<int> pressureBoundaryNodes;
    std::vector<int> wallPressureBoundaryNodes;

    // Read pressure_boundary
    int numBnd = 0;
    vtklb.toAttribute("pressure_boundary");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        const auto pInd = vtklb.template getScalarAttribute<int>();
        numBnd = std::max(numBnd, pInd);
        if ( nodes.isFluidBoundary(nodeNo) ) {
            auto hasSolidNeighbors = [&nodeNo, &nodes, &grid]() -> bool {
                for (auto neighNo: grid.neighbor(nodeNo)) 
                    if ( nodes.isBulkSolid(neighNo) || nodes.isSolidBoundary(neighNo) ) 
                        return true;
                return false;
            };

            if (hasSolidNeighbors()) {
                if (pInd == 0) {
                    wallBoundaryNodes.push_back(nodeNo);
                 } else {
                    wallPressureBoundaryNodes.push_back(nodeNo);
                    wallPressureBoundary_.indicators.push_back(pInd);
                    // pressureBoundaryNodes.push_back(nodeNo);
                    // pressureBoundary_.indicators.push_back(pInd); 
                }
            } else {
                if (pInd > 0) {
                    pressureBoundaryNodes.push_back(nodeNo);
                    pressureBoundary_.indicators.push_back(pInd);
                } else {
                    std::cout << "Node " << nodeNo << " on processor rank " << nodes.getRank(nodeNo);
                    std::cout << " is not a wall node and is not assigned a pressure!" << std::endl;
                    exit(1);
                }
            }
        }  
    }
    // Setup the different boundary nodes.
    wallBoundary_.bnd = Boundary<DXQY>(wallBoundaryNodes, nodes, grid);
    pressureBoundary_.bnd =  Boundary<DXQY>(pressureBoundaryNodes, nodes, grid);
    wallPressureBoundary_.bnd = Boundary<DXQY>(wallPressureBoundaryNodes, nodes, grid);         
    numBoundaries_ = numBnd;

    auto qvalues = readScalarValues<lbBase_t>("q", vtklb, nodes, grid);
    auto normalvec = readVectorValues<lbBase_t>("normal", vtklb, nodes, grid);
    auto neighvec = readVectorValues<int>("neighbor", vtklb, nodes, grid);

    auto fillBoundaryNodes = [&qvalues, &normalvec, &neighvec](DiffusionBoundaryNodes &diffBnd) 
    {
        for (int n = 0; n < diffBnd.bnd.size(); ++n) {
            const int nodeNo = diffBnd.bnd.nodeNo(n);
            diffBnd.qs.push_back(qvalues[nodeNo]);
            diffBnd.normals.push_back(normalvec[nodeNo]);
            std::vector<int> nvec(DXQY::nD);
            for (int d=0; d < DXQY::nD; ++d)  nvec[d] = neighvec[nodeNo][d];
            diffBnd.neighbors.push_back(DXQY::c2q(nvec));
        }        
    };

    fillBoundaryNodes(wallBoundary_);
    fillBoundaryNodes(pressureBoundary_);
    fillBoundaryNodes(wallPressureBoundary_);
}

//                               WetNodeBoundary
//----------------------------------------------------------------------------------- applyBoundaryCondition
template<typename DXQY>
void WetNodeBoundary<DXQY>::applyBoundaryCondition(const int fieldNo, LbField<DXQY> &f, const VectorField<DXQY> & force, const Grid<DXQY> &grid) const
//-----------------------------------------------------------------------------------
{
    //----------------------------------------------------------------------------------- wall boundary
    for (int bndNo = 0; bndNo < wallBoundary_.bnd.size(); ++bndNo) {
        const int nodeNo = wallBoundary_.bnd.nodeNo(bndNo);
        const int alphaNeig = wallBoundary_.neighbors[bndNo];
        const int nodeNeig = grid.neighbor(alphaNeig, nodeNo);
        const auto fNeig = f(fieldNo, nodeNeig);

        // Calculate first moment
        const auto qNode = wallBoundary_.qs[bndNo];
        const std::valarray<lbBase_t> uNeig = DXQY::qSumC(fNeig) + 0.5*force(fieldNo, nodeNeig);
        const auto uNode = qNode*uNeig/(1 +  qNode);

        // Calculate zeroths moment
        lbBase_t lhs = 1.0;
        lbBase_t rhs = f(fieldNo, DXQY::nQNonZero_, nodeNo);

        for (auto &gamma: wallBoundary_.bnd.gamma(bndNo)) {
            rhs += f(fieldNo, gamma, nodeNo);
            rhs += f(fieldNo, DXQY::reverseDirection(gamma), nodeNo);
        }

        const std::valarray<lbBase_t> forceNode = force(fieldNo, nodeNo);
        for (auto &beta: wallBoundary_.bnd.beta(bndNo)) {
            const auto betaHat = DXQY::reverseDirection(beta);
            const auto c_beta = DXQY::c(beta);
            rhs += 2*f(fieldNo, betaHat, nodeNo);
            rhs -= w_[beta]*DXQY::c2Inv*DXQY::dot(c_beta, forceNode);
            lhs -= 2*w_[beta]*DXQY::c2Inv*DXQY::dot(c_beta, uNode);
        }

        const auto rhoNeig = DXQY::qSum(fNeig);
        const auto ccfNeig = DXQY::qSumCCLowTri(fNeig);
        const auto piTrace = DXQY::traceLowTri(ccfNeig);
        const auto uu = DXQY::dot(uNeig, uNeig);
        const auto c2 = DXQY::c2;

        std::valarray<lbBase_t> QPi_neq(DXQY::nQ);
        for (auto &delta: wallBoundary_.bnd.delta(bndNo)) {
            const auto c_delta = DXQY::c(delta);
            const auto ccpi = DXQY::dot(DXQY::contractionLowTriVec(ccfNeig, c_delta), c_delta);
            const auto cc = DXQY::dot(c_delta, c_delta);
            const auto cu = DXQY::dot(c_delta, uNeig);
            QPi_neq[delta] = ccpi-c2*piTrace - rhoNeig*c2*(cc - c2*DXQY::nD) - rhoNeig*(cu*cu - c2*uu);

            rhs += w_[delta]*DXQY::c4Inv*QPi_neq[delta];
            lhs -= 2*w_[delta];
            lhs -=  w_[delta]*DXQY::c4Inv*(cu*cu - c2*uu);
        }
        const auto rhoNode = rhs/lhs;

        for (auto &beta: wallBoundary_.bnd.beta(bndNo)) {           
            const auto betaHat = DXQY::reverseDirection(beta);
            const auto c_beta = DXQY::c(beta);
            const auto cu = DXQY::dot(c_beta, uNode);
            const auto cF = DXQY::dot(c_beta, forceNode);

            f(fieldNo, beta, nodeNo) = f(fieldNo, betaHat, nodeNo) + w_[beta]*DXQY::c2Inv*(2*rhoNode*cu - cF);         
        }

        const auto uuNode = DXQY::dot(uNode, uNode);

        for (auto &delta: wallBoundary_.bnd.delta(bndNo)) {
            const auto deltaHat = DXQY::reverseDirection(delta);
            const auto c_delta = DXQY::c(delta);
            const auto cc = DXQY::dot(c_delta, c_delta);
            const auto cu = DXQY::dot(c_delta, uNode);
            const auto cF = DXQY::dot(c_delta, forceNode);
            
            f(fieldNo, delta, nodeNo) = w_[delta]*rhoNode*(1 + DXQY::c2Inv*cu + DXQY::c4Inv0_5*(cu*cu - c2*uuNode) );
            f(fieldNo, delta, nodeNo) += w_[delta]*DXQY::c4Inv0_5*QPi_neq[delta] - w_[delta]*0.5*DXQY::c2Inv*cF;

            f(fieldNo, deltaHat, nodeNo) = w_[delta]*rhoNode*(1 - DXQY::c2Inv*cu + DXQY::c4Inv0_5*(cu*cu - c2*uuNode) );
            f(fieldNo, deltaHat, nodeNo) += w_[delta]*DXQY::c4Inv0_5*QPi_neq[delta] + w_[delta]*0.5*DXQY::c2Inv*cF;
        }     
    }
    //----------------------------------------------------------------------------------- wall pressure boundary
    for (int bndNo = 0; bndNo < wallPressureBoundary_.bnd.size(); ++bndNo) {
        const int nodeNo = wallPressureBoundary_.bnd.nodeNo(bndNo);
        const int alphaNeig = wallPressureBoundary_.neighbors[bndNo];
        const int nodeNeig = grid.neighbor(alphaNeig, nodeNo);
        const auto fNeig = f(fieldNo, nodeNeig);

        // Calculate first moment
        const auto qNode = wallPressureBoundary_.qs[bndNo];
        const std::valarray<lbBase_t> uNeig = DXQY::qSumC(fNeig) + 0.5*force(fieldNo, nodeNeig);
        std::valarray<lbBase_t>  uNode = qNode*uNeig/(1 +  qNode);
        const std::valarray<lbBase_t> forceNode = force(fieldNo, nodeNo);

        const lbBase_t rhoNode = 1.0;

        lbBase_t wUnknown = 0.0;
        std::valarray<lbBase_t> uUnknown(DXQY::nD);
        uUnknown = 0;

        for (auto &beta: wallPressureBoundary_.bnd.beta(bndNo)) {           
            const auto betaHat = DXQY::reverseDirection(beta);
            const auto c_beta = DXQY::c(beta);
            const auto cu = DXQY::dot(c_beta, uNode);
            const auto cF = DXQY::dot(c_beta, forceNode);

            f(fieldNo, beta, nodeNo) = f(fieldNo, betaHat, nodeNo) + w_[beta]*DXQY::c2Inv*(2*rhoNode*cu - cF);         
            wUnknown += w_[beta];
        }

        const auto rhoNeig = DXQY::qSum(fNeig);
        const auto ccfNeig = DXQY::qSumCCLowTri(fNeig);
        const auto piTrace = DXQY::traceLowTri(ccfNeig);
        const auto uu = DXQY::dot(uNeig, uNeig);
        const auto c2 = DXQY::c2;
        const auto uuNode = DXQY::dot(uNode, uNode);

        for (auto &delta: wallPressureBoundary_.bnd.delta(bndNo)) {
            const auto deltaHat = DXQY::reverseDirection(delta);
            const auto c_delta = DXQY::c(delta);
            const auto ccpi = DXQY::dot(DXQY::contractionLowTriVec(ccfNeig, c_delta), c_delta);
            const auto cc = DXQY::dot(c_delta, c_delta);
            const auto cu = DXQY::dot(c_delta, uNode);
            const auto cF = DXQY::dot(c_delta, forceNode);
            const auto QPi_neq = ccpi-c2*piTrace - rhoNeig*c2*(cc - c2*DXQY::nD) - rhoNeig*(cu*cu - c2*uu);
            
            f(fieldNo, delta, nodeNo) = w_[delta]*rhoNode*(1 + DXQY::c2Inv*cu + DXQY::c4Inv0_5*(cu*cu - c2*uuNode) );
            f(fieldNo, delta, nodeNo) += w_[delta]*DXQY::c4Inv0_5*QPi_neq - w_[delta]*0.5*DXQY::c2Inv*cF;

            f(fieldNo, deltaHat, nodeNo) = w_[delta]*rhoNode*(1 - DXQY::c2Inv*cu + DXQY::c4Inv0_5*(cu*cu - c2*uuNode) );
            f(fieldNo, deltaHat, nodeNo) += w_[delta]*DXQY::c4Inv0_5*QPi_neq + w_[delta]*0.5*DXQY::c2Inv*cF;
            wUnknown += 2*w_[delta];
        }     
        const lbBase_t rhoCorr = (rhoNode - DXQY::qSum(f(fieldNo, nodeNo)))/wUnknown;
        for (auto &beta: wallPressureBoundary_.bnd.beta(bndNo)) {
             f(fieldNo, beta, nodeNo) += w_[beta]*rhoCorr;
        }
        for (auto &delta: wallPressureBoundary_.bnd.delta(bndNo)) {
            const auto deltaHat = DXQY::reverseDirection(delta);
             f(fieldNo, delta, nodeNo) += w_[delta]*rhoCorr;
             f(fieldNo, deltaHat, nodeNo) += w_[delta]*rhoCorr;
        } 
    }

    //----------------------------------------------------------------------------------- pressure boundary
    for (int bndNo = 0; bndNo < pressureBoundary_.bnd.size(); ++bndNo) {
        const int nodeNo = pressureBoundary_.bnd.nodeNo(bndNo);
        const int alphaNeig = pressureBoundary_.neighbors[bndNo];
        const auto nVec = pressureBoundary_.normals[bndNo];
        const int nodeNeig = grid.neighbor(alphaNeig, nodeNo);
        const auto fNeig = f(fieldNo, nodeNeig);
        const lbBase_t rhoNode = 1.0;

        lbBase_t eqA = 0;
        lbBase_t eqB = -rhoNode;
        lbBase_t eqC = 0.5*DXQY::dot(force(fieldNo, nodeNo), nVec);
        for (auto &gamma: pressureBoundary_.bnd.gamma(bndNo)) {
            const auto gammaHat = DXQY::reverseDirection(gamma);
            const auto c = DXQY::c(gamma);
            const auto cn = DXQY::dot(c, nVec);
            eqC += (f(fieldNo, gamma, nodeNo) - f(fieldNo, gammaHat, nodeNo))*cn;
        }

        const std::valarray<lbBase_t> uNeig = DXQY::qSumC(fNeig) + 0.5*force(fieldNo, nodeNeig);
        const auto rhoNeig = DXQY::qSum(fNeig);
        const auto ccfNeig = DXQY::qSumCCLowTri(fNeig);
        const auto piTrace = DXQY::traceLowTri(ccfNeig);
        const auto uu = DXQY::dot(uNeig, uNeig);
        const auto c2 = DXQY::c2;

        for (auto &beta: pressureBoundary_.bnd.beta(bndNo)) {
            const auto betaHat = DXQY::reverseDirection(beta);
            const auto c = DXQY::c(beta);
            const auto cn = DXQY::dot(c, nVec);
            const auto ccpi = DXQY::dot(DXQY::contractionLowTriVec(ccfNeig, c), c);
            const auto cc = DXQY::dot(c, c);
            const auto cu = DXQY::dot(c, uNeig);
            const auto QPi_neq = ccpi-c2*piTrace - rhoNeig*c2*(cc - c2*DXQY::nD) - rhoNeig*(cu*cu - c2*uu);

            eqC += -f(fieldNo, betaHat, nodeNo)*cn;

            eqC += w_[beta]*rhoNode*cn;
            eqC += w_[beta]*DXQY::c4Inv0_5*QPi_neq*cn;
            eqC -= w_[beta]*0.5*DXQY::c2Inv*DXQY::dot(c, force(fieldNo, nodeNo))*cn;

            eqB += w_[beta]*DXQY::c2*rhoNode*cn*cn;

            eqA += w_[beta]*D2Q9::c4Inv0_5*rhoNode*(cn*cn - D2Q9::c2)*cn;
        }
        for (auto &delta: pressureBoundary_.bnd.delta(bndNo)) {
            const auto c = DXQY::c(delta);
            const auto cn = DXQY::dot(c, nVec);
            eqC -= w_[delta]*0.5*DXQY::c2Inv*DXQY::dot(c, force(fieldNo, nodeNo))*cn;
            eqB += w_[delta]*DXQY::c2*rhoNode*cn*cn;
        }

        const lbBase_t uNodeN = -eqC/eqB - (eqA*eqC*eqC)/(eqB*eqB*eqB);
        const std::valarray<lbBase_t> uNode = uNodeN*nVec;

        lbBase_t wUnknown = 0;

        for (auto &beta: pressureBoundary_.bnd.beta(bndNo)) {           
            const auto betaHat = DXQY::reverseDirection(beta);
            const auto c_beta = DXQY::c(beta);
            const auto cu = DXQY::dot(c_beta, uNode);
            const auto cF = DXQY::dot(c_beta, force(fieldNo, nodeNo));

            f(fieldNo, beta, nodeNo) = f(fieldNo, betaHat, nodeNo) + w_[beta]*DXQY::c2Inv*(2*rhoNode*cu - cF);         

            wUnknown += w_[beta];
        }
        const auto uuNode = DXQY::dot(uNode, uNode);

        for (auto &delta: pressureBoundary_.bnd.delta(bndNo)) {
            const auto deltaHat = DXQY::reverseDirection(delta);
            const auto c_delta = DXQY::c(delta);
            const auto ccpi = DXQY::dot(DXQY::contractionLowTriVec(ccfNeig, c_delta), c_delta);
            const auto cc = DXQY::dot(c_delta, c_delta);
            const auto cu = DXQY::dot(c_delta, uNode);
            const auto cF = DXQY::dot(c_delta, force(fieldNo, nodeNo));
            const auto QPi_neq = ccpi-c2*piTrace - rhoNeig*c2*(cc - c2*DXQY::nD) - rhoNeig*(cu*cu - c2*uu);
            
            f(fieldNo, delta, nodeNo) = w_[delta]*rhoNode*(1 + DXQY::c2Inv*cu + DXQY::c4Inv0_5*(cu*cu - c2*uuNode) );
            f(fieldNo, delta, nodeNo) += w_[delta]*DXQY::c4Inv0_5*QPi_neq - w_[delta]*0.5*DXQY::c2Inv*cF;

            f(fieldNo, deltaHat, nodeNo) = w_[delta]*rhoNode*(1 - DXQY::c2Inv*cu + DXQY::c4Inv0_5*(cu*cu - c2*uuNode) );
            f(fieldNo, deltaHat, nodeNo) += w_[delta]*DXQY::c4Inv0_5*QPi_neq + w_[delta]*0.5*DXQY::c2Inv*cF;

            wUnknown += 2*w_[delta];
        }    

        std::valarray<lbBase_t> uCorr = uNode - DXQY::qSumC(f(fieldNo, nodeNo)) - 0.5*force(fieldNo, nodeNo);
        const lbBase_t rhoCorr = (rhoNode - DXQY::qSum(f(fieldNo, nodeNo)))/wUnknown;
        for (auto &beta: pressureBoundary_.bnd.beta(bndNo)) {
             const auto c = DXQY::c(beta);
             f(fieldNo, beta, nodeNo) += w_[beta]*rhoCorr;
        }
        for (auto &delta: pressureBoundary_.bnd.delta(bndNo)) {
            const auto deltaHat = DXQY::reverseDirection(delta);
             f(fieldNo, delta, nodeNo) += w_[delta]*rhoCorr;
             f(fieldNo, deltaHat, nodeNo) += w_[delta]*rhoCorr;
        } 
    }
}

//                               WetNodeBoundary
//----------------------------------------------------------------------------------- readScalarValues
template<typename DXQY>
template<typename T>
std::vector<T> WetNodeBoundary<DXQY>::readScalarValues(const std::string attributeName, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
{
    vtklb.toAttribute(attributeName);
    std::vector<T> ret(grid.size());
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        ret[nodeNo] = vtklb.template getScalarAttribute<T>();
    }

    return ret;
}

//                               WetNodeBoundary
//----------------------------------------------------------------------------------- readVectorValues
template<typename DXQY>
template<typename T>
std::vector<std::valarray<T>> WetNodeBoundary<DXQY>::readVectorValues(const std::string attributeName, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
{
    const std::array<const std::string, 3> index{"x", "y", "z"};
    std::vector<std::valarray<T>> ret(grid.size(), std::valarray<T>(DXQY::nD));

    for (int d = 0; d < DXQY::nD; ++d) {
        vtklb.toAttribute(attributeName + "_" + index[d]);
        for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
        {
            ret[nodeNo][d] = vtklb.template getScalarAttribute<T>();
        }
    }
    return ret;
}



#endif