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
    void applyBoundaryCondition(const int fieldNo, LbField<DXQY> &f, const VectorField<DXQY> & froce, const Grid<DXQY> &grid) const;
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
        const int alphaNeig = wallBoundary_.diffBnd.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(alphaNeig, nodeNo));

        // Calculate first moment
        const auto qNode = wallBoundary_.qs[bndNo];
        const std::valarray<lbBase_t> uNeig = D2Q9::qSumC(fNeig) + 0.5*force(fieldNo, nodeNo);
        const auto uNode = qNode*uNeig/(1 +  qNode);

        lbBase_t lhs = 1.0;
        lbBase_t rhs = 0.0;

        for (auto &gamma: wallBoundary_.bnd.gamma(bndNo)) {
            lhs -= f(fieldNo, gamma, nodeNo);
            lhs -= f(fieldNo, D2Q9::reverseDirection(delta), nodeNo);
        }

        for (auto &beta: wallBoundary_.bnd.beta(bndNo)) {
            const auto betaHat = D2Q9::reverseDirection(beta);
        }
    }
    //----------------------------------------------------------------------------------- wall pressure boundary
    //----------------------------------------------------------------------------------- pressure boundary

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