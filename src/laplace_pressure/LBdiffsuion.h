#ifndef LBDIFFUSION_H
#define LBDIFFUSION_H

#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBvtk.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBboundary.h"

#include <array>

//=====================================================================================
//
//                         D I F F U S I O N S O L V E R
//
//=====================================================================================
template<typename DXQY>
class DiffusionSolver
{
public:
    // Bulk diffusion
    DiffusionSolver(const lbBase_t tau, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);

    int numBoundaries() const {return numBoundaries_;}

    std::valarray<lbBase_t> setF(const lbBase_t & rho) const;

    template<typename T>
    std::valarray<lbBase_t> omegaBGK(const T & f, const lbBase_t & rho) const;

    lbBase_t diffusionCoefficient() const;

    // Boundary condition
    void applyBoundaryCondition(LbField<DXQY> &f, const Grid<DXQY> &grid) const;
    void applyBoundaryCondition(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const;

    // Set forcing
    template<typename T>
    VectorField<DXQY> getForcing(const int fieldNo, const LbField<DXQY> & f, const T & bulkNodes) const;

    // collision
    std::valarray<lbBase_t> collision(const std::valarray<lbBase_t> &fNode) const;
    std::valarray<lbBase_t> collision(const std::valarray<lbBase_t> &fNode, lbBase_t &rho) const;

    // Write forces and pressures
    template<typename T>
    void writeFieldsToFile(const LbField<DXQY> & f, const ScalarField &rho, const T & bulkNodes);

    // Help function BEGIN
    Boundary<DXQY> getWallBoundary() {return wallBoundary_.bnd;}
    Boundary<DXQY> getWallPressureBoundaryNodes() {return wallPressureBoundary_.bnd;}
    Boundary<DXQY> getPressureBoundaryNodes() {return pressureBoundary_.bnd;}
    std::vector<std::valarray<lbBase_t>> getWallNormals() {return wallBoundary_.normals;}
    std::vector<int> getWallNeighbors() {return wallBoundary_.neighbors;}
    // Help function END
private:
    // Boundary diffusion
    void setupBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
    template<typename T>
    std::vector<T> readScalarValues(const std::string attributeName, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
    template<typename T>
    std::vector<std::valarray<T>> readVectorValues(const std::string attributeName, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);

    // Definition of local structures  
    struct DiffusionBoundaryNodes
    {
        std::vector<std::valarray<lbBase_t>> normals;
        std::vector<lbBase_t> qs;
        std::vector<int> indicators;
        std::vector<int> neighbors;
        Boundary<DXQY> bnd;
    };  
    const int size_;
    const lbBase_t tau_;
    const lbBase_t tauInv_;
    const std::valarray<lbBase_t> w_;
    int numBoundaries_;

    DiffusionBoundaryNodes wallBoundary_;
    DiffusionBoundaryNodes pressureBoundary_;
    DiffusionBoundaryNodes wallPressureBoundary_;
};

//                               DiffusionSolver
//----------------------------------------------------------------------------------- DiffusionSolver
template<typename DXQY>
DiffusionSolver<DXQY>::DiffusionSolver(const lbBase_t tau, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
:size_(grid.size()), tau_(tau), tauInv_(1.0/tau), w_(DXQY::w, DXQY::nQ) 
{
    setupBoundaryNodes(vtk, nodes, grid);
    auto qvalues = readScalarValues<lbBase_t>("q", vtk, nodes, grid);
    auto normalvec = readVectorValues<lbBase_t>("normal", vtk, nodes, grid);
    auto neighvec = readVectorValues<int>("neighbor", vtk, nodes, grid);

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

//                               DiffusionSolver
//----------------------------------------------------------------------------------- setF
template<typename DXQY>
std::valarray<lbBase_t> DiffusionSolver<DXQY>::setF(const lbBase_t & rho) const
//-----------------------------------------------------------------------------------
{
    return w_*rho;
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- omegaBGK
template<typename DXQY>
template<typename T>
std::valarray<lbBase_t> DiffusionSolver<DXQY>::omegaBGK(const T & f, const lbBase_t & rho) const
//-----------------------------------------------------------------------------------
{
    return tauInv_*(w_*rho - f);
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- diffusionCoefficient
template<typename DXQY>
lbBase_t DiffusionSolver<DXQY>::diffusionCoefficient() const
//-----------------------------------------------------------------------------------
{
    return DXQY::c2 * (tau_ - 0.5);
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- applyBoundaryCondition
template<typename DXQY>
void DiffusionSolver<DXQY>::applyBoundaryCondition(LbField<DXQY> &f, const Grid<DXQY> &grid) const
//-----------------------------------------------------------------------------------
{
    for (int n=0; n < f.num_fields(); ++n) {
        applyBoundaryCondition(n, f, grid);
    }
}
//                               DiffusionSolver
//----------------------------------------------------------------------------------- applyBoundaryCondition
template<typename DXQY>
void DiffusionSolver<DXQY>::applyBoundaryCondition(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const
//-----------------------------------------------------------------------------------
{
    auto calcpi = [](const std::valarray<lbBase_t> &fNeig) {
        return DXQY::qSumCCLowTri(fNeig) - DXQY::c2 * DXQY::qSum(fNeig) * DXQY::deltaLowTri();
    };
    auto calcjBukl = [&](const std::valarray<lbBase_t> &fNeig, const int alpha, const std::valarray<lbBase_t> &pi) {
        std::valarray<lbBase_t> ret = DXQY::qSumC(fNeig);
        ret +=  DXQY::contractionLowTriVec(pi, DXQY::c(alpha)) / ( (2*tau_ - 1)*DXQY::c2 );
        return ret;
    };
    auto calcjWall = [&](const std::valarray<lbBase_t> &fNeig, const int alpha, const lbBase_t q, const std::valarray<lbBase_t> &nVec, const std::valarray<lbBase_t> &pi)  {
        const auto jNeig = DXQY::qSumC(fNeig); // \vec{j}(\vec{c}_\beta) in article
        std::valarray<lbBase_t> ret = jNeig;
        ret -= ( DXQY::dot(nVec, jNeig)/(1+q) )*nVec;
        const auto piDotC = DXQY::contractionLowTriVec(pi, DXQY::c(alpha));
        ret += ( piDotC - nVec*DXQY::dot(nVec, piDotC) ) / ( (2*tau_ - 1)*DXQY::c2 );
        return ret;
    };
    auto calcPhi = [&](const std::valarray<lbBase_t> &fNode, const std::vector<int> &gamma, const std::vector<int> &beta, const std::vector<int> &delta, const std::valarray<lbBase_t> &jVec, const std::valarray<lbBase_t> &pi) {
        lbBase_t rhs = 0.0;
        rhs += fNode[DXQY::nQNonZero_];
        for (auto & alpha : gamma) {
            rhs += fNode[alpha];
            rhs += fNode[DXQY::reverseDirection(alpha)];
        }
        for (auto alpha: beta) {
            const lbBase_t cj = DXQY::dot(DXQY::c(alpha), jVec);
            rhs += 2*fNode[DXQY::reverseDirection(alpha)] + 2*w_[alpha]*DXQY::c2Inv*cj;
        }
        lbBase_t lhs = 1.0;
        for (auto alpha: delta) {
            const lbBase_t cPic = DXQY::dot(DXQY::c(alpha),DXQY::contractionLowTriVec(pi, DXQY::c(alpha)));
            const lbBase_t piTrace = DXQY::traceLowTri(pi);
            rhs += w_[alpha]*DXQY::c4Inv*(cPic - DXQY::c2*piTrace);
            lhs -= 2*w_[alpha];
        }
        return rhs/lhs;
    };
    auto calcf = [&](const lbBase_t alpha, const lbBase_t phi, const std::valarray<lbBase_t> &jVec, const std::valarray<lbBase_t> &piMat) {
        const lbBase_t cj = DXQY::dot(DXQY::c(alpha), jVec);
        const lbBase_t ccPi = DXQY::dot(DXQY::c(alpha),DXQY::contractionLowTriVec(piMat, DXQY::c(alpha)));
        const lbBase_t iiPi = DXQY::traceLowTri(piMat);

        return w_[alpha]*(phi + DXQY::c2Inv*cj + DXQY::c4Inv0_5*(ccPi - DXQY::c2*iiPi));
    };  

    // ------------------------------------------------------------------------------ WALL BOUNDARY
    DiffusionBoundaryNodes diffBnd = wallBoundary_;
    Boundary<DXQY> bnd = diffBnd.bnd;
    for (int bndNo = 0; bndNo < bnd.size(); ++bndNo) 
    {
        const int nodeNo = bnd.nodeNo(bndNo);
        const int alphaNeig = diffBnd.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(alphaNeig, nodeNo));
        
        // Calculate the second moment
        const std::valarray<lbBase_t> piNode = calcpi(fNeig); 

        // Calculate the first moment
        const auto qNode = diffBnd.qs[bndNo];  // q in article
        const auto nNode = diffBnd.normals[bndNo]; // \vec{n} in article

        auto jNode = calcjWall(fNeig, alphaNeig, qNode, nNode, piNode);
        auto phiNode = calcPhi(f(fieldNo, nodeNo), bnd.gamma(bndNo), bnd.beta(bndNo), bnd.delta(bndNo), jNode, piNode);

        //
        for (auto &beta : bnd.beta(bndNo)) {
            auto betaHat = bnd.dirRev(beta);
            f(fieldNo, beta, nodeNo) = f(fieldNo, betaHat, nodeNo) + 2*DXQY::c2Inv*w_[beta]*DXQY::dot(DXQY::c(beta), jNode);
        }

        for (auto &delta : bnd.delta(bndNo)) {
            f(fieldNo, delta, nodeNo) = calcf(delta, phiNode, jNode, piNode);
            auto deltaHat = bnd.dirRev(delta);
            f(fieldNo, deltaHat, nodeNo) = calcf(deltaHat, phiNode, jNode, piNode);
        }
    }

    // ------------------------------------------------------------------------------ WALL PRESSURE BOUNDARY 
    diffBnd = wallPressureBoundary_;
    bnd = diffBnd.bnd;
    for (int bndNo = 0; bndNo < bnd.size(); ++bndNo) 
    {
        const int nodeNo = bnd.nodeNo(bndNo);
        const int alphaNeig = diffBnd.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(alphaNeig, nodeNo));        
        // Calculate the second moment
        const std::valarray<lbBase_t> piNode = calcpi(fNeig); 
        // Calculate the first moment
        const auto qNode = diffBnd.qs[bndNo];  // q in article
        const auto nNode = diffBnd.normals[bndNo]; // \vec{n} in article
        auto jNode = calcjWall(fNeig, alphaNeig, qNode, nNode, piNode);
        // Calculate the zeroth moment
        auto phiNode = 1.0 * (diffBnd.indicators[bndNo] == (fieldNo+1));

        auto alphaZero = DXQY::nQNonZero_;
        f(fieldNo, alphaZero, nodeNo) = calcf(alphaZero, phiNode, jNode, piNode);
        auto setfs = [&](const std::vector<int> &alphaList) {
            for (auto & alpha: alphaList) {
                f(fieldNo, alpha, nodeNo) = calcf(alpha, phiNode, jNode, piNode);
                auto alphaHat = bnd.dirRev(alpha);
                f(fieldNo, alphaHat, nodeNo) = calcf(alphaHat, phiNode, jNode, piNode);
            }
        };
        setfs(bnd.gamma(bndNo));
        setfs(bnd.beta(bndNo));
        setfs(bnd.delta(bndNo));
    }

    // ------------------------------------------------------------------------------ PRESSURE BOUNDARY 
    diffBnd = pressureBoundary_;
    bnd = diffBnd.bnd;
    for (int bndNo = 0; bndNo < bnd.size(); ++bndNo) 
    {
        const int nodeNo = bnd.nodeNo(bndNo);
        const int alphaNeig = diffBnd.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(alphaNeig, nodeNo));        
        // Calculate the second moment
        const std::valarray<lbBase_t> piNode = calcpi(fNeig); 
        // Calculate the first moment
        const auto qNode = diffBnd.qs[bndNo];  // q in article
        const auto nNode = diffBnd.normals[bndNo]; // \vec{n} in article
        auto jNode = calcjBukl(fNeig, alphaNeig, piNode);
        // Calculate the zeroth moment
        auto phiNode = 1.0 * (diffBnd.indicators[bndNo] == (fieldNo+1));

        auto alphaZero = DXQY::nQNonZero_;
        f(fieldNo, alphaZero, nodeNo) = calcf(alphaZero, phiNode, jNode, piNode);
        auto setfs = [&](const std::vector<int> &alphaList) {
            for (auto & alpha: alphaList) {
                f(fieldNo, alpha, nodeNo) = calcf(alpha, phiNode, jNode, piNode);
                auto alphaHat = bnd.dirRev(alpha);
                f(fieldNo, alphaHat, nodeNo) = calcf(alphaHat, phiNode, jNode, piNode);
            }
        };
        setfs(bnd.gamma(bndNo));
        setfs(bnd.beta(bndNo));
        setfs(bnd.delta(bndNo));
    }
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- setupBoundaryNodes
template<typename DXQY>
void DiffusionSolver<DXQY>::setupBoundaryNodes(LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
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
        // Setup the different boundary nodes.
        wallBoundary_.bnd = Boundary<DXQY>(wallBoundaryNodes, nodes, grid);
        pressureBoundary_.bnd =  Boundary<DXQY>(pressureBoundaryNodes, nodes, grid);
        wallPressureBoundary_.bnd = Boundary<DXQY>(wallPressureBoundaryNodes, nodes, grid);     
    }
    numBoundaries_ = numBnd;
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- readScalarValues
template<typename DXQY>
template<typename T>
std::vector<T> DiffusionSolver<DXQY>::readScalarValues(const std::string attributeName, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
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

//                               DiffusionSolver
//----------------------------------------------------------------------------------- readVectorValues
template<typename DXQY>
template<typename T>
std::vector<std::valarray<T>> DiffusionSolver<DXQY>::readVectorValues(const std::string attributeName, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
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

//                               DiffusionSolver
//----------------------------------------------------------------------------------- getForcing
template<typename DXQY>
template<typename T>
VectorField<DXQY> DiffusionSolver<DXQY>::getForcing(const int fieldNo, const LbField<DXQY> & f, const T & bulkNodes) const
//-----------------------------------------------------------------------------------
{
    VectorField<DXQY> ret(1, f.getNumNodes());

    for (auto & nodeNo: bulkNodes) {
        ret.set(0, nodeNo) = (DXQY::c2Inv * tauInv_) * DXQY::qSumC(f(fieldNo, nodeNo));
    }    
    return ret;
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- collision
template<typename DXQY>    
std::valarray<lbBase_t> DiffusionSolver<DXQY>::collision(const std::valarray<lbBase_t> &fNode) const
//-----------------------------------------------------------------------------------
{
    const lbBase_t rhoNode = calcRho<DXQY>(fNode);
    return fNode + omegaBGK(fNode, rhoNode);
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- collision
template<typename DXQY>    
std::valarray<lbBase_t> DiffusionSolver<DXQY>::collision(const std::valarray<lbBase_t> &fNode, lbBase_t &rho) const
//-----------------------------------------------------------------------------------
{
    const lbBase_t rhoNode = calcRho<DXQY>(fNode);
    rho = rhoNode;
    return fNode + omegaBGK(fNode, rhoNode);
}

//                               DiffusionSolver
//-----------------------------------------------------------------------------------
template<typename DXQY>  
template<typename T>
void DiffusionSolver<DXQY>::writeFieldsToFile(const LbField<DXQY> & f, const ScalarField &rho, const T & bulkNodes)
{

}

#endif