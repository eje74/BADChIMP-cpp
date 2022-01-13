#ifndef LBDIFFUSION_H
#define LBDIFFUSION_H

#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBvtk.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBboundary.h"

#include <array>

template<typename DXQY>
class DiffusionSolver
{
public:
    // Bulk diffusion
    DiffusionSolver(const lbBase_t tau, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
    std::valarray<lbBase_t> setF(const lbBase_t & rho) const;
    template<typename T>
    std::valarray<lbBase_t> omegaBGK(const T & f, const lbBase_t & rho) const;
    lbBase_t diffusionCoefficient() const;

    // Boundary condition
    void applyBoundaryCondition(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const;

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
    const lbBase_t tau_;
    const lbBase_t tauInv_;
    const std::valarray<lbBase_t> w_;

    DiffusionBoundaryNodes wallBoundary_;
    DiffusionBoundaryNodes pressureBoundary_;
    DiffusionBoundaryNodes wallPressureBoundary_;
};

/* *******************************************************************
 *  Bulk diffusion
 * ****************************************************************** */
template<typename DXQY>
DiffusionSolver<DXQY>::DiffusionSolver(const lbBase_t tau, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
:tau_(tau), tauInv_(1.0/tau), w_(DXQY::w, DXQY::nQ) 
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

template<typename DXQY>
std::valarray<lbBase_t> DiffusionSolver<DXQY>::setF(const lbBase_t & rho) const
{
    return w_*rho;
}

template<typename DXQY>
template<typename T>
std::valarray<lbBase_t> DiffusionSolver<DXQY>::omegaBGK(const T & f, const lbBase_t & rho) const
{
    return tauInv_*(w_*rho - f);
}

template<typename DXQY>
lbBase_t DiffusionSolver<DXQY>::diffusionCoefficient() const
{
    return DXQY::c2 * (tau_ - 0.5);
}

/* *******************************************************************
 *  Boundary diffusion
 * ****************************************************************** */
template<typename DXQY>
void DiffusionSolver<DXQY>::applyBoundaryCondition(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const
{
    const DiffusionBoundaryNodes diffBnd = wallBoundary_;
    const Boundary<DXQY> bnd = diffBnd.bnd;

    for (int bndNo = 0; bndNo < bnd.size(); ++bndNo) 
    {
        const int nodeNo = bnd.nodeNo(bndNo);
        const int alphaNeig = diffBnd.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(nodeNo, alphaNeig));
        
        // Calculate the second moment
        const auto piNode = DXQY::qSumCCLowTri(fNeig); // \Pi in article

        // Calculate the first moment
        const auto qNode = diffBnd.qs[bndNo];  // q in article
        const auto nNode = diffBnd.normals[bndNo]; // \vec{n} in article
        const auto jNeig = DXQY::qSumC(fNeig); // \vec{j}(\vec{c}_\beta) in article

        auto jNode = jNeig;
        jNode -= ( DXQY::dot(nNode, jNeig)/(1+qNode) )*nNode;
        const auto piDotC = DXQY::contractionLowTriVec(piNode, DXQY::c(alphaNeig));
        jNode += ( piDotC - nNode*DXQY::dot(nNode, piDotC) ) / ( (2*tau_ - 1)*DXQY::c2 );

        if (bnd.nDelta(bndNo) > 0)  {
            std::cout << "Number of delta should be 0, but is " << bnd.nDelta(bndNo) << std::endl;
            exit(1);
        }    

        // Set the unknown distributions
        for (auto &beta : bnd.beta(bndNo)) {
            auto betaHat = bnd.dirRev(beta);
            f(fieldNo, nodeNo, beta) = f(fieldNo, nodeNo, betaHat) 
        }
    }
}

template<typename DXQY>
void DiffusionSolver<DXQY>::setupBoundaryNodes(LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
{
    std::vector<int> wallBoundaryNodes;
    std::vector<int> pressureBoundaryNodes;
    std::vector<int> wallPressureBoundaryNodes;

    // Read pressure_boundary
    vtklb.toAttribute("pressure_boundary");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        const auto pInd = vtklb.template getScalarAttribute<int>();
        if ( nodes.isFluidBoundary(nodeNo) ) {
            
            bool hasSolidNeighbors = false;
            for (auto neighNo: grid.neighbor(nodeNo)) {
                if ( nodes.isBulkSolid(neighNo) || nodes.isSolidBoundary(neighNo) ) {
                    hasSolidNeighbors = true;
                }
            }

            if (hasSolidNeighbors) {
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
}

template<typename DXQY>
template<typename T>
std::vector<T> DiffusionSolver<DXQY>::readScalarValues(const std::string attributeName, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
{
    vtklb.toAttribute(attributeName);
    std::vector<T> ret(grid.size());
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        ret[nodeNo] = vtklb.template getScalarAttribute<T>();
    }

    return ret;
}

template<typename DXQY>
template<typename T>
std::vector<std::valarray<T>> DiffusionSolver<DXQY>::readVectorValues(const std::string attributeName, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
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