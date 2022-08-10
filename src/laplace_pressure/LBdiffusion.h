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

    int maxBoundaryIndicator() const {return maxPressureInidcator_;}

    std::valarray<lbBase_t> setF(const lbBase_t & rho) const;

    template<typename T>
    std::valarray<lbBase_t> omegaBGK(const T & f, const lbBase_t & rho) const;

    lbBase_t diffusionCoefficient() const;

    // Geometry
    std::vector<int> findBulkNodes(LBvtk<DXQY> &vtk, const Nodes<DXQY> &nodes);

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
    std::valarray<lbBase_t> getNormals(const int nodeNo) {return normalVector(0, nodeNo);}
    lbBase_t getSignedDistance(const int nodeNo) {return signedDistance(0, nodeNo);}
//    void fillBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
    ScalarField fillBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);

    // Boundary diffusion
    ScalarField setupBoundaryNodesTest(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);


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
    int maxPressureInidcator_;

    DiffusionBoundaryNodes wallBoundary_;
    DiffusionBoundaryNodes pressureBoundary_;
    DiffusionBoundaryNodes wallPressureBoundary_;

    // Test varaibles
    VectorField<DXQY> normalVector;
    ScalarField signedDistance;
};

//                               DiffusionSolver
//----------------------------------------------------------------------------------- DiffusionSolver
//                               DiffusionSolver
//----------------------------------------------------------------------------------- fillBoundaryNodes
template<typename DXQY>
//void DiffusionSolver<DXQY>::fillBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
ScalarField DiffusionSolver<DXQY>::fillBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
{

    // All nodes
    auto sd = readScalarValues<lbBase_t>("signed_distance", vtk, nodes, grid);
    for (int nodeNo=0; nodeNo<grid.size(); ++nodeNo)
    {
        signedDistance(0, nodeNo) = sd[nodeNo];
        if ( nodes.isFluid(nodeNo) )
        {
            std::vector<lbBase_t> sdNeigs(DXQY::nQ);
            int  cnt = 0;
            for (const auto &neigNo : grid.neighbor(nodeNo)) {
                sdNeigs[cnt] = sd[neigNo];
                cnt += 1;
            }
            std::valarray<lbBase_t> grad = DXQY::grad(sdNeigs);
            auto norm = std::sqrt(DXQY::dot(grad, grad));
            if (norm < lbBaseEps) 
                norm = 1.0;
            grad /= norm;
            normalVector.set(0, nodeNo) = grad;
        }
    }

    // Pressure and wall-pressure boundary normal
    // Read the boundary normal from the vtklb file
    std::vector<std::string> cartExt{"_x", "_y", "_z"};
    for (int d = 0; d < DXQY::nD; ++d) 
    { 
        vtk.toAttribute("boundary_normal" + cartExt[d]);
        for (int n=0; n < vtk.numSubsetEntries(); ++n) {
            const auto att = vtk.template getSubsetAttribure<lbBase_t>();
            normalVector(0, d, att.nodeNo) = att.val; 
            //std::cout << "node[" << d << "] = " << att.nodeNo << " " << att.val << std::endl;
        }
    }

    /* std::valarray<lbBase_t> givenNorm{1.0, 0.0};
    for (int bndNo = 0; bndNo<pressureBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = pressureBoundary_.bnd.nodeNo(bndNo);   
        normalVector.set(0, nodeNo) = givenNorm;
    }

    for (int bndNo = 0; bndNo<wallPressureBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = wallPressureBoundary_.bnd.nodeNo(bndNo);   
        normalVector.set(0, nodeNo) = givenNorm;
    } */


    ScalarField vtkPrint(3, grid.size());

    // Set the bulk neighbors  and qvalues
    std::vector<lbBase_t> qValues(grid.size());
    std::vector<int> bulkNeigh(grid.size());


    // wall nodes
    for (int bndNo = 0; bndNo<wallBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = wallBoundary_.bnd.nodeNo(bndNo);  
        // Find the lattice direction that best matches the normal direction
        auto ncVec = DXQY::cDotAll(normalVector(0, nodeNo));
        lbBase_t maxVal = 0;
        int bulkDir = -1;
        lbBase_t sVal = 0;

        for (int q=0; q < DXQY::nQ-1; ++q)
        {            
            const int neigNo = grid.neighbor(q, nodeNo);
            const int neigRevNo = grid.neighbor(DXQY::reverseDirection(q), nodeNo);
            if ( nodes.isMyRank(neigNo) && nodes.isSolid(neigRevNo) && (!nodes.isDefault(neigRevNo)) && nodes.isFluid(neigNo) && (ncVec[q]/DXQY::cNorm[q] > maxVal) ) {
                bulkDir = q;
                maxVal = ncVec[q]/DXQY::cNorm[q];
                sVal = sd[nodeNo]/(sd[nodeNo] - sd[neigRevNo]);
            }
        }
        if (bulkDir == -1) {
            std::cout << "Could not find a bulk node for the diffusion boundary condition" << std::endl;
	        std::cout << "@ node no: " << nodeNo << ", pos ";
	        for(const auto &i: grid.pos(nodeNo))
	            std::cout << i << " ";
	        std::cout << std::endl;
            vtkPrint(2, nodeNo) = 1;
	        // exit(1);
        } 
        bulkNeigh[nodeNo] = bulkDir;
        qValues[nodeNo] = sVal;

        vtkPrint(0, nodeNo) = qValues[nodeNo];
        vtkPrint(1, nodeNo) = bulkNeigh[nodeNo];
    }


    // pressure nodes
    // lbBase_t s_from_file = 0;
    for (int bndNo = 0; bndNo<pressureBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = pressureBoundary_.bnd.nodeNo(bndNo);  

        // Find the lattice direction that best matches the normal direction
        auto ncVec = DXQY::cDotAll(normalVector(0, nodeNo));
        lbBase_t maxVal = 0;
        int bulkDir = -1;
        //lbBase_t sVal = 0;
        for (int q=0; q < DXQY::nQ-1; ++q)
        {
            const int neigNo = grid.neighbor(q, nodeNo);
            // const int neigRevNo = grid.neighbor(DXQY::reverseDirection(q), nodeNo);
            if ( nodes.isMyRank(neigNo) && nodes.isFluid(neigNo) && (ncVec[q]/DXQY::cNorm[q] > maxVal) ) {
                bulkDir = q;
                maxVal = ncVec[q]/DXQY::cNorm[q];
                // sVal = s_from_file; // sd[nodeNo]/(sd[nodeNo] - sd[neigRevNo]);
            }
        }
        if (bulkDir == -1) {
            std::cout << "Could not find a bulk node for the diffusion boundary condition" << std::endl;
            vtkPrint(2, nodeNo) = 2;
            // exit(1);
        }
        bulkNeigh[nodeNo] = bulkDir;
        // qValues[nodeNo] = sVal;

        //vtkPrint(0, nodeNo) = qValues[nodeNo];
        vtkPrint(1, nodeNo) = bulkNeigh[nodeNo];
    }


    // wall pressure nodes
    for (int bndNo = 0; bndNo<wallPressureBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = wallPressureBoundary_.bnd.nodeNo(bndNo);  

        // Find the lattice direction that best matches the normal direction
        auto ncVec = DXQY::cDotAll(normalVector(0, nodeNo));
        lbBase_t maxVal = 0;
        int bulkDir = -1;
        // lbBase_t sVal = 0;
        for (int q=0; q < DXQY::nQ-1; ++q)
        {
            const int neigNo = grid.neighbor(q, nodeNo);
            //const int neigRevNo = grid.neighbor(DXQY::reverseDirection(q), nodeNo);
            if ( nodes.isMyRank(neigNo) && nodes.isFluid(neigNo) && (ncVec[q]/DXQY::cNorm[q] > maxVal) ) {
                bulkDir = q;
                maxVal = ncVec[q]/DXQY::cNorm[q];
                // sVal = s_from_file; //  sd[nodeNo]/(sd[nodeNo] - sd[neigRevNo]);
            }
        }
        if (bulkDir == -1) {
            std::cout << "Could not find a bulk node for the diffusion boundary condition" << std::endl;
            vtkPrint(2, nodeNo) = 3;
            // exit(1);
        }
        bulkNeigh[nodeNo] = bulkDir;
        // qValues[nodeNo] = sVal;

        // vtkPrint(0, nodeNo) = qValues[nodeNo];
        vtkPrint(1, nodeNo) = bulkNeigh[nodeNo];
    }


    vtk.toAttribute("boundary_distance");
    for (int n=0; n < vtk.numSubsetEntries(); ++n) {
        const auto att = vtk.template getSubsetAttribure<lbBase_t>();
        qValues[att.nodeNo] = att.val/DXQY::cDotRef(bulkNeigh[att.nodeNo], normalVector(0, att.nodeNo)); 
        vtkPrint(0, att.nodeNo) = qValues[att.nodeNo];
    }

//    auto qvalues = readScalarValues<lbBase_t>("q", vtk, nodes, grid);
    // auto normalvec = readVectorValues<lbBase_t>("normal", vtk, nodes, grid);
//    auto neighvec = readVectorValues<int>("neighbor", vtk, nodes, grid);

    auto fillBoundaryNodes = [this, &qValues, &bulkNeigh](DiffusionBoundaryNodes &diffBnd) 
    {
        for (int n = 0; n < diffBnd.bnd.size(); ++n) {
            const int nodeNo = diffBnd.bnd.nodeNo(n);
            diffBnd.qs.push_back(qValues[nodeNo]);
            diffBnd.normals.push_back(normalVector(0, nodeNo));
//            std::vector<int> nvec(DXQY::nD);
//            for (int d=0; d < DXQY::nD; ++d)  nvec[d] = neighvec[nodeNo][d];
            diffBnd.neighbors.push_back(bulkNeigh[nodeNo]);
        }        
    };

/*    auto fillBoundaryNodes = [this, &qvalues, &neighvec](DiffusionBoundaryNodes &diffBnd) 
    {
        for (int n = 0; n < diffBnd.bnd.size(); ++n) {
            const int nodeNo = diffBnd.bnd.nodeNo(n);
            diffBnd.qs.push_back(qvalues[nodeNo]);
            diffBnd.normals.push_back(normalVector(0, nodeNo));
            std::vector<int> nvec(DXQY::nD);
            for (int d=0; d < DXQY::nD; ++d)  nvec[d] = neighvec[nodeNo][d];
            diffBnd.neighbors.push_back(DXQY::c2q(nvec));
        }        
    }; */

    fillBoundaryNodes(wallBoundary_);
    fillBoundaryNodes(pressureBoundary_);
    fillBoundaryNodes(wallPressureBoundary_);

    return vtkPrint;

}


template<typename DXQY>
DiffusionSolver<DXQY>::DiffusionSolver(const lbBase_t tau, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
:size_(grid.size()), tau_(tau), tauInv_(1.0/tau), w_(DXQY::w, DXQY::nQ), normalVector(1, grid.size()), signedDistance(1, grid.size()) 
{
    
    setupBoundaryNodes(vtk, nodes, grid);
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
//----------------------------------------------------------------------------------- findBulkNodes
template<typename DXQY>
std::vector<int> DiffusionSolver<DXQY>::findBulkNodes(LBvtk<DXQY> &vtk, const Nodes<DXQY> &nodes)
{
    std::vector<int> bulkNodes;

    vtk.toAttribute("pressure_boundary");
    for (int n = vtk.beginNodeNo(); n < vtk.endNodeNo(); ++n) 
    {
        const auto pInd = vtk.template getScalarAttribute<int>();
        if ( nodes.isFluid(n) && nodes.isMyRank(n) && (pInd >= 0) )
            bulkNodes.push_back(n);
    }
    return bulkNodes;    
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
    const DiffusionBoundaryNodes &diffBnd = wallBoundary_;
    const Boundary<DXQY> &bnd = diffBnd.bnd;
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
    const DiffusionBoundaryNodes &diffBndWP = wallPressureBoundary_;
    const Boundary<DXQY> &bndWP = diffBndWP.bnd;
    for (int bndNo = 0; bndNo < bndWP.size(); ++bndNo)
    {
        const int nodeNo = bndWP.nodeNo(bndNo);
        const int alphaNeig = diffBndWP.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(alphaNeig, nodeNo));
        // Calculate the second moment
        const std::valarray<lbBase_t> piNode = calcpi(fNeig);
        // Calculate the first moment
        const auto qNode = diffBndWP.qs[bndNo]; // q in article
        // const auto nNode = diffBndWP.normals[bndNo]; // \vec{n} in article
        auto jNode = calcjBukl(fNeig, alphaNeig, piNode); // calcjWall(fNeig, alphaNeig, qNode, nNode, piNode);
        // Calculate the zeroth moment
        const lbBase_t phiWall = 1.0 * (diffBndWP.indicators[bndNo] == (fieldNo + 1));
        const lbBase_t phiBulk = DXQY::qSum(fNeig);
        const lbBase_t ccPi = DXQY::dot(DXQY::c(alphaNeig), DXQY::contractionLowTriVec(piNode, DXQY::c(alphaNeig)));

        const lbBase_t phiNode = (phiWall + qNode * phiBulk) / (1 + qNode) - 0.5 * qNode * DXQY::c4Inv * ccPi / ((2 * tau_ - 1.0) * tau_);

        auto alphaZero = DXQY::nQNonZero_;
        f(fieldNo, alphaZero, nodeNo) = calcf(alphaZero, phiNode, jNode, piNode);
        auto setfsRev = [&](const std::vector<int> &alphaList)
        {
            for (auto &alpha : alphaList)
            {
                auto alphaHat = bndWP.dirRev(alpha);
                f(fieldNo, alphaHat, nodeNo) = calcf(alphaHat, phiNode, jNode, piNode);
            }
        };
        auto setfs = [&](const std::vector<int> &alphaList)
        {
            for (auto &alpha : alphaList)
            {
                f(fieldNo, alpha, nodeNo) = calcf(alpha, phiNode, jNode, piNode);
            }
        };
        // setfs(bndWP.gamma(bndNo));
        setfs(bndWP.beta(bndNo));
        setfs(bndWP.delta(bndNo));
        setfsRev(bndWP.delta(bndNo));
    }
    
    // ------------------------------------------------------------------------------ PRESSURE BOUNDARY
    const DiffusionBoundaryNodes &diffBndP = pressureBoundary_;
    const Boundary<DXQY>  &bndP = diffBndP.bnd;
    for (int bndNo = 0; bndNo < bndP.size(); ++bndNo)
    {
        const int nodeNo = bndP.nodeNo(bndNo);
        const int alphaNeig = diffBndP.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(alphaNeig, nodeNo));
        // Calculate the second moment
        const std::valarray<lbBase_t> piNode = calcpi(fNeig);
        // Calculate the first moment
        const auto qNode = diffBndP.qs[bndNo];  // q in article
        //const auto nNode = diffBndP.normals[bndNo]; // \vec{n} in article
        auto jNode = calcjBukl(fNeig, alphaNeig, piNode);//calcjWall(fNeig, alphaNeig, qNode, nNode, piNode);
        // Calculate the zeroth moment
        const lbBase_t phiWall = 1.0 * (diffBndP.indicators[bndNo] == (fieldNo+1));
        const lbBase_t phiBulk = DXQY::qSum(fNeig);
        const lbBase_t ccPi = DXQY::dot(DXQY::c(alphaNeig),DXQY::contractionLowTriVec(piNode, DXQY::c(alphaNeig)));

        const lbBase_t phiNode = (phiWall + qNode*phiBulk)/(1 + qNode) - 0.5*qNode*DXQY::c4Inv*ccPi/((2*tau_-1.0)*tau_);

        auto alphaZero = DXQY::nQNonZero_;
        f(fieldNo, alphaZero, nodeNo) = calcf(alphaZero, phiNode, jNode, piNode);

    /*        auto setfs = [&](const std::vector<int> &alphaList) {
                for (auto & alpha: alphaList) {
                    f(fieldNo, alpha, nodeNo) = calcf(alpha, phiNode, jNode, piNode);
                    auto alphaHat = bndP.dirRev(alpha);
                    f(fieldNo, alphaHat, nodeNo) = calcf(alphaHat, phiNode, jNode, piNode);
                }
            };
            setfs(bndP.gamma(bndNo));
            setfs(bndP.beta(bndNo));
            setfs(bndP.delta(bndNo));
    */
        auto setfsRev = [&](const std::vector<int> &alphaList) {
                for (auto & alpha: alphaList) {
                    auto alphaHat = bndP.dirRev(alpha);
                    f(fieldNo, alphaHat, nodeNo) = calcf(alphaHat, phiNode, jNode, piNode);
                }
            };
            auto setfs = [&](const std::vector<int> &alphaList) {
                for (auto & alpha: alphaList) {
                    f(fieldNo, alpha, nodeNo) = calcf(alpha, phiNode, jNode, piNode);
                }
            };
            setfs(bndP.beta(bndNo));
            setfs(bndP.delta(bndNo));
            setfsRev(bndP.delta(bndNo));

    } 
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- setupBoundaryNodesTest
template<typename DXQY>
ScalarField DiffusionSolver<DXQY>::setupBoundaryNodesTest(LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
{
    std::vector<int> wallBoundaryNodes;
    std::vector<int> pressureBoundaryNodes;
    std::vector<int> wallPressureBoundaryNodes;

    ScalarField ret(1, grid.size());

    // Read pressure_boundary
    int maxPressureInidcator = 0;
    vtklb.toAttribute("pressure_boundary");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        const auto pInd = vtklb.template getScalarAttribute<int>();
        maxPressureInidcator = std::max(maxPressureInidcator, pInd);

        ret(0, nodeNo) = 0;

        if ( (nodes.isFluidBoundary(nodeNo)  || (pInd > 0)) && nodes.isFluid(nodeNo) && nodes.isMyRank(nodeNo) && (pInd >= 0)) {
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
                    //exit(1);
                    ret(0, nodeNo) = 1;
                }
            }
        }  
        /* // Setup the different boundary nodes.
        wallBoundary_.bnd = Boundary<DXQY>(wallBoundaryNodes, nodes, grid);
        pressureBoundary_.bnd =  Boundary<DXQY>(pressureBoundaryNodes, nodes, grid);
        wallPressureBoundary_.bnd = Boundary<DXQY>(wallPressureBoundaryNodes, nodes, grid); */  
    }

    return ret;
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
    int maxPressureInidcator = 0;
    vtklb.toAttribute("pressure_boundary");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        const auto pInd = vtklb.template getScalarAttribute<int>();
        maxPressureInidcator = std::max(maxPressureInidcator, pInd);
        if ( (nodes.isFluidBoundary(nodeNo)  || (pInd > 0)) && nodes.isFluid(nodeNo) && nodes.isMyRank(nodeNo) && (pInd >= 0)) {
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
        /* // Setup the different boundary nodes.
        wallBoundary_.bnd = Boundary<DXQY>(wallBoundaryNodes, nodes, grid);
        pressureBoundary_.bnd =  Boundary<DXQY>(pressureBoundaryNodes, nodes, grid);
        wallPressureBoundary_.bnd = Boundary<DXQY>(wallPressureBoundaryNodes, nodes, grid); */  
    }
    // Setup the different boundary nodes.
    wallBoundary_.bnd = Boundary<DXQY>(wallBoundaryNodes, nodes, grid);
    pressureBoundary_.bnd =  Boundary<DXQY>(pressureBoundaryNodes, nodes, grid);
    wallPressureBoundary_.bnd = Boundary<DXQY>(wallPressureBoundaryNodes, nodes, grid);        
    maxPressureInidcator_ = maxPressureInidcator;
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
