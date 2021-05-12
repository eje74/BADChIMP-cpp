#ifndef LBREGULARBOUNDARYBASICII_H
#define LBREGULARBOUNDARYBASICII_H

#include <vector>
#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBgrid.h"
#include "./LBboundarybasic.h"
//  Linear algebra package
#include "../Eigen/Dense"
#include "../Eigen/SVD"


template <typename DXQY>
class RegularBoundaryBasicII
{
public:
    RegularBoundaryBasicII(
        const std::vector<int>  & boundaryNodes,
        const ScalarField       & qDist,
        const VectorField<DXQY> & normals,
        const VectorField<DXQY> & tangents,
        const ScalarField       & rho,
        const VectorField<DXQY> & force,
        const lbBase_t            tau,
        const Nodes<DXQY>       & nodes,
        const Grid<DXQY>        & grid
    );
    
    void apply(
        const int                 fieldNo,
        LbField<DXQY>           & f,
        const ScalarField       & qDist,
        const VectorField<DXQY> & normals,
        const VectorField<DXQY> & tangents,
        const VectorField<DXQY> & force,
        const lbBase_t            tau,
        const Nodes<DXQY>       & nodes,         
        const Grid<DXQY>        & grid
    );

    VectorField<DXQY> bndVel_;
    ScalarField bndRho_;

private:
    BoundaryBasic<DXQY> bnd_;
    std::vector<BDCSVD<MatrixXd>> svd_;
    lbBase_t boundaryMass_;
};



/********************************
 *      class constructor
 ********************************/
template <typename DXQY>
RegularBoundaryBasicII<DXQY>::RegularBoundaryBasicII(
        const std::vector<int>  & boundaryNodes,
        const ScalarField       & qDist,
        const VectorField<DXQY> & normals,
        const VectorField<DXQY> & tangents,
        const ScalarField       & rho,
        const VectorField<DXQY> & force,
        const lbBase_t            tau,
        const Nodes<DXQY>       & nodes,
        const Grid<DXQY>        & grid
):
    bnd_(boundaryNodes, nodes, grid),
    svd_(boundaryNodes.size()),
    bndVel_(1, boundaryNodes.size()),
    bndRho_(1, boundaryNodes.size())
{
    // Initiation for mass conservation
    lbBase_t boudaryMass_ = 0;
    
    // Setup the SVD matrices
    for (int bndNo = 0; bndNo < bnd_.size(); ++bndNo) {
        auto nodeNo = bnd_.nodeNo(bndNo);

        //  Matrix declaration
        int nRows = bnd_.nKnownDirs(bndNo);
        int nColumns = 3;
        MatrixXd m(nRows, nColumns);

        // Fill matrix
        auto ct = DXQY::cDotAll(tangents(0, nodeNo));
        auto cn = DXQY::cDotAll(normals(0, nodeNo));
        auto q = qDist(0, nodeNo);
        
        int rowNo = 0;
        for (auto alpha : bnd_.knownDirs(bndNo)) {
            auto w = DXQY::w[alpha];
            m(rowNo, 0) = 0;
            m(rowNo, 1) = 0;
            m(rowNo, 2) = 0;
            // rho
            m(rowNo, 0) = w;
            // velocity
            m(rowNo, 1) += -(w * q * DXQY::c4Inv / tau) * ct[alpha];
            m(rowNo, 2) += -(w * q * DXQY::c4Inv0_5 / tau) * cn[alpha];
            // strain rate
            m(rowNo, 1) += (w * DXQY::c4Inv   ) * (ct[alpha]*cn[alpha]);
            m(rowNo, 2) += (w * DXQY::c4Inv0_5) * (cn[alpha]*cn[alpha] - DXQY::c2);
            rowNo++;
        }     
 
        // Setup Singular Value Decomposition
        svd_[bndNo] = svd_[bndNo].compute(m, ComputeThinU | ComputeThinV);

        // Check that the boundary conditions has enough data to solve the system.
        if (svd_[bndNo].rank() < nColumns) {
            std::cout << "Error in the boundary conditions: rank (" << svd_[bndNo].rank() << ") is less than number of unknowns (" << nColumns << ")." << std::endl;
            std::cout << m << std::endl;
            for (int q = 0; q < DXQY::nQ; ++q)
                std::cout << cn[q] << " ";
            std::cout << std::endl;
            for (int q = 0; q < DXQY::nQ; ++q)
                std::cout << ct[q] << " ";
            std::cout << std::endl;
            std::cout << q << std::endl;
            std::cout << tau << std::endl;
            exit(1);
        }
        
        // Mass conservation
        boundaryMass_ += rho(0, nodeNo) - 1;
    }
}



/****************************************
 * Apply boundary condition
 * **************************************/
template <typename DXQY>
void RegularBoundaryBasicII<DXQY>::apply(
    const int                 fieldNo,
    LbField<DXQY>           & f,
    const ScalarField       & qDist,
    const VectorField<DXQY> & normals,
    const VectorField<DXQY> & tangents,
    const VectorField<DXQY> & force,
    const lbBase_t            tau,
    const Nodes<DXQY>       & nodes,
    const Grid<DXQY>        & grid
)
{
    // Mass conservation
    lbBase_t changeInBulkMass = 0;
    lbBase_t boundaryMass = 0;
    lbBase_t boundarySize = 0;

    // No slip boundary condition
    for (int bndNo = 0; bndNo < bnd_.size(); bndNo++) {

        auto nodeNo = bnd_.nodeNo(bndNo);

        // Change in bulk mass
        for (int alpha = 0; alpha < DXQY::nQ; ++alpha) {
            auto alphaRev = bnd_.dirRev(alpha);
            auto nodeNeigNo = grid.neighbor(alphaRev, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo)) {
                changeInBulkMass += f(fieldNo, alphaRev, nodeNeigNo) - f(fieldNo, alpha, nodeNo);
            }
        }

        // Setup the right hand side vector
        auto nRows = bnd_.nKnownDirs(bndNo);
        VectorXd rhs(nRows);

        auto cF = DXQY::cDotAll(force(0, 0));
        int rowNo = 0;
        for (auto alpha : bnd_.knownDirs(bndNo)) {
            rhs[rowNo] = f(fieldNo, alpha, nodeNo) + 0.5 * DXQY::w[alpha] * DXQY::c2Inv * cF[alpha];
            rowNo++;
        }
        // Solve system
        auto x = svd_[bndNo].solve(rhs);

        // Calculate new solutions
        auto ct = DXQY::cDotAll(tangents(0, nodeNo));
        auto cn = DXQY::cDotAll(normals(0, nodeNo));
        auto q = qDist(0, nodeNo);

        for ( auto alpha : bnd_.unknownDirs(bndNo) ) {
            auto w = DXQY::w[alpha];
            // Boundry Size
            boundarySize += w;

            lbBase_t tmp = 0;
            tmp += x[0];
            // velocity
            tmp += -(q * DXQY::c4Inv / tau) * ct[alpha] * x[1];
            tmp += -(q * DXQY::c4Inv0_5 / tau) * cn[alpha] * x[2];
            // strain rate
            tmp += (DXQY::c4Inv   ) * (ct[alpha]*cn[alpha]) * x[1];
            tmp += (DXQY::c4Inv0_5) * (cn[alpha]*cn[alpha] - DXQY::c2) * x[2];

            f(fieldNo, alpha, nodeNo) = w * tmp - 0.5 * w * DXQY::c2Inv * cF[alpha];
        }

        // Mass conservation
        lbBase_t nodeMass = 0;
        for (int alpha = 0; alpha < DXQY::nQ; ++alpha)
            nodeMass += f(fieldNo, alpha, nodeNo);
        boundaryMass += nodeMass - 1;
    }

    // Mass conservation
    lbBase_t changeInWallMass = boundaryMass - boundaryMass_;
    lbBase_t changeInGlobalWallMass;
    lbBase_t changeInGlobalBulkMass;
    lbBase_t boundarySizeGlobal;
    
    MPI_Allreduce(&changeInWallMass, &changeInGlobalWallMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&changeInBulkMass, &changeInGlobalBulkMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&boundarySize, &boundarySizeGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    lbBase_t addMass = -(changeInGlobalWallMass + changeInGlobalBulkMass) / boundarySizeGlobal;

    // Update the boundary mass
    boundaryMass_ = boundaryMass + boundarySize * addMass;

    // Adding the lost mass to the boundary
    for (int bndNo = 0; bndNo < bnd_.size(); ++bndNo) {
        int nodeNo = bnd_.nodeNo(bndNo);
        for ( auto alpha: bnd_.unknownDirs(bndNo)) {
            f(fieldNo, alpha, nodeNo) += DXQY::w[alpha]*addMass;
        }
    }
}

#endif
