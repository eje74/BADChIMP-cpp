#ifndef LBREGULARBOUNDARYBASIC_H
#define LBREGULARBOUNDARYBASIC_H

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
class RegularBoundaryBasic
{
public:
    RegularBoundaryBasic(
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

private:
    BoundaryBasic<DXQY> bnd_;
    std::vector<BDCSVD<MatrixXd>> svd_;
    lbBase_t boundaryMass_;
};



/********************************
 *      class constructor
 ********************************/
template <typename DXQY>
RegularBoundaryBasic<DXQY>::RegularBoundaryBasic(
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
    svd_(boundaryNodes.size())
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
            m(rowNo, 0) = w; // rho
            m(rowNo, 1) = ( w * DXQY::c4Inv ) * ( cn[alpha]*ct[alpha] - q*ct[alpha]/tau );
            m(rowNo, 2) = w * DXQY::c4Inv0_5 * ( ct[alpha]*ct[alpha] - DXQY::c2 );
            rowNo++;
        }

        // Setup Singular Value Decomposition
        svd_[bndNo] = svd_[bndNo].compute(m, ComputeThinU | ComputeThinV);

        // Check that the boundary conditions has enough data to solve the system.
        if (svd_[bndNo].rank() < nColumns) {
            std::cout << "Error in the boundary conditions: rank (" << svd_[bndNo].rank() << ") is less than number of unknowns (" << nColumns << ")." << std::endl;
        }
        
        // Mass conservation
        boundaryMass_ += rho(0, nodeNo) - 1;
    }
}



/****************************************
 * Apply boundary condition
 * **************************************/
template <typename DXQY>
void RegularBoundaryBasic<DXQY>::apply(
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
    
    // No slip boundary condition
    for (int bndNo = 0; bndNo < bnd_.size(); bndNo++) {
        // Setup the right hand side vector
        auto nRows = bnd_.nKnownDirs(bndNo);
        VectorXd rhs(nRows);

        auto nodeNo = bnd_.nodeNo(bndNo);
        auto cF = DXQY::cDotAll(force(0, 0));
        int rowNo = 0;
        for (auto alpha : bnd_.knownDirs(bndNo)) {
            rhs[rowNo] = f(fieldNo, alpha, nodeNo) + 0.5 * DXQY::w[alpha] * DXQY::c2Inv * cF[alpha];
            rowNo++;
            
            // Mass conservation
            auto alphaRev = bnd_.dirRev(alpha);
            auto nodeNeigNo = grid.neighbor(alphaRev, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo)) {
                changeInBulkMass += f(fieldNo, alphaRev, nodeNeigNo) - f(fieldNo, alpha, nodeNo);
            }
        }

        // Solve system
        auto x = svd_[bndNo].solve(rhs);

        // Calculate new solutions
        auto ct = DXQY::cDotAll(tangents(0, nodeNo));
        auto cn = DXQY::cDotAll(normals(0, nodeNo));
        auto q = qDist(0, nodeNo);

        for (int alpha = 0; alpha < DXQY::nQ; ++alpha) {
            auto w = DXQY::w[alpha];
            lbBase_t tmp = 0;

            tmp += x[0];
            tmp += DXQY::c4Inv  * ( cn[alpha]*ct[alpha] - q*ct[alpha]/tau ) * x[1];
            tmp += DXQY::c4Inv0_5 * ( ct[alpha]*ct[alpha] - DXQY::c2 ) * x[2];
            f(fieldNo, alpha, nodeNo) = w * tmp - 0.5 * w * DXQY::c2Inv * cF[alpha];
        }
        
        // Mass conservation
        boundaryMass += x[0] - 1;
    }
    
    // Mass conservation
    lbBase_t changeInWallMass = boundaryMass - boundaryMass_;
    lbBase_t changeInGlobalWallMass;
    lbBase_t changeInGlobalBulkMass;
    lbBase_t boundarySize = bnd_.size(); 
    lbBase_t boundarySizeGlobal;

    MPI_Allreduce(&changeInWallMass, &changeInGlobalWallMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&changeInBulkMass, &changeInGlobalBulkMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&boundarySize, &boundarySizeGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    lbBase_t addMass = -(changeInGlobalWallMass + changeInGlobalBulkMass) / boundarySizeGlobal;
    
    // Update the boundary mass
    boundaryMass_ = boundaryMass + bnd_.size() * addMass;

    for (int n = 0; n < bnd_.size(); ++n) {
        int nodeNo = bnd_.nodeNo(n);
        for (int q = 0; q < DXQY::nQ; ++q) {
            f(fieldNo, q, nodeNo) += 0*DXQY::w[q]*addMass;
        }
    }
}


#endif
