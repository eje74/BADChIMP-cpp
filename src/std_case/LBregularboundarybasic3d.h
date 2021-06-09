#ifndef LBREGULARBOUNDARYBASIC3D_H
#define LBREGULARBOUNDARYBASIC3D_H

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
class RegularBoundaryBasic3d
{
public:
    RegularBoundaryBasic3d(
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
    std::vector<Eigen::BDCSVD<Eigen::MatrixXd>> svd_;
    lbBase_t boundaryMass_;
    Eigen::Matrix3f delta_f_;
};



/********************************
 *      class constructor
 ********************************/
template <typename DXQY>
RegularBoundaryBasic3d<DXQY>::RegularBoundaryBasic3d(
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
        // Surface data
        auto ct1 = DXQY::cDotAll(tangents(0, nodeNo));
        auto ct2 = DXQY::cDotAll(tangents(1, nodeNo));
        auto cn = DXQY::cDotAll(normals(0, nodeNo));
        auto q = qDist(0, nodeNo);
        // Matrix size
        int nRows = bnd_.nKnownDirs(bndNo);
        int nColumns = 4;

        //  Setup matrix
        Eigen::MatrixXd m(nRows, nColumns);
        int rowNo = 0;
        for (auto alpha : bnd_.knownDirs(bndNo)) {
                auto w = DXQY::w[alpha];
                for (int j = 0; j < nColumns; ++j)
                    m(rowNo, j) = 0;
                // rho
                m(rowNo, 0) = w;
                // velocity
                m(rowNo, 1) += -(w * q * DXQY::c4Inv / tau) * ct1[alpha];
                m(rowNo, 2) += -(w * q * DXQY::c4Inv / tau) * ct2[alpha];
                m(rowNo, 3) += -(w * q * DXQY::c4Inv0_5 / tau) * cn[alpha];
                // strain rate
                m(rowNo, 1) += (w * DXQY::c4Inv   ) * (ct1[alpha]*cn[alpha]);
                m(rowNo, 2) += (w * DXQY::c4Inv   ) * (ct2[alpha]*cn[alpha]);
                m(rowNo, 3) += (w * DXQY::c4Inv0_5) * (cn[alpha]*cn[alpha] - DXQY::c2);
                rowNo++;            
        }     
 
        // Setup Singular Value Decomposition
        svd_[bndNo] = svd_[bndNo].compute(m, Eigen::ComputeThinU | Eigen::ComputeThinV);

        // Check that the boundary conditions has enough data to solve the system.
        if (svd_[bndNo].rank() < nColumns) {
            std::cout << "Error in the boundary conditions: rank (" << svd_[bndNo].rank() << ") is less than number of unknowns (" << nColumns << ")." << std::endl;
            std::cout << m << std::endl;
            for (int q = 0; q < DXQY::nQ; ++q)
                std::cout << cn[q] << " ";
            std::cout << std::endl;
            for (int q = 0; q < DXQY::nQ; ++q)
                std::cout << ct1[q] << " ";
            std::cout << std::endl;
            for (int q = 0; q < DXQY::nQ; ++q)
                std::cout << ct2[q] << " ";
            std::cout << std::endl;
            std::cout << q << std::endl;
            std::cout << tau << std::endl;
            exit(1);
        }
        
        // Mass conservation
        boundaryMass_ += rho(0, nodeNo) - 1;
        
        // Setup matrix conseravtion of the first moment
        Eigen::Matrix3f wcc;
        wcc << 0,0,0, 0,0,0, 0,0,0;
        if (bnd_.nUnknownDirs(bndNo) > 1) {
            for (auto alpha: bnd_.unknowDirs(bndNo)) {
                auto w = DXQY::w[alpha];
                auto c = DXQY::c(alpha);
                for (int i=0; i < DXQY::nD; ++i) {
                    for (int j=0; j < DXQY::nD; ++j) {
                        wcc(i, j) += w*c[i]*c[j];
                    }
                }                
            }
        }
        // Inverse of wcc
        delta_f_ = wcc.inverse();
        /* EJ: BEGIN HERE */
    }
}



/****************************************
 * Apply boundary condition
 * **************************************/
template <typename DXQY>
void RegularBoundaryBasic3d<DXQY>::apply(
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
    lbBase_t boundarySize = bnd_.size(); 

    // No slip boundary condition
    for (int bndNo = 0; bndNo < bnd_.size(); bndNo++) {
        // Setup the right hand side vector
        auto nRows = bnd_.nKnownDirs(bndNo);
        Eigen::VectorXd rhs(nRows);

        auto nodeNo = bnd_.nodeNo(bndNo);

        // Change in bulk mass
        for (int alpha = 0; alpha < DXQY::nQ; ++alpha) {                        
            auto alphaRev = bnd_.dirRev(alpha);
            auto nodeNeigNo = grid.neighbor(alphaRev, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo)) {
                changeInBulkMass += f(fieldNo, alphaRev, nodeNeigNo) - f(fieldNo, alpha, nodeNo);
            }            
        }

        // Setup right hand side and solve the system of equations
        auto cF = DXQY::cDotAll(force(0, 0));
        int rowNo = 0;
        for (auto alpha : bnd_.knownDirs(bndNo)) {
                rhs[rowNo] = f(fieldNo, alpha, nodeNo) + 0.5 * DXQY::w[alpha] * DXQY::c2Inv * cF[alpha];
                rowNo++;            
        }

        auto x = svd_[bndNo].solve(rhs);

        // Calculate regular solutions
        auto ct1 = DXQY::cDotAll(tangents(0, nodeNo));
        auto ct2 = DXQY::cDotAll(tangents(1, nodeNo));
        auto cn = DXQY::cDotAll(normals(0, nodeNo));
        auto q = qDist(0, nodeNo);

        std::vector<lbBase_t> fReg(DXQY::nQ);
        for ( int alpha = 0; alpha < DXQY::nQ; ++alpha) {
           
            auto w = DXQY::w[alpha];

            lbBase_t tmp = 0;
            tmp += x[0];
            // velocity
            tmp += -(q * DXQY::c4Inv / tau) * ct1[alpha] * x[1];
            tmp += -(q * DXQY::c4Inv / tau) * ct2[alpha] * x[2];
            tmp += -(q * DXQY::c4Inv0_5 / tau) * cn[alpha] * x[3];
            // strain rate
            tmp += (DXQY::c4Inv   ) * (ct1[alpha]*cn[alpha]) * x[1];
            tmp += (DXQY::c4Inv   ) * (ct2[alpha]*cn[alpha]) * x[2];
            tmp += (DXQY::c4Inv0_5) * (cn[alpha]*cn[alpha] - DXQY::c2) * x[3];
            
            fReg[alpha] = w * tmp - 0.5 * w * DXQY::c2Inv * cF[alpha];            
         }
         
         // Find the perturbation from the regular solutions
         std::vector<lbBase_t> df(DXQY::nD, 0);
         for ( auto alpha: bnd_.knownDirs(bndNo) ) 
         {
            for (int d = 0; d  < DXQY::nD; ++d) 
            {
                df[d] += (f(fieldNo, alpha, nodeNo) - fReg[alpha])*DXQY::c(alpha, d);
            }
         }         
         
         if ( bnd_.nUnknownDirs(bndNo) > 1 )
         {
             lbBase_t a11 = 0, a12 = 0, a21 = 0, a22 = 0;
             for (auto alpha : bnd_.unknownDirs(bndNo)) 
             {
                 auto w = DXQY::w[alpha];
                 auto c = DXQY::c(alpha);
                 a11 += w*c[0]*c[0];
                 a12 += w*c[0]*c[1];
                 a21 += w*c[1]*c[0];
                 a22 += w*c[1]*c[1];
             }
                  
            std::vector<lbBase_t> deltaF(2);
            lbBase_t det = a11*a22 - a12*a21;
            deltaF[0] = (-df[0]*a22 + df[1]*a12)/det;
            deltaF[1] = (-df[1]*a11 + df[0]*a21)/det;

            for (auto alpha : bnd_.unknownDirs(bndNo)) 
            {
                 auto w = DXQY::w[alpha];
                 auto c = DXQY::c(alpha);             
                f(fieldNo, alpha, nodeNo) = fReg[alpha] + w*c[0]*deltaF[0] + w*c[1]*deltaF[1];            
            }
         } else { // One unknown

             // Minimize the error using the unknow direction
             lbBase_t dfUnknowndir = 0;
             auto alphas = bnd_.unknownDirs(bndNo);
             for (int d = 0; d < DXQY::nD; ++d)
                 dfUnknowndir -= df[d]*DXQY::c(alphas[0], d);
             dfUnknowndir /= DXQY::cNorm[alphas[0]]*DXQY::cNorm[alphas[0]];
             f(fieldNo, alphas[0], nodeNo) = fReg[alphas[0]] + dfUnknowndir;
             
             // Find the perturbation from the regular solutions
            df[0] += dfUnknowndir * DXQY::c(alphas[0], 0);
            df[1] += dfUnknowndir * DXQY::c(alphas[0], 1);         
        

             lbBase_t a11 = 0, a12 = 0, a21 = 0, a22 = 0;
             for ( int alpha = 0; alpha < DXQY::nQ; ++alpha ) 
             {
                auto alphaRev = DXQY::reverseDirection(alpha);
                auto neigNode = grid.neighbor(alphaRev, nodeNo);
                if ( nodes.isFluidBoundary(neigNode) || nodes.isSolidBoundary(neigNode) )
                 {
                     auto w = DXQY::w[alpha];
                     auto c = DXQY::c(alpha);
                     a11 += w*c[0]*c[0];
                     a12 += w*c[0]*c[1];
                     a21 += w*c[1]*c[0];
                     a22 += w*c[1]*c[1];                     
                 }
             }
                  
            std::vector<lbBase_t> deltaF(2);
            lbBase_t det = a11*a22 - a12*a21;
            deltaF[0] = (-df[0]*a22 + df[1]*a12)/det;
            deltaF[1] = (-df[1]*a11 + df[0]*a21)/det;

            for ( int alpha = 0; alpha < DXQY::nQ; ++alpha ) 
            {
                 auto alphaRev = DXQY::reverseDirection(alpha);
                 auto neigNode = grid.neighbor(alphaRev, nodeNo);
                 if ( nodes.isFluidBoundary(neigNode) || nodes.isSolidBoundary(neigNode) )
                 {
                     auto w = DXQY::w[alpha];
                     auto c = DXQY::c(alpha);             
                     f(fieldNo, alpha, nodeNo) += w*c[0]*deltaF[0] + w*c[1]*deltaF[1];
                 }            
            }          
         }
        
        // Update the boundary mass
        lbBase_t nodeMass = 0;
        // Mass conservation
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

    for (int bndNo = 0; bndNo < bnd_.size(); ++bndNo) 
    {
        int nodeNo = bnd_.nodeNo(bndNo);
        for ( int alpha = 0; alpha < DXQY::nQ; ++alpha ) {
            f(fieldNo, alpha, nodeNo) += DXQY::w[alpha]*addMass;
        }
    }
}


#endif
