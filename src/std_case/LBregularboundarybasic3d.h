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

/* template<typename DXQY>
Eigen::Matrix<float, DXQY::nD, DXQY::nD> matrixConservationMomentFull()
{
    Eigen::Matrix<float, DXQY::nD, DXQY::nD> m;
    for (int i=0; i < DXQY::nD; ++i) {
        for (int j=0; j < DXQY::nD; ++j) {
            m(i, j) = 0;
            for (int q=0; q<DXQY::nQ; ++q) {
                m(i, j) = DXQY::w[q]*DXQY::c(q, i)*DXQY::c(q, j);
            }
        }
    }
} */

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
    typedef Eigen::Matrix<double, DXQY::nD, DXQY::nD> MatrixCons;
    typedef Eigen::MatrixXd MatrixReg;
    typedef std::vector<lbBase_t> lbVec;

    MatrixCons setupMatrixConservation_();
    template <typename L>
    MatrixCons setupMatrixConservation_(const L &unknownDirections);
    template <typename Tv, typename Tq> 
    MatrixReg setupMatrixRegular_(
        const lbBase_t &q, 
        const Tv &cn, 
        const Tv &ct1, 
        const Tv &ct2, 
        const lbBase_t &tau, 
        const Tq &knownDirections
    );
    template <typename Tf, typename TcF, typename TkD>
    lbVec setRightHandSideRegular(
    const Tf  &f,
    const TcF &cF,
    const TkD &knownDirections
    );
    template <typename Tv, typename Tm>
    lbVec regularsDist(
        const lbBase_t &q, 
        const Tv &cn, 
        const Tv &ct1, 
        const Tv &ct2, 
        const Tv &cF,
        const lbBase_t &tau, 
        const Tm &macroVars
    );
    template <typename T1, typename T2, typename T3>
    lbVec setRightHandSideConservation(
        const T1 &f,
        const T2 &freg,
        const T3 &knownDirections
    );
    template <typename T1, typename T2>
    lbVec setRightHandSideConservation(
        const T1 &f,
        const T2 &freg
    );


    BoundaryBasic<DXQY> bnd_;
    std::vector<Eigen::BDCSVD<MatrixReg>> matrixRegularSolution_;
    std::vector<Eigen::BDCSVD<MatrixReg>> matrixConservation_;
    std::vector<Eigen::BDCSVD<MatrixReg>> matrixConservationAll_;
    std::vector<bool> solveMatrixConservationAll_;
    lbBase_t boundaryMass_;
    Eigen::Matrix3f delta_f_;
};

/********************************
 *      Private functions
 ********************************/
template<typename DXQY>
typename RegularBoundaryBasic3d<DXQY>::MatrixCons RegularBoundaryBasic3d<DXQY>::setupMatrixConservation_()
{
    MatrixCons m;
    for (int i=0; i < DXQY::nD; ++i) {
        for (int j=0; j < DXQY::nD; ++j) {
            m(i, j) = 0;
            for (int q=0; q<DXQY::nQ; ++q) {
                m(i, j) = DXQY::w[q]*DXQY::c(q, i)*DXQY::c(q, j);
            }
        }
    }
    return m;
}

template <typename DXQY>
template <typename L>
typename RegularBoundaryBasic3d<DXQY>::MatrixCons RegularBoundaryBasic3d<DXQY>::setupMatrixConservation_(const L &unknownDirections)
{
    MatrixCons m;
    for (int i=0; i < DXQY::nD; ++i) {
        for (int j=0; j < DXQY::nD; ++j) {
            m(i, j) = 0;
            for (auto q: unknownDirections) {
                m(i, j) = DXQY::w[q]*DXQY::c(q, i)*DXQY::c(q, j);
            }
        }
    }
    return m;
}

template <typename DXQY>
template <typename Tv, typename Tq> 
typename RegularBoundaryBasic3d<DXQY>::MatrixReg RegularBoundaryBasic3d<DXQY>::setupMatrixRegular_(
        const lbBase_t &q, 
        const Tv &cn, 
        const Tv &ct1, 
        const Tv &ct2, 
        const lbBase_t &tau, 
        const Tq &knownDirections
)
{
    int nRows = knownDirections.size();
    int nColumns = 4;

    //  Setup matrix
    MatrixReg m(nRows, nColumns);
    int rowNo = 0;
    for (auto alpha : knownDirections) {
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

    return m;
}

template <typename DXQY>
template <typename Tf, typename TcF, typename TkD>
typename RegularBoundaryBasic3d<DXQY>::lbVec RegularBoundaryBasic3d<DXQY>::setRightHandSideRegular(
    const Tf  &f,
    const TcF &cF,
    const TkD &knownDirections
)
{
    // Setup the right hand side vector
    auto nRows = knownDirections.size();
    lbVec rhs(nRows);
    int rowNo = 0;
    for (auto alpha : knownDirections) {
            rhs[rowNo] = f[alpha] + 0.5 * DXQY::w[alpha] * DXQY::c2Inv * cF[alpha];
            rowNo++;            
    }
    return rhs;
}


template <typename DXQY>
template <typename Tv, typename Tm>
typename RegularBoundaryBasic3d<DXQY>::lbVec RegularBoundaryBasic3d<DXQY>::regularsDist(
    const lbBase_t &q, 
    const Tv &cn, 
    const Tv &ct1, 
    const Tv &ct2, 
    const Tv &cF,
    const lbBase_t &tau, 
    const Tm &macroVars
)
{
    lbVec fReg(DXQY::nQ);
    for (int alpha = 0; alpha < DXQY::nQ; ++alpha) {
        auto w = DXQY::w[alpha];

        lbBase_t tmp = 0;
        tmp += macroVars[0];
        // velocity
        tmp += -(q * DXQY::c4Inv / tau) * ct1[alpha] * macroVars[1];
        tmp += -(q * DXQY::c4Inv / tau) * ct2[alpha] * macroVars[2];
        tmp += -(q * DXQY::c4Inv0_5 / tau) * cn[alpha] * macroVars[3];
        // strain rate
        tmp += (DXQY::c4Inv   ) * (ct1[alpha]*cn[alpha]) * macroVars[1];
        tmp += (DXQY::c4Inv   ) * (ct2[alpha]*cn[alpha]) * macroVars[2];
        tmp += (DXQY::c4Inv0_5) * (cn[alpha]*cn[alpha] - DXQY::c2) * macroVars[3];
        
        fReg[alpha] = w * tmp - 0.5 * w * DXQY::c2Inv * cF[alpha];            
    }
    return fReg;
}

template <typename DXQY>
template <typename T1, typename T2, typename T3>
typename RegularBoundaryBasic3d<DXQY>::lbVec RegularBoundaryBasic3d<DXQY>::setRightHandSideConservation(
    const T1 &f,
    const T2 &freg,
    const T3 &knownDirections
)
{
    lbVec rhs(DXQY::nD);

    for(int d = 0; d < DXQY::nD; ++d) {
        rhs[d] = 0;
        for (auto alpha: knownDirections) {
            rhs[d] += (freg[alpha] - f[alpha])*DXQY::c(alpha, d); 
        }
    } 

    return rhs;
}

template <typename DXQY>
template <typename T1, typename T2>
typename RegularBoundaryBasic3d<DXQY>::lbVec RegularBoundaryBasic3d<DXQY>::setRightHandSideConservation(
    const T1 &f,
    const T2 &freg
)
{
    lbVec rhs(DXQY::nD);

    for(int d = 0; d < DXQY::nD; ++d) {
        rhs[d] = 0;
        for (int alpha=0; alpha < DXQY::nQ; ++alpha) {
            rhs[d] += (freg[alpha] - f[alpha])*DXQY::c(alpha, d); 
        }
    } 

    return rhs;
}


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
    matrixRegularSolution_(boundaryNodes.size()),
    matrixConservation_(boundaryNodes.size()),
    matrixConservationAll_(boundaryNodes.size()),
    solveMatrixConservationAll_(boundaryNodes.size())
{
    // Initiation for mass conservation
    lbBase_t boudaryMass_ = 0;
    
    // Setup the SVD matrices
    for (int bndNo = 0; bndNo < bnd_.size(); ++bndNo) {
        auto nodeNo = bnd_.nodeNo(bndNo);
        // Mass conservation
        boundaryMass_ += rho(0, nodeNo) - 1;

        auto ct1 = DXQY::cDotAll(tangents(0, nodeNo));
        auto ct2 = DXQY::cDotAll(tangents(1, nodeNo));
        auto cn = DXQY::cDotAll(normals(0, nodeNo));
        auto q = qDist(0, nodeNo);
        
       /*
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
        }  */    
 
        // Setup Singular Value Decomposition
        auto m = setupMatrixRegular_(q, cn, ct1, ct2, tau, bnd_.knownDirs(bndNo));
        matrixRegularSolution_[bndNo] = matrixRegularSolution_[bndNo].compute(m, Eigen::ComputeThinU | Eigen::ComputeThinV);

        // Check that the boundary conditions has enough data to solve the system.
        if (matrixRegularSolution_[bndNo].rank() < m.cols()) {
            std::cout << "Error in the boundary conditions: rank (" << matrixRegularSolution_[bndNo].rank() << ") is less than number of unknowns (" << m.cols() << ")." << std::endl;
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
                
        // Impuls conservation
        auto m1 = setupMatrixConservation_(bnd_.unknownDirs(bndNo));
        matrixConservation_[bndNo] = matrixConservation_[bndNo].compute(m1, Eigen::ComputeThinU | Eigen::ComputeThinV);
        solveMatrixConservationAll_[bndNo] = false;
        if (matrixConservation_[bndNo].rank() < m1.cols())
        {
            auto m2 = setupMatrixConservation_();
            matrixConservationAll_[bndNo] = matrixConservationAll_[bndNo].compute(m2, Eigen::ComputeThinU | Eigen::ComputeThinV);
            solveMatrixConservationAll_[bndNo] = true;
        }
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
        auto nodeNo = bnd_.nodeNo(bndNo);

        // Change in bulk mass
        for (int alpha = 0; alpha < DXQY::nQ; ++alpha) {                        
            auto alphaRev = bnd_.dirRev(alpha);
            auto nodeNeigNo = grid.neighbor(alphaRev, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo)) {
                changeInBulkMass += f(fieldNo, alphaRev, nodeNeigNo) - f(fieldNo, alpha, nodeNo);
            }            
        }

        /* // Setup right hand side and solve the system of equations
        // Setup the right hand side vector
        auto nRows = bnd_.nKnownDirs(bndNo);
        Eigen::VectorXd rhs(nRows);
        auto cF = DXQY::cDotAll(force(0, 0));
        int rowNo = 0;
        for (auto alpha : bnd_.knownDirs(bndNo)) {
                rhs[rowNo] = f(fieldNo, alpha, nodeNo) + 0.5 * DXQY::w[alpha] * DXQY::c2Inv * cF[alpha];
                rowNo++;            
        } */


        // Calculate regular solutions
        auto q = qDist(0, nodeNo);
        auto cn = DXQY::cDotAll(normals(0, nodeNo));
        auto ct1 = DXQY::cDotAll(tangents(0, nodeNo));
        auto ct2 = DXQY::cDotAll(tangents(1, nodeNo));
        auto cF = DXQY::cDotAll(force(0, 0));

        auto rhs = setRightHandSideRegular(cF, bnd_.knownDirs(bndNo));
        auto x = matrixRegularSolution_[bndNo].solve(rhs);

        auto fReg = regularsDist(q, cn, ct1, ct2, cF, tau, x);

         // Calculate the right hand side for conservation
         auto df = matrixConservation_[bndNo].solve(setRightHandSideConservation(f(fieldNo, nodeNo), fReg, bnd_.knownDirs(bndNo)));
         for (auto alpha: bnd_.unknownDirs(bndNo)){
             f(fieldNo, alpha, nodeNo) = fReg[alpha];
             for (int d=0; d < DXQY::nD; ++d) {
                f(fieldNo, alpha, nodeNo) += DXQY::w[alpha]*DXQY::c(alpha, d)*df(d);
             }
         }

         if (solveMatrixConservationAll_[bndNo]) 
         {
             auto dfAll = matrixConservationAll_[bndNo].solve(setRightHandSideConservation(f(fieldNo, nodeNo)));
            for (int alpha=0; alpha < DXQY::nD; ++alpha){
                f(fieldNo, alpha, nodeNo) = fReg[alpha];
                for (int d=0; d < DXQY::nD; ++d) {
                    f(fieldNo, alpha, nodeNo) += DXQY::w[alpha]*DXQY::c(alpha, d)*dfAll(d);
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
