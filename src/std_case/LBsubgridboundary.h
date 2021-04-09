#ifndef LBSUBGRIDBOUNDARY_H
#define LBSUBGRIDBOUNDARY_H

#include <math.h>

#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBboundary.h"
#include "../lbsolver/LBgrid.h"
#include "../lbsolver/LBfield.h"

//  Linear  package
#include "../Eigen/Dense"
#include "../Eigen/SVD"

using namespace Eigen;

template <typename DXQY>
class OneNodeSubGridBnd : public BoundaryExtended<DXQY>
{
public:
    OneNodeSubGridBnd(const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const ScalarField &qDist, const VectorField<DXQY> &normals, const VectorField<DXQY> &tangents, const ScalarField &rho, const VectorField<DXQY> &force, lbBase_t tau);
    ~OneNodeSubGridBnd() {}
    void apply(const int fieldNo, LbField<DXQY> &f, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid);
private:
    template <typename MT>
    void matrixFillfa(const int row, MT & m, const int alpha);
    template <typename MT, typename T>
    void matrixFillSnn(const int row, MT & m, const lbBase_t qDist, T & normal, const lbBase_t un_wall, const lbBase_t Snn);
    template <typename MT, typename T>
    void matrixFillSnt(const int row, MT & m, const lbBase_t qDist, T & normal, T & tangent, const lbBase_t ut_wall,
                       const lbBase_t dun_dt, const lbBase_t Fn,const lbBase_t Ft, const lbBase_t rho0, const lbBase_t tau);
    template <typename MT, typename T>
    void matrixFillStt(const int row, MT & m, T & tangent, const lbBase_t Stt, const lbBase_t Ft, const lbBase_t tau);
    template <typename T>
    inline lbBase_t calcfalpha(const int alpha, const T &x) const;

    std::vector<BDCSVD<MatrixXd>> svd_;
    lbBase_t boundaryMass_;
};

template <typename DXQY>
OneNodeSubGridBnd<DXQY>::OneNodeSubGridBnd(const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const ScalarField &qDist, const VectorField<DXQY> &normals,  const VectorField<DXQY> &tangents, const ScalarField &rho,  const VectorField<DXQY> &force, const lbBase_t tau)
    : BoundaryExtended<DXQY>(bndNodes, nodes, grid), svd_(bndNodes.size())
{
    MatrixXd m;
    int nCol = 1 + DXQY::nD + (DXQY::nD*(DXQY::nD + 1))/2;
    // Mass at boundary
    boundaryMass_ = 0;
    for (int n=0; n < this->size(); ++n) {
        int nRow = (1 + this->nBeta(n) + 2*(this->nGamma(n))) + 3;
        // Setup matrix
        m.resize(nRow, nCol);
        int row = 0;
        // -- Setup equations from fa
        // -- -- beta_reversed
        for (auto beta: this->beta(n)) {
            int betaRev = this->dirRev(beta);
            matrixFillfa(row++, m, betaRev);
        }
        // -- -- gamma and gamma_reversed
        for (auto gamma: this->gamma(n)) {
            int gammaRev = this->dirRev(gamma);
            matrixFillfa(row++, m, gamma);
            matrixFillfa(row++, m, gammaRev);
        }
        // -- -- zero velocity
        int zeroVelDir = DXQY::nQ-1;
        matrixFillfa(row++, m, zeroVelDir);
        int nodeNo = this->nodeNo(n);
        // Setup equations from Snn
        lbBase_t q = qDist(0, nodeNo);
        std::valarray<lbBase_t> nVec = normals(0, nodeNo);
        lbBase_t unWall = 0.0;
        lbBase_t snn = 0.0;
        matrixFillSnn(row++, m, q, nVec, unWall, snn);
        // Setup equations from Snt
        std::valarray<lbBase_t> tVec = tangents(0, nodeNo);
        lbBase_t utWall = 0.0;
        lbBase_t dundt = 0.0;
        int nforce = 0;
        if (force.getNumNodes() > 1)
            nforce = nodeNo;
        std::valarray<lbBase_t> F = force(0, nforce);
        lbBase_t Fn = 0;
        lbBase_t Ft = 0;
        for (int i = 0; i < DXQY::nD; ++i) {
            Fn += nVec[i]*F[i];
            Ft += tVec[i]*F[i];
        }
        lbBase_t rho0 = 1.0;
        matrixFillSnt(row++, m, q, nVec, tVec, utWall, dundt, Fn, Ft, rho0, tau);
        // Setup equations frmo Stt
        lbBase_t Stt = 0;
        matrixFillStt(row++, m, tVec, Stt, Ft, tau);
        // Setup Singular Value Decomposition
        svd_[n] = svd_[n].compute(m, ComputeThinU | ComputeThinV);
        // Check that the boundary conditions has enough data to solve the system.
        if (svd_[n].rank() < nCol) {
            std::cout << "Error in the boundary conditions: rank (" << svd_[n].rank() << ") is less than number of unknowns (" << nCol << ")." << std::endl;
        }
        // Find the total mass
        boundaryMass_ += rho(0, nodeNo) - 1;
    }
}

template <typename DXQY>
template <typename MT>
void OneNodeSubGridBnd<DXQY>::matrixFillfa(const int row, MT &m, const int a)
{
    int col = 0;
    m(row, col++) = DXQY::w[a];  // rho
    for (int i=0; i < DXQY::nD; ++i) {  // j_i
        m(row, col++) = DXQY::w[a] * DXQY::c(a, i) * DXQY::c2Inv;
    }
    for (int i=0; i < DXQY::nD; ++i) {  // Pi_ij
        m(row, col++) = DXQY::w[a] * (DXQY::c(a, i)*DXQY::c(a, i) - DXQY::c2) * DXQY::c4Inv0_5;
        for (int j = i+1; j < DXQY::nD; ++j) {
            m(row, col++) = DXQY::w[a] * DXQY::c(a, i)*DXQY::c(a, j) * DXQY::c4Inv;
        }
    }
}

template <typename DXQY>
template <typename MT, typename T>
void OneNodeSubGridBnd<DXQY>::matrixFillSnn(const int row, MT &m, const lbBase_t q, T & n, const lbBase_t un, const lbBase_t Snn)
{
    int col = 0;
    m(row, col++) = q*Snn + un;
    for (int i = 0; i < DXQY::nD; ++i) {
        m(row, col++) = -n[i];
    }
    for (int i = 0; i < DXQY::nD; ++i) {
        for (int j = i; j < DXQY::nD; ++j) {
            m(row, col++) = 0;
        }
    }
}

template <typename DXQY>
template <typename MT, typename T>
void OneNodeSubGridBnd<DXQY>::matrixFillSnt(const int row, MT & m, const lbBase_t q, T & n, T & t, const lbBase_t ut,
        const lbBase_t dudt, const lbBase_t Fn, const lbBase_t Ft, const lbBase_t rho0, const lbBase_t tau)
{
    int col = 0;
    m(row, col++) = q*dudt - ut;
    lbBase_t k = 1.0/(2*DXQY::c2*tau*rho0);
    for (int i = 0; i < DXQY::nD; ++i) {
        m(row, col++) = t[i] + q*k*(t[i]*Fn + n[i]*Ft);
    }
    lbBase_t k2 = DXQY::c2Inv/tau;
    for (int i = 0; i < DXQY::nD; ++i) {
        m(row, col++) = k2*q*n[i]*t[i];
        for (int j = i+1; j < DXQY::nD; ++j) {
            m(row, col++) = q*k2*(t[i]*n[j] + t[j]*n[i]);
        }
    }
}

template <typename DXQY>

template <typename MT, typename T>
void OneNodeSubGridBnd<DXQY>::matrixFillStt(const int row, MT & m, T & t, const lbBase_t Stt,
        const lbBase_t Ft, const lbBase_t tau)
{
    int col = 0;
    m(row, col++) = Stt;
    lbBase_t k1 = 0.5*DXQY::c2Inv/tau;
    for (int i=0; i<DXQY::nD; ++i) {
        m(row, col++) = k1*t[i]*Ft;
    }
    lbBase_t k2 = DXQY::c2Inv/tau;
    for (int i=0; i<DXQY::nD; ++i) {
        m(row, col++) = k1*t[i]*t[i];
        for (int j=i+1; j<DXQY::nD; ++j) {
            m(row, col++) = k2*t[i]*t[j];
        }
    }
}

template <typename DXQY>
template <typename T>
inline lbBase_t OneNodeSubGridBnd<DXQY>::calcfalpha(const int a, const T &x) const
{
    lbBase_t ret = 0;
    int n = 0;
    ret += x[n++]; // Rho
    for (int i = 0; i < DXQY::nD; ++i) {  // \rho u
        ret += DXQY::c2Inv*DXQY::c(a, i)*x[n++];
    }
    for (int i = 0; i < DXQY::nD; ++i) {
        ret += DXQY::c4Inv0_5*(DXQY::c(a, i)*DXQY::c(a, i)-DXQY::c2)*x[n++];
        for (int j = i+1; j < DXQY::nD; ++j) {
            ret += DXQY::c4Inv0_5*DXQY::c(a, i)*DXQY::c(a, j)*x[n++];
        }
    }
    return D2Q9::w[a]*ret;
}

template <typename DXQY>
void OneNodeSubGridBnd<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
{
    lbBase_t predictedBoundaryMass = 0;
    lbBase_t deltaBulkMass = 0;
    for (int n = 0; n < this->size(); ++n) {
        int nRows = svd_[n].rows();
        VectorXd rhs(nRows);
        // Set the right hand values. Rememver to go through directions in the same order
        // as for the matrix.
        int row = 0;
        int nodeNo = this->nodeNo(n);
        // Known directions
        // -- -- beta_reversed
        for (auto beta: this->beta(n)) {
            int betaRev = this->dirRev(beta);
            rhs[row++] = f(fieldNo, betaRev, nodeNo);
            int nodeNeigNo = grid.neighbor(beta, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo)) {
                deltaBulkMass += f(fieldNo, beta, nodeNeigNo) - f(fieldNo, betaRev, nodeNo);
            }
        }
        // -- -- gamma and gamma_reversed
        for (auto gamma: this->gamma(n)) {
            int gammaRev = this->dirRev(gamma);
            rhs[row++] = f(fieldNo, gamma, nodeNo);
            rhs[row++] = f(fieldNo, gammaRev, nodeNo);
            int nodeNeigNo = grid.neighbor(gammaRev, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo))
                deltaBulkMass += f(fieldNo, gammaRev, nodeNeigNo) - f(fieldNo, gamma, nodeNo);
            nodeNeigNo = grid.neighbor(gamma, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo))
                deltaBulkMass += f(fieldNo, gamma, nodeNeigNo) - f(fieldNo, gammaRev, nodeNo);
        }
        // -- -- zero velocity
        int zeroVelDir = DXQY::nQ-1;
        rhs[row++] = f(fieldNo, zeroVelDir, nodeNo);
        // Fill the remaining elements with zeros
        for (; row < nRows;)
            rhs[row++] = 0;
        // Solve system of equations
        VectorXd x = svd_[n].solve(rhs);
        predictedBoundaryMass += x[0] - 1;
        // Set new distributions
        for (int q = 0; q < DXQY::nQ; ++q) {
            f(fieldNo, q, nodeNo) = calcfalpha(q, x);
        }
    }
    lbBase_t deltaWallMass = predictedBoundaryMass - boundaryMass_;
    lbBase_t deltaWallMassGlobal;
    lbBase_t deltaBulkMassGlobal;
    lbBase_t boundarySize = this->size(); 
    lbBase_t boundarySizeGlobal;

    MPI_Allreduce(&deltaWallMass, &deltaWallMassGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&deltaBulkMass, &deltaBulkMassGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&boundarySize, &boundarySizeGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    lbBase_t addMass = -(deltaWallMassGlobal + deltaBulkMassGlobal) / boundarySizeGlobal;
    
    // Update the boundary mass
    boundaryMass_ = predictedBoundaryMass + this->size() * addMass;

    for (int n = 0; n < this->size(); ++n) {
        int nodeNo = this->nodeNo(n);
        for (int q = 0; q < DXQY::nQ; ++q) {
            f(fieldNo, q, nodeNo) += DXQY::w[q]*addMass;
        }
    }
}


#endif
