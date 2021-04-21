#ifndef LBSUBGRIDBOUNDARYD_H
#define LBSUBGRIDBOUNDARYD_H

/***************************************************
 * The dynamic version of of the regularized boundary condition
 *
 * We assume that the
 *
 ***************************************************/

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
class OneNodeSubGridBndDyn : public BoundaryExtended<DXQY>
{
public:
    OneNodeSubGridBndDyn(const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const ScalarField &qDist, const VectorField<DXQY> &normals,
                         const VectorField<DXQY> &tangents, const ScalarField &rho, const VectorField<DXQY> &force, lbBase_t tau);
    ~OneNodeSubGridBndDyn() {}
    void apply(const int fieldNo, LbField<DXQY> &f, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const ScalarField &qDist, const VectorField<DXQY> &normals,
               const VectorField<DXQY> &tangents, const VectorField<DXQY> &force, const lbBase_t tau);
    void applyNonSlip(const int fieldNo, LbField<DXQY> &f, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const ScalarField &qDists, const VectorField<DXQY> &normals,  const VectorField<DXQY> &tangents, const VectorField<DXQY> &force, const lbBase_t tau);
              
private:
    template <typename MT>
    void matrixFillfa(const int row, MT & m, const int alpha);
    template <typename MT, typename T>
    void matrixFillfaNonSlip(const int row, MT & m, const int alpha, const lbBase_t qDist, const T & normal, const T & tangent, const lbBase_t tau);    
    template <typename MT, typename T>
    void matrixFillSnn(const int row, MT & m, const lbBase_t qDist, T & normal, const lbBase_t un_wall, const lbBase_t Snn);
    template <typename MT, typename T>
    void matrixFillSnt(const int row, MT & m, const lbBase_t qDist, T & normal, T & tangent, const lbBase_t ut_wall,
                       const lbBase_t dun_dt, const lbBase_t Fn,const lbBase_t Ft, const lbBase_t rho0, const lbBase_t tau);
    template <typename MT, typename T>
    void matrixFillStt(const int row, MT & m, T & tangent, const lbBase_t Stt, const lbBase_t Ft, const lbBase_t tau);
    template <typename T>
    inline lbBase_t calcfalpha(const int alpha, const T &x) const;
    template <typename T, typename TV>
    inline lbBase_t calcfalphaNonSlip(const int a, const T &x, const lbBase_t & qDist, const TV &normal, const TV &tangent, const VectorField<DXQY> &force , const lbBase_t & tau) const;
    std::vector<double> rhoWall_;
    lbBase_t boundaryMass_; // Total boundary mass
};

template <typename DXQY>
OneNodeSubGridBndDyn<DXQY>::OneNodeSubGridBndDyn(const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const ScalarField &qDist, const VectorField<DXQY> &normals,  const VectorField<DXQY> &tangents, const ScalarField &rho,  const VectorField<DXQY> &force, const lbBase_t tau)
    : BoundaryExtended<DXQY>(bndNodes, nodes, grid), rhoWall_(bndNodes.size())
{
    boundaryMass_ = 0; // Mass at boundary
    for (int n=0; n < this->size(); ++n) {
        int nodeNo = this->nodeNo(n);
        rhoWall_[n] = rho(0, nodeNo);
        boundaryMass_ += rho(0, nodeNo) - 1;
    }
}

template <typename DXQY>
template <typename MT>
void OneNodeSubGridBndDyn<DXQY>::matrixFillfa(const int row, MT &m, const int a)
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
void OneNodeSubGridBndDyn<DXQY>::matrixFillfaNonSlip(const int row, MT &m, const int a, const lbBase_t qDist, const T & normal, const T & tangent, const lbBase_t tau)
// Assumed a stationary non-slip wall 
{
    int col = 0;
    m(row, col++) = DXQY::w[a];  // rho
    lbBase_t c_dot_t = DXQY::cDot(a, &tangent[0]);
    for (int i=0; i < DXQY::nD; ++i) {  // Pi_ij
        m(row, col) = DXQY::w[a] * (DXQY::c(a, i)*DXQY::c(a, i) - DXQY::c2) * DXQY::c4Inv0_5;
        m(row, col++) -= DXQY::w[a]  * qDist * c_dot_t * normal[i] * tangent[i] * DXQY::c4Inv / tau;
        // m(row, col++) -= DXQY::w[a]  * qDist * (DXQY::c(a, i)* normal[i] )* (DXQY::c4Inv / tau);
        for (int j = i+1; j < DXQY::nD; ++j) {
            m(row, col) = DXQY::w[a] * DXQY::c(a, i)*DXQY::c(a, j) * DXQY::c4Inv;
            m(row, col++) -= DXQY::w[a] * qDist * c_dot_t * (normal[i] * tangent[j] + normal[j] * tangent[i]) * DXQY::c4Inv / tau;            
//            m(row, col++) -= DXQY::w[a]  * qDist * (DXQY::c(a, i)* normal[j] +  DXQY::c(a, j)* normal[i])* (DXQY::c4Inv0_5 / tau);
        }
    }
}

template <typename DXQY>
template <typename MT, typename T>
void OneNodeSubGridBndDyn<DXQY>::matrixFillSnn(const int row, MT &m, const lbBase_t q, T & n, const lbBase_t un, const lbBase_t Snn)
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
void OneNodeSubGridBndDyn<DXQY>::matrixFillSnt(const int row, MT & m, const lbBase_t q, T & n, T & t, const lbBase_t ut,
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
void OneNodeSubGridBndDyn<DXQY>::matrixFillStt(const int row, MT & m, T & t, const lbBase_t Stt,
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
inline lbBase_t OneNodeSubGridBndDyn<DXQY>::calcfalpha(const int a, const T &x) const
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
    return DXQY::w[a]*ret;
}

template <typename DXQY>
template <typename T, typename TV>
inline lbBase_t OneNodeSubGridBndDyn<DXQY>::calcfalphaNonSlip(const int a, const T &x, const lbBase_t & qDist, const TV &normal, const TV &tangent, const VectorField<DXQY> &force, const lbBase_t & tau) const
{
    lbBase_t ret = 0;
    int n = 0;
    
    // Calculate the velocity
    n = 1;
    lbBase_t Pi_nt = 0;
    for (int i = 0; i < DXQY::nD; ++i) {
        Pi_nt += x[n++]*normal[i]*tangent[i];
        for (int j = i+1; j < DXQY::nD; ++j) {
            Pi_nt += (normal[i]*tangent[j] + normal[j]*tangent[i]) * x[n++];
        }
    }    
    lbBase_t rhoVel[DXQY::nD];
    for (int i = 0; i < DXQY::nD; ++i)
    {
        rhoVel[i] = -qDist * Pi_nt * tangent[i] * D2Q9::c2Inv / tau;
    }
    // Calculate distribution
    n=0;
    ret += x[n++]; // Rho
    // \rho u
    ret += DXQY::c2Inv*DXQY::cDot(a, rhoVel);
    // rho uu
    for (int i = 0; i < DXQY::nD; ++i) {
        ret += DXQY::c4Inv0_5*(DXQY::c(a, i)*DXQY::c(a, i)-DXQY::c2)*x[n++];
        for (int j = i+1; j < DXQY::nD; ++j) {
            ret += DXQY::c4Inv0_5*DXQY::c(a, i)*DXQY::c(a, j)*x[n++];
        }
    }
    
    // Add force term
    auto cF = DXQY::cDotAll(force(0, 0));
    ret -= 0.5*cF[a]*D2Q9::c2Inv;
    return D2Q9::w[a]*ret;
}

template <typename DXQY>
void OneNodeSubGridBndDyn<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const ScalarField &qDist, const VectorField<DXQY> &normals,  const VectorField<DXQY> &tangents, const VectorField<DXQY> &force, const lbBase_t tau)
{
    lbBase_t predictedBoundaryMass = 0;
    lbBase_t deltaBulkMass = 0;

    int nCol = 1 + DXQY::nD + (DXQY::nD*(DXQY::nD + 1))/2;
    for (int n = 0; n < this->size(); ++n) {
        int nRow = (1 + this->nBeta(n) + 2*(this->nGamma(n))) + 3;
        MatrixXd m(nRow, nCol);
        VectorXd rhs(nRow);
        // Set the right hand values. Rememver to go through directions in the same order
        // as for the matrix.
        int row = 0;
        int nodeNo = this->nodeNo(n);
        // Known directions
        // -- -- beta_reversed
        for (auto beta: this->beta(n)) {
            int betaRev = this->dirRev(beta);
            matrixFillfa(row, m, betaRev);
            rhs[row++] = f(fieldNo, betaRev, nodeNo);
            // Mass conservation
            int nodeNeigNo = grid.neighbor(beta, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo)) {
                deltaBulkMass += f(fieldNo, beta, nodeNeigNo) - f(fieldNo, betaRev, nodeNo);
            }
        }
        // -- -- gamma and gamma_reversed
        for (auto gamma: this->gamma(n)) {
            int gammaRev = this->dirRev(gamma);
            matrixFillfa(row, m, gamma);
            rhs[row++] = f(fieldNo, gamma, nodeNo);
            matrixFillfa(row, m, gammaRev);
            rhs[row++] = f(fieldNo, gammaRev, nodeNo);
            // Mass conservation
            int nodeNeigNo = grid.neighbor(gammaRev, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo))
                deltaBulkMass += f(fieldNo, gammaRev, nodeNeigNo) - f(fieldNo, gamma, nodeNo);
            nodeNeigNo = grid.neighbor(gamma, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo))
                deltaBulkMass += f(fieldNo, gamma, nodeNeigNo) - f(fieldNo, gammaRev, nodeNo);
        }
        // -- -- zero velocity
        int zeroVelDir = DXQY::nQ-1;
        matrixFillfa(row, m, zeroVelDir);
        rhs[row++] = f(fieldNo, zeroVelDir, nodeNo);
        // Macroscopic boundary conditions
        // -- Setup equations from Snn
        lbBase_t q = qDist(0, nodeNo);
        std::valarray<lbBase_t> nVec = normals(0, nodeNo);
        lbBase_t unWall = 0.0;
        lbBase_t snn = 0.0;
        matrixFillSnn(row, m, q, nVec, unWall, snn);
        rhs[row++] = 0;
        // -- Setup equations from Snt
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
        lbBase_t rho0 = rhoWall_[n];
        matrixFillSnt(row, m, q, nVec, tVec, utWall, dundt, Fn, Ft, rho0, tau);
        rhs[row++] = 0;
        // -- Setup equations frmo Stt
        lbBase_t Stt = 0;
        matrixFillStt(row, m, tVec, Stt, Ft, tau);
        rhs[row++] = 0;

        // Solve system of equations
        JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
        VectorXd x = svd.solve(rhs);

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
        rhoWall_[n] = 0;
        for (int q = 0; q < DXQY::nQ; ++q) {
            f(fieldNo, q, nodeNo) += DXQY::w[q]*addMass;
            rhoWall_[n] += f(fieldNo, q, nodeNo);
        }
    }
}

template <typename DXQY>
void OneNodeSubGridBndDyn<DXQY>::applyNonSlip(const int fieldNo, LbField<DXQY> &f, const Nodes<DXQY> &nodes, 
const Grid<DXQY> &grid, const ScalarField &qDists, const VectorField<DXQY> &normals,  const VectorField<DXQY> &tangents, const VectorField<DXQY> &force, const lbBase_t tau)
{
    lbBase_t predictedBoundaryMass = 0;
    lbBase_t deltaBulkMass = 0;

    int nCol = 1 + (DXQY::nD*(DXQY::nD + 1))/2;
    for (int n = 0; n < this->size(); ++n) {
        int nRow = (1 + this->nBeta(n) + 2*(this->nGamma(n)));
        MatrixXd m(nRow, nCol);
        VectorXd rhs(nRow);
        // Set the right hand values. Rememver to go through directions in the same order
        // as for the matrix.
        int row = 0;
        int nodeNo = this->nodeNo(n);
        lbBase_t qDist = qDists(0, nodeNo);
        std::valarray<lbBase_t> normal = normals(0, nodeNo);
        std::valarray<lbBase_t> tangent = tangents(0, nodeNo);
        // Known directions
        // -- -- beta_reversed
        auto cF = DXQY::cDotAll(force(0, 0));
        for (auto beta: this->beta(n)) {
            int betaRev = this->dirRev(beta);
            
            matrixFillfaNonSlip(row, m, betaRev, qDist, normal, tangent, tau);
            rhs[row++] = f(fieldNo, betaRev, nodeNo) + 0.5*DXQY::w[betaRev]*DXQY::c2Inv*cF[betaRev];
            // Mass conservation
            int nodeNeigNo = grid.neighbor(beta, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo)) {
                deltaBulkMass += f(fieldNo, beta, nodeNeigNo) - f(fieldNo, betaRev, nodeNo);
            }
        }
        // -- -- gamma and gamma_reversed
        for (auto gamma: this->gamma(n)) {
            int gammaRev = this->dirRev(gamma);
            matrixFillfaNonSlip(row, m, gamma, qDist, normal, tangent, tau);
            rhs[row++] = f(fieldNo, gamma, nodeNo)  + 0.5*DXQY::w[gamma]*DXQY::c2Inv*cF[gamma];
            matrixFillfaNonSlip(row, m, gammaRev, qDist, normal, tangent, tau);
            rhs[row++] = f(fieldNo, gammaRev, nodeNo)  + 0.5*DXQY::w[gammaRev]*DXQY::c2Inv*cF[gammaRev];
            // Mass conservation
            int nodeNeigNo = grid.neighbor(gammaRev, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo))
                deltaBulkMass += f(fieldNo, gammaRev, nodeNeigNo) - f(fieldNo, gamma, nodeNo);
            nodeNeigNo = grid.neighbor(gamma, nodeNo);
            if (nodes.isBulkFluid(nodeNeigNo))
                deltaBulkMass += f(fieldNo, gamma, nodeNeigNo) - f(fieldNo, gammaRev, nodeNo);
        }
        // -- -- zero velocity
        int zeroVelDir = DXQY::nQ-1;
        matrixFillfaNonSlip(row, m, zeroVelDir, qDist, normal, tangent, tau);
        rhs[row++] = f(fieldNo, zeroVelDir, nodeNo);

        // Solve system of equations
        JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
        VectorXd x = svd.solve(rhs);

        predictedBoundaryMass += x[0] - 1;
        // Set new distributions
        for (int q = 0; q < DXQY::nQ; ++q) {
            f(fieldNo, q, nodeNo) = calcfalphaNonSlip(q, x, qDist, normal, tangent, force, tau);
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
        rhoWall_[n] = 0;
        for (int q = 0; q < DXQY::nQ; ++q) {
            f(fieldNo, q, nodeNo) += 0*DXQY::w[q]*addMass;
            rhoWall_[n] += f(fieldNo, q, nodeNo);
        }
    }
}



#endif
