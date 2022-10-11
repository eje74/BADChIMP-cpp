#ifndef LBPRESSUREBOUNDARY_H
#define LBPRESSUREBOUNDARY_H

#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBvtk.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBboundary.h"

#include <array>
#include <set>


//=====================================================================================
//
//                  R E G U L A R I Z E D   D I S T R I B U T I O N
//
//=====================================================================================
template<int DIM>
inline lbBase_t lowerDiagCcCont(const std::valarray<lbBase_t> &m, const std::valarray<lbBase_t> & c)
{
    std::cout << "Error in template specialization" << std::endl;
    exit(1);
    return 0;
}

template<>
inline lbBase_t lowerDiagCcCont<2>(const std::valarray<lbBase_t> &m, const std::valarray<lbBase_t> & c)
{
    return c[0]*c[0]*m[0] + 2*c[0]*c[1]*m[1] + c[1]*c[1]*m[2];
}

template<>
inline lbBase_t lowerDiagCcCont<3>(const std::valarray<lbBase_t> &m, const std::valarray<lbBase_t> & c)
{
    return c[0]*c[0]*m[0] + 2*c[1]*c[0]*m[1] + c[1]*c[1]*m[2] + 2*c[2]*c[0]*m[3] + 2*c[2]*c[1]*m[4] + c[2]*c[2]*m[5];
}

template<int DIM>
inline lbBase_t lowerDiagTrace(const std::valarray<lbBase_t> &m)
{
    std::cout << "Error in template specialization" << std::endl;
    exit(1);
    return 0;
}

template<>
inline lbBase_t lowerDiagTrace<2>(const std::valarray<lbBase_t> &m)
{
    return m[0] + m[2];   
}

template<>
inline lbBase_t lowerDiagTrace<3>(const std::valarray<lbBase_t> &m)
{
    return m[0] + m[2] + m[5];   
}

template<typename DXQY>
std::valarray<lbBase_t> fRegularized(const std::valarray<lbBase_t> &f, const lbBase_t dRhoNode)
{
    const lbBase_t M = DXQY::qSum(f) + dRhoNode;
    const auto Mi = DXQY::qSumC(f);
    const auto Mij = DXQY::qSumCCLowTri(f);

    std::valarray<lbBase_t> fReg(DXQY::nQ);
    const lbBase_t c2traceM = DXQY::c2*lowerDiagTrace<DXQY::nD>(Mij);
    for (int q=0; q<DXQY::nQ; ++q) {
        const lbBase_t Qdelta = DXQY::cNorm[q]*DXQY::cNorm[q] - DXQY::nD*DXQY::c2;
        const lbBase_t cM = DXQY::cDotRef(q, Mi)*DXQY::c2Inv;        
        const lbBase_t QM = DXQY::c4Inv0_5*(lowerDiagCcCont<DXQY::nD>(Mij, DXQY::cValarray(q)) - c2traceM - DXQY::c2*Qdelta*M);
        fReg[q] = DXQY::w[q]*(M + cM + QM);
    }

    return fReg;
}

template<typename DXQY>
std::valarray<lbBase_t> fRegularizedNormVel(const std::valarray<lbBase_t> &f, const std::valarray<lbBase_t> &normVec, const std::valarray<lbBase_t> &F)
{
    const lbBase_t M = DXQY::qSum(f);
    const auto MiOrg = DXQY::qSumC(f);
    const auto Mij = DXQY::qSumCCLowTri(f);

    const std::valarray<lbBase_t> Mi = DXQY::dot(MiOrg,normVec)*normVec + 0.5*(DXQY::dot(F, normVec)*normVec - F);

    std::valarray<lbBase_t> fReg(DXQY::nQ);
    const lbBase_t c2traceM = DXQY::c2*lowerDiagTrace<DXQY::nD>(Mij);
    for (int q=0; q<DXQY::nQ; ++q) {
        const lbBase_t Qdelta = DXQY::cNorm[q]*DXQY::cNorm[q] - DXQY::nD*DXQY::c2;
        const lbBase_t cM = DXQY::cDotRef(q, Mi)*DXQY::c2Inv;        
        const lbBase_t QM = DXQY::c4Inv0_5*(lowerDiagCcCont<DXQY::nD>(Mij, DXQY::cValarray(q)) - c2traceM - DXQY::c2*Qdelta*M);
        fReg[q] = DXQY::w[q]*(M + cM + QM);
    }

    return fReg;
}




template <typename DXQY>
int addPressureBoundary(const std::string &attr, LBvtk<DXQY> &vtklb, Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
/*
 *  Adds pressure boundaries to the nodes-object.
 *  Returns the number of pressure boundaries.
 *
 *  Reads the "attr" attribute from the vtklb-file.
 *  - Assums integer entries:
 *       0 : do nothing
 *      -1 : set to solid boundary
 *     > 0 : set to fluid boundary if it is a fluid node.
 *           set the nodes tag equal to the attribute.
 *
 */
{
  int maxBoundaryIndicatorLocal = 0;
  vtklb.toAttribute(attr);
  for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++)
  {
    const auto pInd = vtklb.template getScalarAttribute<int>();
    maxBoundaryIndicatorLocal = std::max(maxBoundaryIndicatorLocal, pInd);
    if (pInd == -1)
    { // Solid boundary
      nodes.addNodeType(1, nodeNo);
      nodes.addNodeTag(pInd, nodeNo);
    }
    else
    {
      nodes.addNodeTag(0, nodeNo);
    }
  }

  vtklb.toAttribute(attr);
  for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++)
  {
    const auto pInd = vtklb.template getScalarAttribute<int>();

    if (pInd > 0)
    { // Fluid boundary
      if (nodes.isFluid(nodeNo))
      {
        for (int q = 0; q < DXQY::nQ; q++)
        {
          if (nodes.getTag(grid.neighbor(q, nodeNo)) == -1)
          {
            nodes.addNodeType(2, nodeNo);
            nodes.addNodeTag(pInd, nodeNo);
            break;
          }
        }
      }
    }
  }

  // Maximum boundary indicator
  int maxBoundaryIndicator;
  MPI_Allreduce(&maxBoundaryIndicatorLocal, &maxBoundaryIndicator, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return maxBoundaryIndicator;
}

template <typename DXQY>
class fluidPressureBoundary
{
public:
  fluidPressureBoundary() {}

  fluidPressureBoundary(const std::vector<int> boundaryNodes, Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
  {

    boundary_ = Boundary<DXQY>(boundaryNodes, nodes, grid);
  }

  template <typename T>
  void apply(const int fieldNum, LbField<DXQY> &f, const VectorField<DXQY> &bndNorm, const VectorField<DXQY> &force, const T &rho_bnd, const VectorField<DXQY> &vel, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
  {

    for (const auto &bn : boundary_())
    {
      const int nodeNo = bn.nodeNo();
      const int tag = nodes.getTag(nodeNo);
      if (tag > 0)
      {
        const std::valarray<lbBase_t> v = vel(fieldNum, nodeNo);
        // const std::valarray <lbBase_t> v = DXQY::dot(vel(fieldNum, nodeNo), bndNorm(0, tag-1))*bndNorm(0, tag-1);
        lbBase_t vSq = DXQY::dot(v, v);
        for (const auto &q : bn.unknown())
        {
          int qrev = DXQY::reverseDirection(q);
          int neighNo = grid.neighbor(qrev, nodeNo);
          if (nodes.getTag(neighNo) == -1)
          {
            lbBase_t cDotv = DXQY::cDotRef(qrev, v);
            f(fieldNum, q, nodeNo) = -f(fieldNum, qrev, neighNo) + 2 * DXQY::w[q] * rho_bnd[tag - 1] * (1 + 0.5 * DXQY::c4Inv * (cDotv * cDotv - DXQY::c2 * vSq));
          }
          else
          {
            f(fieldNum, q, nodeNo) = f(fieldNum, qrev, neighNo);
          }
        }
        // f.set(fieldNum, nodeNo) = fRegularized<DXQY>(f(fieldNum, nodeNo), 0.0);
        const std::valarray<lbBase_t> surfNorm = bndNorm(0, tag-1);
        const std::valarray<lbBase_t> forceNode = force(0, nodeNo);
        f.set(fieldNum, nodeNo) = fRegularizedNormVel<DXQY>(f(fieldNum, nodeNo), surfNorm, forceNode);  
      }
      else if (tag == 0)
      {
        for (const auto &q : bn.unknown())
        {
          int qrev = DXQY::reverseDirection(q);
          f(fieldNum, q, nodeNo) = f(fieldNum, qrev, grid.neighbor(qrev, nodeNo));
        }
      }
      else
      {
        std::cout << "ERROR: fluidPRessureBoundary.apply() does not recognize tag =" << tag << std::endl;
      }
    }
  }

template <typename T>
  void applyTest(const int fieldNum, LbField<DXQY> &f, const VectorField<DXQY> &bndNorm, const T &rho_bnd, const VectorField<DXQY> &vel, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, ScalarField & testNumSolid, VectorField<DXQY> & testDirSolid)
  {

    for (const auto &bn : boundary_())
    {
      const int nodeNo = bn.nodeNo();
      const int tag = nodes.getTag(nodeNo);
      if (tag > 0)
      {
        const std::valarray<lbBase_t> v = vel(fieldNum, nodeNo);
        // const std::valarray <lbBase_t> v = DXQY::dot(vel(fieldNum, nodeNo), bndNorm(0, tag-1))*bndNorm(0, tag-1);
        lbBase_t vSq = DXQY::dot(v, v);
        std::valarray<lbBase_t> cSum0(0.0, DXQY::nQ);
        std::valarray<lbBase_t> cSum1(0.0, DXQY::nQ);
        for (const auto &q : bn.unknown())
        {
          int qrev = DXQY::reverseDirection(q);
          int neighNo = grid.neighbor(qrev, nodeNo);
          cSum0 += DXQY::cValarray(q);
          testNumSolid(0, nodeNo) += 1;
          if (nodes.getTag(neighNo) == -1)
          {
            lbBase_t cDotv = DXQY::cDotRef(qrev, v);
            f(fieldNum, q, nodeNo) = -f(fieldNum, qrev, neighNo) + 2 * DXQY::w[q] * rho_bnd[tag - 1] * (1 + 0.5 * DXQY::c4Inv * (cDotv * cDotv - DXQY::c2 * vSq));
            cSum1 += DXQY::cValarray(q);
            testNumSolid(1, nodeNo) += 1;
          }
          else
          {
            f(fieldNum, q, nodeNo) = f(fieldNum, qrev, neighNo);
          }
        }
        testDirSolid.set(0, nodeNo) = cSum0;
        testDirSolid.set(1, nodeNo) = cSum1;
        f.set(fieldNum, nodeNo) = fRegularized<DXQY>(f(fieldNum, nodeNo), 0.0);


      }
      else if (tag == 0)
      {
        for (const auto &q : bn.unknown())
        {
          int qrev = DXQY::reverseDirection(q);
          f(fieldNum, q, nodeNo) = f(fieldNum, qrev, grid.neighbor(qrev, nodeNo));
        }
      }
      else
      {
        std::cout << "ERROR: fluidPRessureBoundary.apply() does not recognize tag =" << tag << std::endl;
      }
    }
  }


private:
  Boundary<DXQY> boundary_;
};



#endif
