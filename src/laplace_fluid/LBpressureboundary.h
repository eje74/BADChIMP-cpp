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
  void apply(const int fieldNum, LbField<DXQY> &f, const VectorField<DXQY> &bndNorm, const T &rho_bnd, const VectorField<DXQY> &vel, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
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
