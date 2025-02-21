#ifndef LBBOUNDARYINTERPOLATION_H
#define LBBOUNDARYINTERPOLATION_H

#include "../LBSOLVER.h"
#include "../IO.h"
#include "LBrans.h"
// #include <chrono>
#include <numeric>


struct InterpolationElement
{
    int nodeNo;
    double gamma, gamma2;
    double surfaceWeight;
    double E, yp, ycut_off;
    double uWallMin, uWallMax;
    double lowC0, lowC1;
    std::vector<int> pnts;
    std::vector<double> wa, wb, wc;
    std::vector<double> normal;  
};
//-------------------------------------------------------------------------------------     boundaryMap [declaration]
template<typename DXQY>
std::vector<int> boundaryMap(const Boundary<DXQY> &boundary, Grid<DXQY> &grid);
//-------------------------------------------------------------------------------------     readBoundaryNodeEntry [declaration]
template<typename DXQY>
InterpolationElement readBoundaryNodeEntry(std::string &line, Grid<DXQY> &grid);

template<typename DXQY>
void readBoundaryNodeFile(
    const std::string &fileName,
    std::vector<InterpolationElement> &boundaryInterpolation,
    const Boundary<DXQY> &boundary,
    lbBase_t viscocity0,
    lbBase_t kappa,
    Nodes<DXQY> &nodes,
    Grid<DXQY> &grid,
    LBvtk<DXQY> &vtklb)
{
  std::ifstream ifs(fileName);
  std::string line;

  if (!ifs.is_open())
  {
    std::cerr << "Unable to open file " << fileName << std::endl;
    MPI_Finalize();
    exit(1);    
  } else {
    auto bndLabel = boundaryMap(boundary, grid);

    while (std::getline(ifs, line))
    {
      InterpolationElement ret = readBoundaryNodeEntry(line, grid);
      auto cn = DXQY::cDotAll(ret.normal);
      int bndNo = bndLabel[ret.nodeNo];
      double boundaryWeight = 0.0;
      for (auto q: boundary.unknown(bndNo)) {
        if (cn[q] < 0) {
          std::cout << cn[q] << " error" << std::endl;
          exit(1);
        }

        boundaryWeight += 6.0*DXQY::w[q]*cn[q];
      }
      ret.surfaceWeight = boundaryWeight;

      // Add local boundary constants
      // ret.uWallMin = 10.0*viscosity0_*std::log(ret.E*10.0)/(kappa_*ret.yp);
      // ret.uWallMax = 300.0*viscosity0_*std::log(ret.E*300.0)/(kappa_*ret.yp);
      // ret.lowC0 = 1.0/kappa_;
      // ret.lowC1 = ret.E*yp_/viscosity0_; 


      boundaryInterpolation[bndNo] = ret;

    }
  }

  std::cout << "READ BOUNDARY NODE FILE " << std::endl;
}

//-------------------------------------------------------------------------------------     boundaryMap [definition]
template<typename DXQY>
std::vector<int> boundaryMap(const Boundary<DXQY> &boundary, Grid<DXQY> &grid)
{
  std::vector<int> ret(grid.size(), 0);
  for (int bndNo=0; bndNo < boundary.size(); ++bndNo)
  {
    const int nodeNo = boundary.nodeNo(bndNo);
    ret[nodeNo] = bndNo;
  }
  return ret;
}

template<typename DXQY>
int readNodeNumber(std::stringstream &iss, Grid<DXQY> &grid)
{
  std::string word;
  int i, j;

  iss >> word;
  i = std::stoi(word);
  iss >> word;
  j = std::stoi(word);

  int tmp[] = {i, j};
  std::vector<int> ij(tmp, tmp + DXQY::nD);

  return grid.nodeNo(ij); 
}

double readDouble(std::stringstream &iss)
{
  std::string word;
  iss >> word;
  return std::stod(word);
}

template<typename DXQY>
std::vector<int> readVectorNodeNo(int numNodes, std::stringstream &iss, Grid<DXQY> &grid)
{
  std::vector<int> ret(numNodes, 0);
  for (int i=0; i < numNodes; ++i)
  {
    ret[i] = readNodeNumber(iss, grid);
  }
  return ret;
}

std::vector<double> readVectorDouble(int numVals, std::stringstream &iss)
{
  std::vector<double> ret(numVals, 0.0);
  for (auto &v: ret)
    v = readDouble(iss);

  return ret; 
}

template<typename DXQY>
InterpolationElement readBoundaryNodeEntry(std::string &line, Grid<DXQY> &grid)
{
  InterpolationElement ret;
  std::stringstream iss(line);

  ret.nodeNo = readNodeNumber(iss, grid);
  ret.gamma  = readDouble(iss);
  ret.pnts = readVectorNodeNo(4, iss, grid);
  ret.gamma2 = readDouble(iss);
  ret.wa = readVectorDouble(ret.pnts.size(), iss);
  ret.wb = readVectorDouble(ret.pnts.size(), iss);
  ret.wc = readVectorDouble(ret.pnts.size(), iss);
  ret.normal = readVectorDouble(DXQY::nD, iss);
  ret.E = readDouble(iss);
  ret.yp = readDouble(iss);
  ret.ycut_off = readDouble(iss);

  return ret;
}



#endif