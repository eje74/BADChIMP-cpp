#ifndef LBBOUNDARYINTERPOLATION_H
#define LBBOUNDARYINTERPOLATION_H

#include "../LBSOLVER.h"
#include "../IO.h"
// #include "LBrans.h"
// #include <chrono>
#include <numeric>


struct InterpolationElement
{
    int nodeNo;
    double gamma, gamma2;
    double surfaceWeight;
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
      auto ret = readBoundaryNodeEntry(line, grid);
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

      boundaryInterpolation[bndNo] = ret;
      //std::cout << bndNo << std::endl;    

      //InterpolationElement &el = boundaryInterpolation[bndNo];        

      //el.gamma = 
      /* auto posBnd = grid.pos(ret.nodeNo);
      double xb, yb;
      xb = posBnd[0] + ret.gamma2*ret.normal[0];
      yb = posBnd[1] + ret.gamma2*ret.normal[1];
      double x, y;
      x = 0.0;
      y = 0.0;
      for (int i=0; i < 4; ++i) {
        auto pos = grid.pos(ret.pnts[i]);
        x += pos[0]*ret.wb[i];
        y += pos[1]*ret.wb[i];
      }
      if ( std::abs(xb - x) > 1e-12  || std::abs(yb - y) > 1e-12 ) 
        std::cout << xb - x << "  " << yb - y << std::endl;
        */
      /* std::cout << ret.nodeNo << " "; 
      std::cout << ret.gamma << " ";
      for (auto &p: ret.pnt)
        std::cout << p << " ";
      std::cout << ret.gamma2 << " ";
      for (auto &w: ret.wa)
        std::cout << w << " ";
      for (auto &w: ret.wb)
        std::cout << w << " ";
      for (auto &w: ret.wc)
        std::cout << w << " ";
      for (auto &n: ret.normal)
        std::cout << n << " ";
      std::cout << std::endl; */
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

  return ret;
}



#endif