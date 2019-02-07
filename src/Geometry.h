/*
 * Geometry.h
 *
 *  Created on: Feb 7, 2019
 *      Author: janlv
 */

#ifndef SRC_GEOMETRY_H_
#define SRC_GEOMETRY_H_
#include <vector>
#include "Input.h"
#include "MPI.h"

class Geo {
public:
  Geo(std::vector<char>& geo, MPI& mpi)
  {
    set_nodes(geo, mpi);
  }

private:
  std::vector<char> nodes;
  void set_nodes(const std::vector<char>& geo, MPI& mpi);
  const int get_pos(const std::vector<int>& ind, const std::vector<int>& stride) {
    return ind[0] + ind[1]*stride[0] + ind[2]*stride[0]*stride[1]; }
};

//------------------------------------
//
//------------------------------------
void Geo::set_nodes(std::vector<char>& geo, MPI& mpi) {
  std::vector<int> n = mpi.get_local_system_size();
  std::vector<int> i(n.size());
  i[2] = n.size() - 2; // 2D:z=0, 3D:z=1 to skip periodic/ghost rim
  do {
    for(i[1]=1; i[1]<n[1]-1; ++i[1]) {
      for(i[0]=1; i[0]<n[0]-1; ++i[0]) {
        int a = mpi.get_local_pos(i);
        int b = mpi.get_geofile_pos(i);
        //std::cout << a << "," << b << ":" << geo_in[b] << std::endl;
        nodes[a] = geo[b];
      }
    }
    ++i[2];
  } while (i[2]<n[2]-1);
}





#endif /* SRC_GEOMETRY_H_ */
