/*
 * Geometry.h
 *
 *  Created on: Feb 7, 2019
 *      Author: janlv
 */

#ifndef SRC_GEO_H_
#define SRC_GEO_H_
#include <vector>
#include "Input.h"
#include "Mpi.h"

class Geo {
private:
  std::vector<int> N_;
  std::vector<int> n_;
  std::vector<int> lb_, ub_;
  int rank_ = 0;
  std::vector<char> nodes;
  double dx = 0.0;
  int dim = 0;

public:
  Geo(const std::string& filename, Mpi& mpi) : rank_(mpi.get_rank())
  {
    Input geo(filename); //geo.print();
    dim = geo["dim"].ncols();
    N_.resize(dim);
    n_.resize(dim);
    ub_.resize(dim);
    lb_.assign(dim, 0);
    N_ = n_ = ub_ = geo["dim"];
    dx = geo["res"];

    //set_limits(mpi);
    //add_ghost_nodes();
    //set_nodes(geo);
  }
  void print_limits() {
    std::cout << rank_ << ": N = " << N_ << ", n = " << n_ << ", lb = " << lb_ << ", ub = " << ub_ << std::endl;
  }
  inline const std::vector<int>& get_lower_bounds() const {return lb_;}
  inline const std::vector<int>& get_upper_bounds() const {return ub_;}
  inline const int get_global_num_elements() const {return prod(N_);}
  inline const int get_local_num_elements() const {return prod(n_);};
  inline const std::vector<int>& get_local_size() const {return n_;};
  inline const std::vector<int>& get_global_size() const {return N_;};
  inline const size_t get_dim() const {return n_.size();}
  inline const int get_pos(const std::vector<int>& ind, const std::vector<int>& stride) {
    if (ind.size()>2)
      return ind[0] + ind[1]*stride[0] + ind[2]*stride[0]*stride[1];
    else
      return ind[0] + ind[1]*stride[0];
  }
  //inline const int global_to_local_pos(std::vector<int>& ind) { return(get_pos(ind-lb_+1, n_)); }
  inline const int local_to_global_pos(const std::vector<int>& ind) { return(get_pos(ind+lb_-1, N_)); }
  inline const int get_local_pos(const std::vector<int>& ind) { return(get_pos(ind, n_)); }
  inline const int get_geofile_pos(const std::vector<int>& ind) { return(get_pos(ind+lb_-2, N_-2)); }

private:
  void set_limits(Mpi& mpi);
  void set_nodes(const std::vector<char>& geo);
  void add_ghost_nodes();
};




//------------------------------------
//
//------------------------------------
//void Mpi::set_geometry(std::vector<char>& geo_out, const std::vector<char>& geo_in) {
//  std::vector<int> geo_n = N_ - 2;
//  std::vector<int> i(get_dim());
//  i[2] = n_.size() - 2; // 2D:z=0, 3D:z=1 to skip periodic/ghost rim
//  do {
//    for(i[1]=1; i[1]<n_[1]-1; ++i[1]) {
//      for(i[0]=1; i[0]<n_[0]-1; ++i[0]) {
//        int a = get_pos(i, n_);
//        int b = get_pos(i+lb_-2, geo_n);
//        //std::cout << a << "," << b << ":" << geo_in[b] << std::endl;
//        geo_out[a] = geo_in[b];
//      }
//    }
//    ++i[2];
//  } while (i[2]<n_[2]-1);
//}




#endif /* SRC_GEO_H_ */
