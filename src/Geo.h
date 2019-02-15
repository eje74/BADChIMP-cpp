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
public:
  struct Limits {
    std::vector<int> n_, lb_, ub_;
    //Limits(int size) : N_(size), n_(size), lb_(size), ub_(size){};
    void assign(std::vector<int>& dim) {
      //size_t size = dim.size();
      n_.assign(dim.begin(), dim.end());
      lb_.assign(dim.size(), 0);
      ub_.assign(dim.begin(), dim.end());
    }
    inline const std::vector<int>& get_lower_bounds() const {return lb_;}
    inline const std::vector<int>& get_upper_bounds() const {return ub_;}
    //inline const std::vector<int>& get_lower_bounds_without_ghosts() const {return(lb_-1*num_ghost_);}
    //inline const std::vector<int>& get_upper_bounds_without_ghosts() const {return(ub_-1*num_ghost_);}
    inline const int get_num_elements() const {return prod(n_);}
    inline const std::vector<int>& get_size() const {return n_;};
    inline const size_t get_dim() const {return n_.size();}
    void print_limits() {
      std::cout << "n = " << n_ << ", lb = " << lb_ << ", ub = " << ub_;
    }
  };
  Limits local_;
  Limits global_;
  int*** labels_ = nullptr;

private:
  int num_ghost_ = 0;
  int rank_ = 0;
  std::vector<char> nodes_;
  double dx_ = 0.0;
  //size_t dim_ = 0;

public:
  Geo(const std::string& filename, Mpi& mpi) : rank_(mpi.get_rank())
  {
    Input geo(filename); //geo.print();
    std::vector<int> dim = geo["dim"];
    local_.assign(dim);
    global_.assign(dim);
    dx_ = geo["res"];

    //set_limits(mpi);
    //add_ghost_nodes(1);
    //set_nodes(geo);
  }
  void print_limits() {
    std::cout << "("<< rank_ << ") ghosts: " << num_ghost_ << ", global: ";
    global_.print_limits();
    std::cout << ", local: ";
    local_.print_limits();
    std::cout << std::endl;
  };
  //void set_labels(int*** labels) {labels_ = labels;};
  inline const int get_pos(const std::vector<int>& ind, const std::vector<int>& stride) {
    if (ind.size()>2)
      return ind[0] + ind[1]*stride[0] + ind[2]*stride[0]*stride[1];
    else
      return ind[0] + ind[1]*stride[0];
  }
  //inline const int global_to_local_pos(std::vector<int>& ind) { return(get_pos(ind-lb_+1, n_)); }
  inline const int local_to_global_pos(const std::vector<int>& ind) { return(get_pos(ind+local_.lb_-1, global_.n_)); }
  inline const int get_local_pos(const std::vector<int>& ind) { return(get_pos(ind, local_.n_)); }
  inline const int get_geofile_pos(const std::vector<int>& ind) { return(get_pos(ind+local_.lb_-2, global_.n_-2)); }
  //inline const std::vector<int>& get_global_size_without_ghosts() const { return(global_.n_-2*num_ghost_); };
  inline const int get_num_ghosts(){ return(num_ghost_); };

private:
  void set_limits(Mpi& mpi);
  void set_nodes(const std::vector<char>& geo);
  void add_ghost_nodes(int num);
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
