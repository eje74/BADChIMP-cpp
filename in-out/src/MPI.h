/*
 * mpi.h
 *
 *  Created on: Nov 19, 2018
 *      Author: janlv
 */

#ifndef SRC_MPI_H_
#define SRC_MPI_H_
#include <vector>
#include <iostream>
#include <mpi.h>
//#include <iomanip>
#include "vector_func.h"


//------------------------------------
//
//------------------------------------
class MPI {
public:
  std::vector<int> N_;
  std::vector<int> n_;
  std::vector<int> lb_, ub_;
private:
  std::vector<int> rank_ind_;
  std::vector<int> procs_;
  int nr_procs_ = 0, max_rank_ = 0;
  int rank_ = 0;
  //int gn_ = 2;

public:
  MPI(int *argc, char ***argv, const std::vector<int> &n, const std::vector<int> &procs);
  //~MPI() { MPI_Finalize(); };
  void end() { MPI_Finalize(); };
  void set_rank_ind();
  void set_N_n_lb_ub();
  void add_ghost_nodes();
  void print();
  const int get_global_size() const {return N_[0]*N_[1]*N_[2];};
  const int get_local_size() const {return prod(n_);}
  const std::vector<int>& get_local_system_size() const {return n_;};
  const std::vector<int>& get_global_system_size() const {return N_;};
  const int get_rank() const {return rank_;};
  const int get_max_rank() const {return max_rank_;};
  const size_t get_dim() const {return n_.size();}
  void set_geometry(std::vector<char>& geo_out, const std::vector<char>& geo_in);
  //const int global_to_local_pos(const std::vector<int>& N) { return(get_pos(N-lb_+1, n_)); }
  //const int local_to_global_pos(const std::vector<int>& n) { return(get_pos(n+lb_-1, N_)); }
  const int get_pos(const std::vector<int>& ind, const std::vector<int>& stride) {
    return ind[0] + ind[1]*stride[0] + ind[2]*stride[0]*stride[1]; }
  const int get_local_pos(const std::vector<int>& ind) { return get_pos(ind, n_); }
  const int get_global_pos(const std::vector<int>& ind) { return get_pos(ind + lb_-2, N_-2); }

};





#endif /* SRC_MPI_H_ */
