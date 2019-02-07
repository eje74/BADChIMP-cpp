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
  inline const int get_global_size() const {return prod(N_);}
  inline const int get_local_size() const {return prod(n_);};
  inline const std::vector<int>& get_local_system_size() const {return n_;};
  inline const int get_rank() const {return rank_;};
  inline const int get_max_rank() const {return max_rank_;};
  inline const size_t get_dim() const {return n_.size();}
  void set_geometry(std::vector<char>& geo_out, const std::vector<char>& geo_in);
  inline const int get_pos(const std::vector<int>& ind, const std::vector<int>& stride) {
    return ind[0] + ind[1]*stride[0] + ind[2]*stride[0]*stride[1]; }
  inline const int global_to_local_pos(const std::vector<int>& ind) { return(get_pos(ind-lb_+1, n_)); }
  inline const int local_to_global_pos(const std::vector<int>& ind) { return(get_pos(ind+lb_-1, N_)); }
  inline const int get_local_pos(const std::vector<int>& ind) { return(get_pos(ind, n_)); }
  inline const int get_geofile_pos(const std::vector<int>& ind) { return(get_pos(ind+lb_-2, N_-2)); }

};





#endif /* SRC_MPI_H_ */
