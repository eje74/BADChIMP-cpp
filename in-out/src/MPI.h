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

public:
  MPI(int *argc, char ***argv, const std::vector<int> &n, const std::vector<int> &procs);
  //~MPI() { MPI_Finalize(); };
  void end() { MPI_Finalize(); };
  void set_rank_ind();
  void set_N_n_lb_ub();
  void add_ghost_nodes() {for (auto& i:n_) i+=2;};
  void print();
  const int get_global_total_size() const {return N_[0]*N_[1]*N_[2];};
  const int get_local_total_size() const {return n_[0]*n_[1]*n_[2];};
  const std::vector<int>& get_local_system_size() const {return n_;};
  const int get_rank() const {return rank_;};
  const int get_max_rank() const {return max_rank_;};
};





#endif /* SRC_MPI_H_ */
