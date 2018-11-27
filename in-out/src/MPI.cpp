/*
 * MPI.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: janlv
 */
#include "MPI.h"

//------------------------------------
// << operator for std::vector
//------------------------------------
template <typename T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &v) {
  os << "(";
  std::string sep = "";
  for (const auto& i : v) {
    os << sep << i;
    sep = ",";
  }
  os << ")";
  return os;
}

//------------------------------------
//
//------------------------------------
MPI::MPI(int *argc, char ***argv, const std::vector<int> &n, const std::vector<int> &procs)
  : n_(n),
    procs_(procs)
{
  MPI_Init(argc, argv);                      // start up _MPI_
  MPI_Comm_size(MPI_COMM_WORLD, &nr_procs_); // number of processes
  max_rank_ = nr_procs_ - 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);     // process rank

  set_rank_ind();
  set_N_n_lb_ub();
  add_ghost_nodes();
  //std::cout << "MPI: " << rank_ << "/" << nr_procs_ << std::endl;
};


//------------------------------------
//
//------------------------------------
void MPI::set_rank_ind() {
  rank_ind_.resize(3);
  rank_ind_[0] = rank_%procs_[0];
  rank_ind_[1] = ((rank_ - rank_ind_[0]) / procs_[0]) % procs_[1];
  rank_ind_[2] = (int)(rank_ / (procs_[1] * procs_[0]));
}

//------------------------------------
// Load distribution
//------------------------------------
void MPI::set_N_n_lb_ub() {
  N_.resize(n_.size());
  lb_.resize(n_.size());
  ub_.resize(n_.size());
  for (size_t i = 0; i < N_.size(); ++i) {
    N_[i] = n_[i] + 2;
    int n_tmp = n_[i] / procs_[i];
    int n_rest = n_[i] % procs_[i];
    lb_[i] = rank_ind_[i] * n_tmp;
    lb_[i] += (rank_ind_[i] <= n_rest) ? rank_ind_[i] : n_rest;
    if (rank_ind_[i] + 1 <= n_rest)
      ++n_tmp;
    n_[i] = n_tmp;
    lb_[i] += 1;
    ub_[i] = lb_[i] + n_[i] - 1;
  }
}

//------------------------------------
//
//------------------------------------
void MPI::print() {
  std::cout << rank_ << ": N = " << N_ << ", n = " << n_ << ", lb = " << lb_ << ", ub = " << ub_ << std::endl;
}

