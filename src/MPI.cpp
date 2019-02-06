/*
 * MPI.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: janlv
 */
#include "MPI.h"


//------------------------------------
//
//------------------------------------
MPI::MPI(int *argc, char ***argv, const std::vector<int> &n, const std::vector<int> &procs)
  : n_(n),
    procs_(procs)
{
  MPI_Init(argc, argv);                      // start up _MPI_
  MPI_Comm_size(MPI_COMM_WORLD, &nr_procs_); // number of processes
  if (procs[0]*procs[1]*procs[2] != nr_procs_) {
    std::cerr << std::endl << "   ERROR: mpi-vector in input-file (" << procs
        << ") different from command-line processes (-n " << nr_procs_ << "), aborting..." << std::endl << std::endl;
    end();
    exit(-1);
  }
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
void MPI::add_ghost_nodes() {
  for (auto& i:n_)
    i+=2;
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
void MPI::set_geometry(std::vector<char>& geo_out, const std::vector<char>& geo_in) {
  std::vector<int> geo_n = N_ - 2;
  std::vector<int> i(get_dim());
  i[2] = n_.size() - 2; // 2D:z=0, 3D:z=1 to skip periodic/ghost rim
  do {
    for(i[1]=1; i[1]<n_[1]-1; ++i[1]) {
      for(i[0]=1; i[0]<n_[0]-1; ++i[0]) {
        int a = get_pos(i, n_);
        int b = get_pos(i+lb_-2, geo_n);
        //std::cout << a << "," << b << ":" << geo_in[b] << std::endl;
        geo_out[a] = geo_in[b];
      }
    }
    ++i[2];
  } while (i[2]<n_[2]-1);
}


//------------------------------------
//
//------------------------------------
void MPI::print() {
  std::cout << rank_ << ": N = " << N_ << ", n = " << n_ << ", lb = " << lb_ << ", ub = " << ub_ << std::endl;
}

