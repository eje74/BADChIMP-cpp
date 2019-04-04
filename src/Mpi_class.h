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
#include <algorithm>
#include <mpi.h>

#include "Input.h"
//#include <iomanip>
#include "vector_func.h"


//------------------------------------
//
//------------------------------------
class Mpi {
private:
  std::vector<int> rank_ind_;
  std::vector<int> procs_;                // number of processes in x-, y-, and z-direction
  std::vector<int> send_data_, recv_data_;  // values to send from and receive to this process
  std::vector<int> send_size_, recv_size_;  // size of data-chunks in send_data and recv_data
  std::vector<int> to_proc_, from_proc_;    // list of processes to send send_data to (receive recv_data from)
  int nr_procs_ = 0, max_rank_ = 0;
  int rank_ = -1;
  int num_ghost = 0;
  //int *argc_ = nullptr;
  //char ***argv_ = nullptr;
  int running_ = 0;
  int dim_ = 0;
  //int gn_ = 2;

public:
  //MPI(int *argc, char ***argv, const std::vector<int> &n, const std::vector<int> &procs);
  //  Mpi(int *argc, char ***argv, const std::vector<int> &procs)
  //  : procs_(procs), nr_procs_(prod(procs)), max_rank_(nr_procs_-1),
  //    argc_(argc), argv_(argv), dim_(procs.size()) { };
  Mpi() { };
  Mpi(const std::vector<int> &procs) : procs_(procs), nr_procs_(prod(procs)), max_rank_(nr_procs_-1), dim_(procs.size()) { };
  void start(int *argc, char ***argv, const std::string& node_file, const std::string& rank_file);
  void end() { MPI_Finalize(); }
//  void setup() {
//    Input node_labels("node_labels.mpi"); //node_labels.print();
//    std::vector<int> labels = node_labels["label"];
//    Input rank_map("rank.mpi"); //rank_map.print();
//    std::vector<int> ranks = rank_map["rank"];
//    std::cout << "MPI LABELS: " << labels << std::endl;
//    std::cout << "MPI RANKS: " << ranks << std::endl;
//  }
  void set_rank_ind();
  //void set_N_n_lb_ub();
  void add_ghost_nodes();
  void print() const { std::cout << rank_ << "/" << max_rank_ << std::endl;};
  inline int get_rank() const {return rank_;}
  inline int get_max_rank() const {return max_rank_;}
  inline int get_num_procs(const int i) const {return procs_[i];}
  inline int get_rank_ind(const int i) const {return rank_ind_[i];}
  inline int get_num_ghosts() const {return num_ghost;}
  void print_error_abort(const std::string& err_msg) {
    if (rank_==0) {
      std::cerr << std::endl << "ERROR! " << err_msg << ", aborting... " << std::endl;
    }
    end();
    exit(-1);
  }
  inline int get_rank(const std::vector<int>& pos) {
    return rank_;
  }
};

  //  inline std::vector<int> get_process_position(const int rank) {
//    std::vector<int> rank_ind(dim_);
//    rank_ind_[0] = rank%procs_[0];
//    rank_ind_[1] = ((rank - rank_ind_[0]) / procs_[0]) % procs_[1];
//    rank_ind_[2] = (int)(rank / (procs_[1] * procs_[0]));
//    return rank_ind;
//  }
//  inline std::vector<int> get_local_limits(const int rank, const std::vector<int>& N) {
//    std::vector<int> limits(2*dim_, 0);
//    limits
//    std::vector<int> pros_pos = get_process_position(rank);
//    for (int i = 0; i < dim_; ++i) {
//      int n_tmp = local_.n_[i] / nr_procs_; //mpi.procs_[i];
//      int n_rest = local_.n_[i] % nprocs; // mpi.get_num_procs(i); //mpi.procs_[i];
//    //lb_[i] = mpi.rank_ind_[i] * n_tmp;
//    local_.lb_[i] = rank_ind * n_tmp;
//    local_.lb_[i] += (rank_ind <= n_rest) ? rank_ind : n_rest;
//    if (rank_ind + 1 <= n_rest)
//      ++n_tmp;
//    local_.n_[i] = n_tmp;
//    local_.lb_[i] += num_ghost_;
//    local_.ub_[i] = local_.lb_[i] + local_.n_[i] - num_ghost_;
//          }
//  }





#endif /* SRC_MPI_H_ */
