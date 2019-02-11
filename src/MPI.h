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
private:
  std::vector<int> rank_ind_;
  std::vector<int> procs_;
  int nr_procs_ = 0, max_rank_ = 0;
  int rank_ = 0;
  //int gn_ = 2;

public:
  //MPI(int *argc, char ***argv, const std::vector<int> &n, const std::vector<int> &procs);
  MPI(int *argc, char ***argv, const std::vector<int> &procs);
  //~MPI() { MPI_Finalize(); };
  void end() { MPI_Finalize(); };
  void set_rank_ind();
  //void set_N_n_lb_ub();
  void add_ghost_nodes();
  void print();
  inline const int get_rank() const {return rank_;};
  inline const int get_max_rank() const {return max_rank_;};
  inline const int get_num_procs(const int i) const {return procs_[i];};
  inline const int get_rank_ind(const int i) const {return rank_ind_[i];};

};





#endif /* SRC_MPI_H_ */
