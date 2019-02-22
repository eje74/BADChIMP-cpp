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
class Mpi {
private:
  std::vector<int> rank_ind_;
  std::vector<int> procs_;
  int nr_procs_ = 0, max_rank_ = 0;
  int rank_ = 0;
  int num_ghost = 0;
  //int gn_ = 2;

public:
  //MPI(int *argc, char ***argv, const std::vector<int> &n, const std::vector<int> &procs);
  Mpi(int *argc, char ***argv, const std::vector<int> &procs);
  //~Mpi() { MPI_Finalize(); };
  void end() { MPI_Finalize(); }
  void set_rank_ind();
  //void set_N_n_lb_ub();
  void add_ghost_nodes();
  void print();
  inline int get_rank() const {return rank_;}
  inline int get_max_rank() const {return max_rank_;}
  inline int get_num_procs(const int i) const {return procs_[i];}
  inline int get_rank_ind(const int i) const {return rank_ind_[i];}
  inline int get_num_ghosts() const {return num_ghost;}

};





#endif /* SRC_MPI_H_ */
