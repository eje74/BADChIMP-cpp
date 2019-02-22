/*
 * output.h
 *
 *  Created on: 27. nov. 2017
 *      Author: janlv
 */

#ifndef SRC_OUTPUT_H_
#define SRC_OUTPUT_H_
#include <unordered_map>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <bitset>
//#include <sys/types.h>
#include <sys/stat.h>
#if defined (_WIN32)
#include <direct.h>
#else
#include <unistd.h>
#endif
#include "Mpi.h"
#include "Geo.h"
//#include "global.h"
#include "defines.h"
#include "vector_func.h"

//struct Node;


//--------------------------------
//
//--------------------------------
class Variable {
public:
  std::string name;
  void *data_pointer = nullptr;
  int datasize = 0;
  int dim = 0;
  int data_stride = 0; // step between data-nodes
  int nbytes = 0;
  int offset = 0;
  std::string datatype;
  std::bitset<NUM_MODE> write_node;
  std::string data_array;

  // constructor
  Variable(const std::string &_name, void *_data_pointer, int _datasize, int _dim, int _data_stride, const std::vector<int> &_n)
  : name(_name),
    data_pointer(_data_pointer),
    datasize(_datasize),
    dim(_dim),
    data_stride(_data_stride)
  {
    datatype = "Float" + std::to_string(sizeof(OUTPUT_DTYPE)*8);
    set_write_nodes({FLUID, WALL, SOLID});
    set_nbytes(_n);
    set_data_array();
  }

  void set_nbytes(const std::vector<int> &n);
  void set_offset(const Variable &prev_var);
  void set_write_nodes(const std::vector<int> &_node_mode_to_write);
  void set_data_array();
  template <typename T>
  T get_data(int nn, int dim) const {
    int element = data_stride*nn + dim;
    //int element = nn + dim;
    if (datasize==1) {
      char *char_pointer = static_cast<char*>(data_pointer);
      return T(char_pointer[element]);
    } else if (datasize==8) {
      double *double_pointer = static_cast<double*>(data_pointer);
      return T(double_pointer[element]);
    } else {
      printf("\nERROR in Variable::get_data: unrecognized datasize: %zu\n", static_cast<size_t>(datasize));
      MPI_Finalize();
      return -1;
    }
  }

};

//--------------------------------------------
//
//--------------------------------------------
class File {
protected:
  std::string name_;
  std::ofstream file_;
  std::string filename_;
  std::vector<std::string> folders_;
  std::string path_;
  std::string extension_;
  //int rank_ = 0, max_rank_ = 0;
  int num_ghosts_ = 0;
  //Mpi mpi_;
  //std::string path_filename_;

public:
  int nwrite_ = 0;
  File(const std::string &name, const std::vector<std::string> &folders , const std::string &extension, const Mpi &mpi)
: name_(name),
  folders_(folders),
  extension_(extension),
  num_ghosts_(mpi.get_num_ghosts())
  //rank_(mpi.get_rank()),
  //max_rank_(mpi.get_max_rank()),
{
    for (const auto& f : folders_) {
      path_ += f;
      make_dir(path_);
    }
}
  const std::string& get_filename() const {return filename_;}

private:
  void make_dir(std::string &dir);
};

//--------------------------------------------
//
//--------------------------------------------
class VTI_file : public File {
private:
  std::vector<Variable> variables_;
  std::vector<int> n;
  std::string type, extent;
  std::string cell_data_string_;
  //Node *node_ = nullptr;

public:
  //  VTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, Node *node)
  //  : File(_name, {_path, "vti/"}, ".vti", _mpi), n(_geo.get_local_size()), node_(geo.nodes) { set_extent(_geo); }
  VTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, const Mpi& _mpi)
  : File(_name, {_path, "vti/"}, ".vti", _mpi), n(_geo.get_n()) { set_extent(_geo); }

  void set_extent(const Geo &geo);
  //void write(int*** labels);
  //void write(int*** labels,const Variable& var);
  void write_data(int*** labels);
  void write_header();
  void write_footer_and_close();
  void set_filename_and_open(const int rank);
  //std::string get_filename();
  const std::string get_piece_extent_string() const;
  const std::vector<int>& get_system_size() const {return n;}
  const std::string get_cell_data_string() const {return cell_data_string_;}
  void set_cell_data_string();
  void add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
      const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides);
  const std::vector<Variable>& get_variables() const {return variables_;}
};

//--------------------------------------------
//
//--------------------------------------------
class PVTI_file : public File {
private:
  std::string whole_extent;
  long file_position = 0;

public:
  // constructor
  PVTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, const Mpi &_mpi) :
    File(_name, {_path}, ".pvti", _mpi) { set_whole_extent(_geo); }

  std::string get_timestring();
  //void write(const double time, const VTI_file &vti);
  void set_whole_extent(const Geo& geo);
  void MPI_write_piece(const std::string& piece_extent_string, const int rank);
  void write_header(const double time, const VTI_file &vti);
  void write_footer_and_close();
  void set_filename();
  void open();
  void set_position_and_close();
};

//--------------------------------------------
//
//--------------------------------------------
class Outfile {
private:
  PVTI_file pvti_file_;
  VTI_file vti_file_;
  std::vector<Variable> variables_;
  std::vector<int> n;
  //int ***labels_ = nullptr;

public:
  Outfile(std::string &_path, const std::string &_name, const Geo& geo, const Mpi &mpi)
  : pvti_file_(PVTI_file(_path, _name, geo, mpi)),
    vti_file_ (VTI_file (_path, _name, geo, mpi)) { }
    //labels_(geo.labels_) { };

  //void set_cell_data();
  void add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
      const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides);
//  void write(const double time) {
//    vti_file_.write(labels_);
//    pvti_file_.write(time, vti_file_);
//  };
//  void write(int*** labels, const double time){
//    vti_file_.write(labels);
//    pvti_file_.write(time, vti_file_);
//  };
  const std::string& get_filename() const {return pvti_file_.get_filename();}
  VTI_file& get_vti_file() {return(vti_file_);}
  PVTI_file& get_pvti_file() {return(pvti_file_);}

  friend class Variable;
};

//--------------------------------
//
//--------------------------------
class Output {
private:
  std::string path;
  std::vector<Outfile> file;
  std::unordered_map<std::string, int> get_index;
  Mpi mpi_;
  Geo geo_;
  //double time_= 0.0 ;
  int*** labels_ = nullptr;
  int nwrite_ = 0;
  int rank_ = 0;
  int max_rank_ = 0;

public:
  // Constructor
  Output(const std::string _path, const Mpi &mpi, const Geo& geo)
  : path(_path),
    mpi_(mpi),
    geo_(geo),
    labels_(geo.labels_),
    rank_(mpi.get_rank()),
    max_rank_(mpi.get_max_rank()){ }

  Outfile& operator[](const std::string &_name) { return file[get_index[_name]]; }
  Outfile& add_file(const std::string &_name);

  //void set_time(const double time) { time_ = time; }

  void write_all(const double time) {
    for (auto& of:file) {
      write(of, time);
    }
  }

/*
  void push_back_variable() {
    for (std::size_t i=0; i<dims.size(); ++i) {
      variables_.emplace_back(names[i], data_ptrs[i], datasizes[i], dims[i], data_strides[i], n-2*num_ghosts_);
      std::size_t end = variables_.size()-1;
      if (end>0)
        // set the offset in bytes between variables in the same .vti-file
        variables_[end].set_offset(variables_[end-1]);
    }
    set_cell_data_string();
  }
*/

  void write(Outfile& outfile, const double time) {
    VTI_file& vti = outfile.get_vti_file();
    PVTI_file& pvti = outfile.get_pvti_file();
    vti.set_filename_and_open(rank_);
    vti.write_header();
    vti.write_data(labels_);
    vti.write_footer_and_close();
    pvti.set_filename();
    if (rank_ == 0) {
      pvti.open();
      pvti.write_header(time, vti);
      pvti.set_position_and_close();
    }
    pvti.MPI_write_piece(vti.get_piece_extent_string(), rank_);
    if (rank_ == max_rank_) {
      pvti.write_footer_and_close();
    }
    ++(vti.nwrite_);
    ++(pvti.nwrite_);
    ++nwrite_;
  }

  void write(const std::string& var_name, const double time) {
    write(file[get_index[var_name]], time);
  }
  const std::string& get_filename(const std::string& var_name) {
    return(file[get_index[var_name]].get_filename());
  }

  friend class Outfile;
};


#endif /* SRC_OUTPUT_H_ */
