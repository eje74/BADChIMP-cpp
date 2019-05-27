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
#include "Mpi_class.h"
#include "Geo.h"
//#include "global.h"
#include "defines.h"
#include "vector_func.h"
#include "LBgrid.h"
//struct Node;

struct Cells;

//--------------------------------
//
//--------------------------------
//template <typename T>
class Variable {
public:
  std::string name;
  void *data_pointer_ = nullptr;
  const std::vector<double> *data_ptr_;
  const std::vector<std::vector<double>::iterator> data_iter_;
  int datasize = 0;
  int dim = 0;
  int data_stride = 0; // step between data-nodes
  int nbytes = 0;
  int offset = 0;
  int num_nodes = 0;
  std::string datatype;
  std::bitset<NUM_MODE> write_node;
  std::string data_array;

public:
  // constructor
    Variable(const std::string &_name, const std::vector<std::vector<double>::iterator> &data_itr, int _datasize, int _dim)
  : name(_name),
    data_iter_(data_itr),
    datasize(_datasize),
    dim(_dim)
  {
    datatype = "Float" + std::to_string(sizeof(OUTPUT_DTYPE)*8);
    //set_write_nodes({FLUID, WALL, SOLID});
    //set_nbytes(_n);
    set_data_array();
    //std::cout << data_ptr->at(10) << std::endl;
  }

    //  Variable(const std::string &_name, void *data_pointer, int _datasize, int _dim, int _data_stride, const std::vector<int> &_n)
//  : name(_name),
//    data_pointer_(data_pointer),
//    datasize(_datasize),
//    dim(_dim),
//    data_stride(_data_stride)
//  {
//    datatype = "Float" + std::to_string(sizeof(OUTPUT_DTYPE)*8);
//    set_write_nodes({FLUID, WALL, SOLID});
//    set_nbytes(_n);
//    set_data_array();
//  }

//  Variable(const std::string &_name, void *data_pointer, int _datasize, int _dim, int _data_stride)
//  : name(_name),
//    data_pointer_(data_pointer),
//    datasize(_datasize),
//    dim(_dim),
//    data_stride(_data_stride)
//  {
//    datatype = "Float" + std::to_string(sizeof(OUTPUT_DTYPE)*8);
//    set_write_nodes({FLUID, WALL, SOLID});
//    //set_nbytes(_n);
//    set_data_array();
//  }
//
//  Variable(const std::string &_name, const std::vector<double> *data_ptr, int _datasize, int _dim, int _data_stride)
//  : name(_name),
//    data_ptr_(data_ptr),
//    datasize(_datasize),
//    dim(_dim),
//    data_stride(_data_stride)
//  {
//    datatype = "Float" + std::to_string(sizeof(OUTPUT_DTYPE)*8);
//    set_write_nodes({FLUID, WALL, SOLID});
//    //set_nbytes(_n);
//    set_data_array();
//    //std::cout << data_ptr->at(10) << std::endl;
//  }

//  Variable(const std::string &_name, std::vector<double> &data, const std::vector<int> &index, int _datasize, int _dim)
//  : name(_name),
//    data_iter_(index.size()),
//    datasize(_datasize),
//    dim(_dim)
//  {
//    std::vector<double>::iterator start = data.begin();
//    std::cout << std::vector<double>(data.begin(), data.begin()+20) << std::endl;
//    std::cout << std::vector<double>(index.begin(), index.begin()+20) << std::endl;
//    for (int i=0; i<index.size(); ++i) {
//      if (index[i] < data.size()) {
//        data_iter_[i] = data.begin() + index[i];
//      } else {
//        std::cerr << "ERROR in Variables: data index exceeds data size, index = " << index[i] << ", size = " << data.size() << std::endl;
//      }
//    }
//    datatype = "Float" + std::to_string(sizeof(OUTPUT_DTYPE)*8);
//    set_write_nodes({FLUID, WALL, SOLID});
//    //set_nbytes(_n);
//    set_data_array();
//    //std::cout << data_ptr->at(10) << std::endl;
//  }


  void set_nbytes(const std::vector<int> &n);
  void set_offset(const Variable &prev_var);
  void set_write_nodes(const std::vector<int> &_node_mode_to_write);
  void set_data_array();
  template <typename T>
  T get_data(int nn, int dim) const {
    int element = data_stride*nn + dim;
    //int element = nn + dim;
    if (datasize==1) {
      char *char_pointer = static_cast<char*>(data_pointer_);
      return T(char_pointer[element]);
    } else if (datasize==8) {
      double *double_pointer = static_cast<double*>(data_pointer_);
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
  //int num_ghosts_ = 0;
  //Mpi mpi_;
  //std::string path_filename_;

public:
  int nwrite_ = 0;
  File(const std::string &name, const std::vector<std::string> &folders, const std::string &extension)
: name_(name),
  folders_(folders),
  extension_(extension)
  //num_ghosts_(num_ghosts)
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
// VTK Serial XML file for ImageData (structured grid)
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
  //VTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, const Mpi& _mpi)
  //: File(_name, {_path, "vti/"}, ".vti", _mpi), n(_geo.get_n()) { set_extent(_geo); }
  VTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, const int num_ghosts)
  : File(_name, {_path, "vti/"}, ".vti"), n(_geo.get_n()) { set_extent(_geo, num_ghosts); }

  void set_extent(const Geo &geo, const int num_ghosts);
  //void write(int*** labels);
  //void write(int*** labels,const Variable& var);
  void write_data(int*** labels, const int num_ghosts);
  void write_data();
  void write_header();
  void write_footer_and_close();
  void set_filename_and_open(const int rank);
  //std::string get_filename();
  const std::string get_piece_extent_string() const;
  const std::vector<int>& get_system_size() const {return n;}
  const std::string get_cell_data_string() const {return cell_data_string_;}
  void set_cell_data_string();
  void add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
      const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides, const int num_ghosts);
  const std::vector<Variable>& get_variables() const {return variables_;}
};

//--------------------------------------------
// VTK Serial XML file for Unstructured Grid data
//
//  cell_type is set to VTK_VOXEL (=11)
//
//                6        7
//                 o--------o
//              4 /|     5 /|
//               o--------o |       z
//               | | 2    | | 3         y
//               | o--------o       |
//               |/       |/        | /
//               o--------o         |/____ x
//              0        1
//
//--------------------------------------------
class VTU_file : public File {
private:
  std::vector<Variable> variables_;
  //std::vector<int> n;
  std::string type; //, extent;
  std::string cell_data_string_;
  std::string scalar_string_, vector_string_;
  //std::vector<double> points_;
  //std::vector<int> connectivity_;
  int cell_type_ = 11; // VTK_VOXEL
  //int num_points_ = 0, num_cells_ = 0;
  //Node *node_ = nullptr;
  int num_scalar_=0, num_vector_=0;

public:
  //  VTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, Node *node)
  //  : File(_name, {_path, "vti/"}, ".vti", _mpi), n(_geo.get_local_size()), node_(geo.nodes) { set_extent(_geo); }
  //VTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, const Mpi& _mpi)
  //: File(_name, {_path, "vti/"}, ".vti", _mpi), n(_geo.get_n()) { set_extent(_geo); }
  VTU_file(const std::string &_path, const std::string &_name)
  : File(_name, {_path, "vtu/"}, ".vtu") { }

  std::vector<Variable>& get_variables() {return variables_;};

  void set_extent(const Geo &geo, const int num_ghosts);
  //void write(int*** labels);
  //void write(int*** labels,const Variable& var);
  void write_data(int*** labels, const int num_ghosts);
  void write_data();
  void write_header(int, int);
  void write_footer_and_close();
  void set_filename_and_open(const int rank);
  //std::string get_filename();
  const std::string get_piece_string() const;
  //const std::vector<int>& get_system_size() const {return n;}
  const std::string get_cell_data_string() const {return cell_data_string_;}
  void set_cell_data_string();
  void update_cell_data_string(Variable &var);
  void add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
      const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides);
  void add_variable(const std::string &name, const std::vector<double> *data_ptr,
    const int datasize, const int dim, const int data_stride);
  const std::vector<Variable>& get_variables() const {return variables_;}
  void write_grid(std::vector<double> &points, std::vector<int> &connectivity, int num_cells, int);

  //friend class Variable;
};

//--------------------------------------------
// VTK Parallel XML file for ImageData (structured grid)
//--------------------------------------------
class PVTI_file : public File {
private:
  std::string whole_extent;
  long file_position = 0;

public:
  // constructor
//  PVTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, const Mpi &_mpi) :
//    File(_name, {_path}, ".pvti", _mpi) { set_whole_extent(_geo); }
  PVTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, const int num_ghosts) :
    File(_name, {_path}, ".pvti") { set_whole_extent(_geo, num_ghosts); }

  std::string get_timestring();
  //void write(const double time, const VTI_file &vti);
  void set_whole_extent(const Geo& geo, const int num_ghosts);
  void MPI_write_piece(const std::string& piece_extent_string, const int rank);
  void write_header(const double time, const VTI_file &vti);
  void write_footer_and_close();
  void set_filename();
  void open();
  void set_position_and_close();
};

//--------------------------------------------
// VTK Parallel XML file for ImageData (structured grid)
//--------------------------------------------
class PVTU_file : public File {
private:
  std::string whole_extent;
  long file_position = 0;

public:
  // constructor
//  PVTI_file(const std::string &_path, const std::string &_name, const Geo &_geo, const Mpi &_mpi) :
//    File(_name, {_path}, ".pvti", _mpi) { set_whole_extent(_geo); }
  PVTU_file(const std::string &_path, const std::string &_name) :
    File(_name, {_path}, ".pvtu") { }

  std::string get_timestring();
  //void write(const double time, const VTI_file &vti);
  void set_whole_extent(const Geo& geo, const int num_ghosts);
  void MPI_write_piece(const std::string& piece_extent_string, const int rank);
  void write_header(const double time, const VTU_file &vtu);
  void write_footer_and_close();
  void set_filename();
  void open();
  void set_position_and_close();
};


//--------------------------------------------
//
//--------------------------------------------
class Outfile {
protected:
  //PVTI_file pvti_file_;
  //VTI_file vti_file_;
  PVTU_file pvtu_file_;
  VTU_file vtu_file_;
private:
  std::vector<Variable> variables_;
  std::vector<int> n;
  //int num_ghosts_ = 0;
  //int ***labels_ = nullptr;

public:
//  Outfile(std::string &_path, const std::string &_name, const Geo& geo, const int num_ghosts)
//  : pvti_file_(PVTI_file(_path, _name, geo, num_ghosts)),
//    vti_file_ (VTI_file (_path, _name, geo, num_ghosts)),
//    pvtu_file_(PVTU_file(_path, _name)),
//    vtu_file_ (VTU_file (_path, _name)),
//    num_ghosts_(num_ghosts) { }

  Outfile(std::string &_path, const std::string &_name)
  : pvtu_file_(PVTU_file(_path, _name)),
    vtu_file_ (VTU_file (_path, _name)) { }

  //void set_cell_data();
//  void add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
//      const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides);
//  void add_variables(const std::vector<std::string> &names, const std::vector<std::vector<double>&> &data_refs,
//      const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides);
//  void add_variable(const std::string &name, const std::vector<double> data_ptr,
//    const int datasize, const int dim, const int data_stride) {
//    vtu_file_.get_variables().emplace_back(name, data_ptr, datasize, dim, data_stride);
//  }
//  void add_variable(const std::string &name, std::vector<double> &data, const std::vector<int> &index, const int datasize, const int dim) {
//    vtu_file_.get_variables().emplace_back(name, data, index, datasize, dim);
//  }
  void add_variable(const std::string &name, const std::vector<std::vector<double>::iterator> &data_itr, const int datasize, const int dim) {
    vtu_file_.get_variables().emplace_back(name, data_itr, datasize, dim);
    vtu_file_.update_cell_data_string(vtu_file_.get_variables().back());
  }

//  void write(const double time) {
//    vti_file_.write(labels_);
//    pvti_file_.write(time, vti_file_);
//  };
//  void write(int*** labels, const double time){
//    vti_file_.write(labels);
//    pvti_file_.write(time, vti_file_);
//  };
//  const std::string& get_filename() const {return pvti_file_.get_filename();}
//  VTI_file& get_vti_file() {return(vti_file_);}
//  PVTI_file& get_pvti_file() {return(pvti_file_);}
  const std::string& get_filename() const {return pvtu_file_.get_filename();}
  VTU_file& get_vtu_file() {return(vtu_file_);}
  PVTU_file& get_pvtu_file() {return(pvtu_file_);}

  //friend class Variable;
};


//--------------------------------
//
//--------------------------------
class Output {
private:
  std::string path;
  std::vector<Outfile> outfiles_;
  std::unordered_map<std::string, int> get_index;
  // VTK_VOXEL
  const std::vector<std::vector<int>> cell_points_3D_ = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
  // VTK_PIXEL
  const std::vector<std::vector<int>> cell_points_2D_ = {{0,0}, {1,0}, {0,1}, {1,1}};
  int nwrite_ = 0;
  int rank_ = 0;
  int max_rank_ = 0;
  std::vector<int> dim_;
  std::vector<double> points_;
  std::vector<int> connectivity_;
  int num_cell_points_ = 0;
  int num_points_ = 0, num_cells_ = 0;
  double cell_size_ = 1.0;

public:
  // Constructor
  //Output(const std::string _path, const Mpi &mpi, const Geo& geo)
//  Output(int dim, const std::string _path, const int rank, const int max_rank, int *p0, int *p1, std::vector<int> connectivity) //, int *c0, int *c1)
//  : path(_path),
//    //mpi_(mpi),
//    //geo_(geo),
//    points_(p0,p1),
//    connectivity_(connectivity),
//    //labels_(geo.labels_),
//    rank_(rank),
//    max_rank_(max_rank),
//    dim_(dim)
//  {
//    num_points_ = int(points_.size()/dim_);
//    num_cells_ = int(connectivity_.size()/num_cell_points_);
//    //if (rank==0)
//    //  std::cout << "Points: " << points_ << std::endl;
//  }
//  Output(int dim, const std::string _path, const int rank, const int max_rank, std::vector<double>& points, std::vector<int>& connectivity) //, int *c0, int *c1)
//  : path(_path),
//    //mpi_(mpi),
//    //geo_(geo),
//    points_(points),
//    connectivity_(connectivity),
//    //labels_(geo.labels_),
//    rank_(rank),
//    max_rank_(max_rank),
//    dim_(dim)
//  {
//    num_points_ = int(points_.size()/dim_.size());
//    num_cells_ = int(connectivity_.size()/num_cell_points_);
//    //if (rank==0)
//    //  std::cout << "Points: " << points_ << std::endl;
//  }

  Output(std::vector<int> dim, const std::string _path, const int rank, const int max_rank, std::vector<std::vector<int>>& node_pos)
  : path(_path),
    rank_(rank),
    max_rank_(max_rank),
    dim_(dim)
  {
    set_points_and_connectivity(node_pos);

    num_points_ = int(points_.size()/dim_.size());
    num_cells_ = int(connectivity_.size()/num_cell_points_);
  }

  void set_points_and_connectivity(std::vector<std::vector<int>>& node_pos);
  Outfile& operator[](const std::string &_name) { return outfiles_[get_index[_name]]; }
  Outfile& add_file(const std::string &_name);

  //void set_time(const double time) { time_ = time; }

//  void write_all(const double time) {
//    for (auto& of:file) {
//      write(of, time);
//    }
//  }

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

  //void write_unstructured(Outfile& outfile, const double time) {
  void write(Outfile& outfile, const double time) {
    VTU_file& vtu = outfile.get_vtu_file();
    PVTU_file& pvtu = outfile.get_pvtu_file();
    vtu.set_filename_and_open(rank_);
    vtu.write_header(num_points_, num_cells_);
    vtu.write_grid(points_, connectivity_, num_cells_, num_cell_points_);
    vtu.write_data();
    //vti.write_data();
    vtu.write_footer_and_close();
    pvtu.set_filename();
    if (rank_ == 0) {
      pvtu.open();
      pvtu.write_header(time, vtu);
      pvtu.set_position_and_close();
    }
    pvtu.MPI_write_piece(vtu.get_piece_string(), rank_);
    if (rank_ == max_rank_) {
      pvtu.write_footer_and_close();
    }
    ++(vtu.nwrite_);
    ++(pvtu.nwrite_);
    ++nwrite_;
  }

//  void write(Outfile& outfile, const double time) {
//    VTI_file& vti = outfile.get_vti_file();
//    PVTI_file& pvti = outfile.get_pvti_file();
//    vti.set_filename_and_open(rank_);
//    vti.write_header();
//    vti.write_data(labels_, num_ghosts_);
//    //vti.write_data();
//    vti.write_footer_and_close();
//    pvti.set_filename();
//    if (rank_ == 0) {
//      pvti.open();
//      pvti.write_header(time, vti);
//      pvti.set_position_and_close();
//    }
//    pvti.MPI_write_piece(vti.get_piece_extent_string(), rank_);
//    if (rank_ == max_rank_) {
//      pvti.write_footer_and_close();
//    }
//    ++(vti.nwrite_);
//    ++(pvti.nwrite_);
//    ++nwrite_;
//  }

  void write(const std::string& var_name, const double time) {
    write(outfiles_[get_index[var_name]], time);
  }

  const std::string& get_filename(const std::string& var_name) {
    return(outfiles_[get_index[var_name]].get_filename());
  }

  friend class Outfile;
};


#endif /* SRC_OUTPUT_H_ */
