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
//#include "vector_func.h"
#include "LBgrid.h"
//struct Node;

//struct Cells;

//--------------------------------------------
// Grid class for VTK Unstructured Grid data
//
//  cell_type is set to VTK_VOXEL (=11) in 3D
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
// and VTK_PIXEL in 2D
//--------------------------------------------
class VTKGrid {
private:
  // VTK_VOXEL
  const std::vector<std::vector<int>> cell_points_3D_ = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
  const int cell_type_3D_ = 11;
  // VTK_PIXEL
  const std::vector<std::vector<int>> cell_points_2D_ = {{0,0}, {1,0}, {0,1}, {1,1}};
  const int cell_type_2D_ = 8;
  std::vector<int> dim_;
  std::vector<double> points_;
  std::vector<int> connectivity_;
  int num_cell_points_ = 0;
  int num_cells_ = 0;
  int num_points_ = 0;
  double cell_edge_ = 1.0;
  int cell_type_ = 0;

public:
  VTKGrid(std::vector<int> dim, std::vector<std::vector<int>>& node_pos) : dim_(dim) {
    cell_type_ = (dim_.size()>2)? cell_type_3D_ : cell_type_2D_;
    calculate_points_and_connectivity(node_pos);
    set_num_points();
    set_num_cells();
  };
  const std::string get_point_list() const { return vector_as_string_list(points_); }
  const std::string get_connectivity_list() const { return vector_as_string_list(connectivity_); }
  const std::string get_offset_list() {return incremental_list<int>(num_cell_points_, connectivity_.size(), num_cell_points_); }
  const std::string get_types_list() { return repeated_list<int>(cell_type_, num_cells_); }
  int get_num_points() const { return num_points_ ; };
  int get_num_cells() const { return num_cells_ ; };

private:
  void set_num_points() { num_points_ = int(points_.size()/dim_.size()); };
  void set_num_cells() { num_cells_ = int(connectivity_.size()/num_cell_points_); };
  void calculate_points_and_connectivity(std::vector<std::vector<int>>& node_pos) {
    // set 2D or 3D cell
    std::vector<std::vector<int>> cell_points = (dim_.size()>2)? cell_points_3D_ : cell_points_2D_;
    num_cell_points_ = cell_points.size();
    connectivity_.reserve(cell_points.size()*node_pos.size());

    std::vector<int> point_index(prod(dim_+1), -1);
    for (const auto& n:node_pos) {
      for (const auto& c:cell_points) {
        std::vector<int> p = n+c;
        // give cell-points a unique index
        int idx = p[0] + p[1]*dim_[0] + p[2]*dim_[0]*dim_[1];
        if (point_index[idx]<0) {
          points_.insert(points_.end(), std::begin(p), std::end(p));
          point_index[idx] = int(points_.size()/dim_.size())-1;
        }
        connectivity_.push_back(point_index[idx]);
      }
    }
    // offset cell (corner) points by half the cell-size
    points_ = points_ - 0.5*cell_edge_;
  }
};

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

public:
  int nwrite_ = 0;
  File(const std::string &name, const std::vector<std::string> &folders, const std::string &extension)
: name_(name),
  folders_(folders),
  extension_(extension)
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
  int cell_type_ = 11; // VTK_VOXEL
  int num_scalar_=0, num_vector_=0;

public:
  VTU_file(const std::string &_path, const std::string &_name)
  : File(_name, {_path, "vtu/"}, ".vtu") { }

  std::vector<Variable>& get_variables() {return variables_;};

  void set_extent(const Geo &geo, const int num_ghosts);
  void write_data(int*** labels, const int num_ghosts);
  void write_data();
  void write_header(int, int);
  void write_header(VTKGrid& grid);
  void write_footer_and_close();
  void set_filename_and_open(const int rank);
  const std::string get_piece_string() const;
  const std::string get_cell_data_string() const {return cell_data_string_;}
  void set_cell_data_string();
  void update_cell_data_string(Variable &var);
  void add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
      const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides);
  void add_variable(const std::string &name, const std::vector<double> *data_ptr,
    const int datasize, const int dim, const int data_stride);
  const std::vector<Variable>& get_variables() const {return variables_;}
  void write_grid(std::vector<double> &points, std::vector<int> &connectivity, int num_cells, int);
  void write_grid(VTKGrid& grid);
};


//--------------------------------------------
// VTK Parallel XML file for Unstructured Grid
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
  PVTU_file pvtu_file_;
  VTU_file vtu_file_;
private:
  std::vector<Variable> variables_;
  std::vector< std::pair<std::string,int> > var_dim_;
  //std::vector<int> n;

public:
  Outfile(std::string &_path, const std::string &_name)
  : pvtu_file_(PVTU_file(_path, _name)),
    vtu_file_ (VTU_file (_path, _name)) { }

  void add_variable(const std::string &name, const std::vector<std::vector<double>::iterator> &data_itr, const int datasize, const int dim) {
    vtu_file_.get_variables().emplace_back(name, data_itr, datasize, dim);
    vtu_file_.update_cell_data_string(vtu_file_.get_variables().back());
  }

  const std::string& get_filename() const {return pvtu_file_.get_filename();}
  VTU_file& get_vtu_file() {return(vtu_file_);}
  PVTU_file& get_pvtu_file() {return(pvtu_file_);}
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
  VTKGrid grid_;

public:
  // Constructor

  Output(std::vector<int> dim, const std::string _path, const int rank, const int max_rank, std::vector<std::vector<int>>& node_pos)
  : path(_path),
    rank_(rank),
    max_rank_(max_rank),
    grid_(dim, node_pos),
    dim_(dim)
  {
    set_points_and_connectivity(node_pos);

    num_points_ = int(points_.size()/dim_.size());
    num_cells_ = int(connectivity_.size()/num_cell_points_);
    //grid_.set_num_points();
    //grid_.set_num_cells();
  }

  void set_points_and_connectivity(std::vector<std::vector<int>>& node_pos);
  Outfile& operator[](const std::string &_name) { return outfiles_[get_index[_name]]; }
  Outfile& add_file(const std::string &_name);

  void write(Outfile& outfile, const double time) {
    VTU_file& vtu = outfile.get_vtu_file();
    PVTU_file& pvtu = outfile.get_pvtu_file();
    vtu.set_filename_and_open(rank_);
    vtu.write_header(grid_);
    vtu.write_grid(grid_);
    //vtu.write_header(num_points_, num_cells_);
    //vtu.write_grid(points_, connectivity_, num_cells_, num_cell_points_);
    vtu.write_data();
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


  void write(const std::string& var_name, const double time) {
    write(outfiles_[get_index[var_name]], time);
  }

  const std::string& get_filename(const std::string& var_name) {
    return(outfiles_[get_index[var_name]].get_filename());
  }

  friend class Outfile;
};


#endif /* SRC_OUTPUT_H_ */
