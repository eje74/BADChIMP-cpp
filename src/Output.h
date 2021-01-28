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
//#include "Mpi_class.h"
#include <mpi.h>
#include "Geo.h"
//#include "global.h"
#include "defines.h"
//#include "vector_func.h"
#include "LBgrid.h"
#include "LBnodes.h"
#include "LBfield.h"
#include "LBvtk.h"

//template <typename T1> using datatype = std::valarray<T1>;
using data_array = std::valarray<double>;

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
    // set cell to 2D or 3D
    cell_type_ = (dim_.size()>2)? cell_type_3D_ : cell_type_2D_;
    calculate_points_and_connectivity(node_pos);
    set_num_points();
    set_num_cells();
  };
  const std::string get_point_list() const { return vector_as_string_list(points_); }
  const std::string get_connectivity_list() const { return vector_as_string_list(connectivity_); }
  const std::string get_offset_list() const {return incremental_list<int>(num_cell_points_, connectivity_.size(), num_cell_points_); }
  const std::string get_types_list() const { return repeated_list<int>(cell_type_, num_cells_); }
  int get_num_points() const { return num_points_ ; };
  int get_num_cells() const { return num_cells_ ; };

private:
  void set_num_points() { num_points_ = int(points_.size()/dim_.size()); };
  void set_num_cells() { num_cells_ = int(connectivity_.size()/num_cell_points_); };
  void calculate_points_and_connectivity(std::vector<std::vector<int>>& node_pos);
};

//--------------------------------
//
//--------------------------------
//template <typename T>
class Variable {
public:
  std::string name;
  const data_array& data_;
  const std::vector<int> index_;
  int dim = 0;
  std::string type = "Float64";

public:
//  // constructor
    Variable(const std::string &_name, const data_array& data, const std::vector<int>& index, int _dim)
  : name(_name),
    data_(data),
    index_(index),
    dim(_dim)
  {
      //type = "Float" + std::to_string(sizeof(OUTPUT_DTYPE)*8);
  }
};

//--------------------------------------------
// Base class for the .vtu and .pvtu file classes
//--------------------------------------------
class File {
protected:
  std::string name_;
  std::ofstream file_;
  std::string filename_;
  std::vector<std::string> folders_;
  std::string path_;
  std::string extension_;
  const int precision_ = 5;

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
  void open() {file_.open(path_+filename_, std::ios::out);}
  const std::string& get_filename() const {return filename_;}

private:
  void make_dir(std::string &dir);
};


//--------------------------------------------
// VTK Serial XML file for Unstructured Grid data
//--------------------------------------------
class VTU_file : public File {
private:
  std::string cell_data_string_;
  std::string scalar_string_, vector_string_;
  int num_scalar_=0, num_vector_=0;

public:
  VTU_file(const std::string &_path, const std::string &_name) : File(_name, {_path, "vtu/"}, ".vtu") { }
  void write_data(const std::vector<Variable>& variables);
  void write_header(const VTKGrid& grid);
  void write_footer_and_close();
  void set_filename_and_open(const int rank);
  const std::string get_piece_string() const;
  void update_cell_data_string(Variable &var);
};


//--------------------------------------------
// VTK Parallel XML file for Unstructured Grid
//--------------------------------------------
class PVTU_file : public File {
private:
  long file_position = 0;

public:
  // constructor
  PVTU_file(const std::string &_path, const std::string &_name) : File(_name, {_path}, ".pvtu") { }

  std::string get_timestring();
  void MPI_write_piece(const std::string& piece_extent_string, const int rank);
  void write_header(const double time, const std::vector<Variable> &varlist);
  void write_footer_and_close();
  void set_filename();
  void set_position_and_close();
};


//--------------------------------------------
//
//--------------------------------------------
class Outfile {
private:
  PVTU_file pvtu_file_;
  VTU_file vtu_file_;
  std::vector<Variable> variables_;

public:
  Outfile(std::string &_path, const std::string &_name)
  : pvtu_file_(PVTU_file(_path, _name)),
    vtu_file_ (VTU_file (_path, _name)) { }

//  void add_variable(const std::string &name, const std::vector<double> &data, const std::vector<int> index, const int dim) {
  void add_variable(const std::string &name, const data_array &data, const std::vector<int> index, const int dim) {
    variables_.emplace_back(name, data, index, dim);
    vtu_file_.update_cell_data_string(variables_.back());
  }
  const std::vector<Variable>& get_variables() const {return variables_;}
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
  int nwrite_ = 0;
  int rank_ = 0;
  int max_rank_ = 0;
  VTKGrid grid_;

public:
  // Constructor
  Output(std::vector<int> dim, const std::string _path, const int rank, const int num_procs, std::vector<std::vector<int>>& node_pos)
  : path(_path),
    rank_(rank),
    max_rank_(num_procs-1),
    grid_(dim, node_pos)
    { }

  Outfile& operator[](const std::string &_name) { return outfiles_[get_index[_name]]; }
  Outfile& add_file(const std::string &_name);

  void write(Outfile& outfile, const double time);
  inline void write(const std::string& var_name, const double time) {
    write(outfiles_[get_index[var_name]], time);
  }
  const std::string& get_filename(const std::string& var_name) {
    return(outfiles_[get_index[var_name]].get_pvtu_file().get_filename());
  }

};

/* ********************************************************************** *
 *                                                                        *
 *          FUNCTIONS                                                     *
 *                                                                        *
 *                                                                        *
 * ********************************************************************** */
template<typename DXQY, typename T>
void outputStdVector(const std::string &fieldName, const std::vector<T> &scalarField, const std::string &outputDir, 
    const int &myRank, const int &nProcs, const Grid<DXQY> &grid, const LBvtk<DXQY> &vtklb)
{
    auto node_pos_all = grid.getNodePos(vtklb.beginNodeNo(), vtklb.endNodeNo());
    auto globalDim = vtklb.getGlobaDimensions(); // Set as default to 3 dimensions, as prescribed by the Output class

    // Setup allNodes
    Output output(globalDim, outputDir, myRank, nProcs, node_pos_all);

    ScalarField val(1, grid.size());
    std::vector<int> allNodes(vtklb.endNodeNo()-vtklb.beginNodeNo());
    int cnt = 0;
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        val(0, nodeNo) = scalarField[nodeNo];
        allNodes[cnt] = nodeNo;
        cnt++;
    }

    // Write field to file
    output.add_file(fieldName);
    output[fieldName].add_variable(fieldName, val.get_data(), val.get_field_index(0, allNodes), 1); 
    output.write(fieldName, 0);
}

template<typename DXQY>
void outputGeometry(const std::string &fileName, const std::string &outputDir, const int &myRank, const int &nProcs, 
    const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const LBvtk<DXQY> &vtklb)
{
    
    std::vector<int> val(grid.size());
    val[0] = -1;
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        val[nodeNo] = nodes.isSolid(nodeNo) ? 1 : 0;        
    }
    
    outputStdVector(fileName, val, outputDir, myRank, nProcs, grid, vtklb);

}


#endif /* SRC_OUTPUT_H_ */
