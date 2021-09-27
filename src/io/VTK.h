
#ifndef SRC_VTK_H_
#define SRC_VTK_H_
#include <unordered_map>
#include <map>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <bitset>
#include <algorithm>
#include <mpi.h>
//#include <sys/types.h>
#include <sys/stat.h>
//#if defined (_WIN32)
//#include <direct.h>
//#else
#include <unistd.h>
//#endif
#include "vector_func.h"

namespace VTK {
  
  using node_dtype = float;
  using output_dtype = double;

  static constexpr int BINARY = 0;
  static constexpr int ASCII = 1;

  //--------------------------------------------
  //  Get datatype name
  //  Adapted from: https://stackoverflow.com/questions/4484982/how-to-convert-typename-t-to-string-in-c
  //--------------------------------------------
  template <typename T>
  struct datatype {
    static const char* name() {
      return typeid(T).name();
    } 
  };
  template <>
  struct datatype<int> {
    static const char* name() {
      return "Int32";
    } 
  };
  template <>
  struct datatype<unsigned int> {
    static const char* name() {
      return "UInt32";
    } 
  };
  template <>
  struct datatype<float> {
    static const char* name() {
      return "Float32";
    } 
  };
  template <>
  struct datatype<double> {
    static const char* name() {
      return "Float64";
    } 
  };


  //--------------------------------------------
  // 1: VTK_VERTEX
  //            o
  //--------------------------------------------
  struct vertex {
    static constexpr int type = 1;
    static constexpr int dim = 1;
    static constexpr int n = 1;
    static constexpr std::array<std::array<int, dim>, n> points = { {{0}} };
  };
  constexpr std::array<std::array<int, vertex::dim>, vertex::n> vertex::points;
  constexpr int vertex::type;
  //--------------------------------------------
  // 3: VTK_LINE
  //            o
  //             \
  //              \
  //               \
  //                o
  //-------------------------------------------- 
  struct line {
    static constexpr int type = 3;
    static constexpr int dim = 1;
    static constexpr int n = 2;
    static constexpr std::array<std::array<int, dim>, n> points = { {{0}, {1}} };
  };
  constexpr std::array<std::array<int, line::dim>, line::n> line::points;
  constexpr int line::type;
  //-------------------------------------------- 
  //  8: VTK_PIXEL
  //
  //             2 o--------o 3
  //               |        | 
  //               |        |         y
  //               |        |         | 
  //             0 o--------o 1       |____ x
  //                      
  //-------------------------------------------- 
  struct pixel {
    static constexpr int type = 8;
    static constexpr int dim = 2;
    static constexpr int n = 4;
    static constexpr std::array<std::array<int, dim>, n> points = {{{0,0}, {1,0}, {0,1}, {1,1}}};
  };
  constexpr std::array<std::array<int, pixel::dim>, pixel::n> pixel::points;
  constexpr int pixel::type;
  //-------------------------------------------- 
  //  9: VTK_QUAD
  //                   2 
  //                  o 
  //                 /  \
  //                /    \
  //               /      \
  //            3 o        o 1
  //               \      /
  //                \    /  
  //                 \  /
  //               0  o     
  //-------------------------------------------- 
  struct quad {
    static constexpr int type = 9;
    static constexpr int dim = 2;
    static constexpr int n = 4;
    static constexpr std::array<std::array<int, dim>, n> points = {{{0,0}, {1,0}, {1,1}, {0,1}}};
  };
  constexpr std::array<std::array<int, quad::dim>, quad::n> quad::points;
  constexpr int quad::type;
  //--------------------------------------------
  //  11: VTK_VOXEL
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
  //--------------------------------------------
  struct voxel {
    static constexpr int type = 11;
    static constexpr int dim = 3;
    static constexpr int n = 8;
    static constexpr std::array<std::array<int, dim>, n> points = {{{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}}};
  };
  constexpr std::array<std::array<int, voxel::dim>, voxel::n> voxel::points;
  constexpr int voxel::type;

  // ============================================
  //                M E S H
  // ============================================
  template <typename T>
  class Mesh
  {
  private:
    std::vector<node_dtype> points_;
    std::vector<int> conn_;
    std::vector<int> offsets_;
    std::vector<int> types_;
    //std::vector<int> size_;
    //std::vector<node_dtype> max_;
    //std::vector<node_dtype> min_;
    node_dtype unit_;
     //int n_points_ = 0;

  public:
    //--------------------------------------------
    Mesh(std::vector<node_dtype> &points, std::vector<int> &conn, node_dtype unit = 1) 
      : points_(points), conn_(conn), offsets_(), types_(), unit_(unit)
    //--------------------------------------------
    {
      offsets_ = linspace<int>(T::n, conn_.size(), T::n);
      types_ = std::vector<int>(num(), T::type);
    } 

    //--------------------------------------------
    Mesh(std::vector<std::vector<node_dtype>> &nodes, node_dtype unit = 1) 
      : points_(), conn_(), offsets_(), types_(), unit_(unit)
    //--------------------------------------------
    {
      calc_points_and_conn(nodes);
      offsets_ = linspace<int>(T::n, conn_.size(), T::n);
      types_ = std::vector<int>(num(), T::type);
      std::cout << "type, dim, n: " << T::type << ", " << T::dim << ", " << T::n << std::endl;
      std::cout << "points: "; print_points();
      std::cout << "connectivity: " << conn_ << std::endl;
      std::cout << "offsets: " << offsets_ << std::endl;
      std::cout << "types: " << types_ << std::endl;
    }
    const std::vector<node_dtype>& points() const { return points_; }
    const std::vector<int>& connectivity() const { return conn_; }
    const std::vector<int>& offsets() const { return offsets_; }
    const std::vector<int>& types() const { return types_; }
    const int num() const { return int(conn_.size()/T::n); }
    void print_points() 
    { 
      std::cout << "(";
      int n = points_.size();
      for (int i=0; i<n; ++i) {
          //std::cout << points_[i] << (( i<n-1 && (i+1)%T::dim==0) ? "),(" : ",");
          std::cout << points_[i] << (( i<n-1 && (i+1)%3==0) ? "),(" : ",");
      }
      std::cout << ")" << std::endl;
    }

  private:
    //--------------------------------------------
    std::vector<int> get_stride(std::vector<std::vector<node_dtype>> &nodes)
    //--------------------------------------------
    {
      std::vector<node_dtype> min_n(T::dim, 9e9), max_n(T::dim, -1);
      for (const auto &node : nodes) {
        for (auto i=0; i < T::dim; i++) {
          if (node[i] > max_n[i])
            max_n[i] = node[i];
          if (node[i] < min_n[i])
            min_n[i] = node[i];
        }
      }
      std::vector<int> stride_vec(T::dim, 1); 
      for (auto i = 1; i < T::dim; i++) {
        int size = int(max_n[i] - min_n[i]) + 2;
        stride_vec[i] = stride_vec[i-1] * size;
      }
      return stride_vec;
    }

    //--------------------------------------------
    void calc_points_and_conn(std::vector<std::vector<node_dtype>> &nodes) 
    //--------------------------------------------
    {
      // Loop over nodes and calculate corner points and their unique index 
      auto stride = get_stride(nodes);
      std::vector<double> pts, index;
      pts.reserve(T::n * nodes.size() * T::dim);
      index.reserve(T::n * nodes.size());
      for (const auto &node : nodes) {
        for (const auto &cell_point : T::points) {
          int idx = 0;
          for (auto i=0; i < T::dim; i++) {
            auto p = node[i] + cell_point[i];
            idx += p*stride[i];
            pts.push_back(p);
          }
          index.push_back(idx);
        }
      }
      //std::cout << index << std::endl;
      // Create unique index_list 
      auto minmax = std::minmax_element(index.begin(), index.end());
      auto min = *minmax.first;
      auto max = *minmax.second;
      //std::cout << index << std::endl;
      std::vector<int> index_list(max - min + 1, -1);
      conn_.reserve(T::n * nodes.size());
      int npts = -1;
      for (int i=0; i<index.size(); ++i) {
        int idx = index[i] - min;
        //std::cout << idx << std::endl;
        //std::cout << index_list << std::endl;
        if (index_list[idx] < 0) {
          ++npts;
          auto pos = pts.begin() + T::dim*i;
          points_.insert(points_.end(), pos, pos + T::dim);
          // Add 0 for dimensions less than 3
          for (int k=3-T::dim; k>0; --k) {
            points_.push_back(0.5*unit_);
          }
          index_list[idx] = npts;
        }
        //std::cout << idx << std::endl;
        //std::cout << index_list << std::endl;
        conn_.push_back(index_list[idx]);
      }
      // Shift position so that the node-point is in the centre of the cell
      for (auto& p : points_)
        p -= 0.5*unit_;
    }
  };

  //============================================
  //            D A T A A R R A Y
  //============================================
  class DataArray {
    private:
    std::vector<std::string> data_formats_ = {"appended", "ascii"};
    std::map<std::string, std::string> tag_var_  = { 
      {"name", "Name"}, {"type", "type"}, {"dim", "NumberOfComponents"}, {"format", "format"}, {"offset", "offset"}
    };
    std::map<std::string, std::string> tags_;
    int format_ = 0;

    public:
    //--------------------------------------------
    DataArray(const std::string& name, int dim, const int format, const std::string& type="", int offset=-1) : tags_(), format_(format) 
    {
      set_tag("name", name);
      set_tag("type", type);
      //set_tag("dim", (dim>1) ? 3 : -1);
      set_tag("dim", dim);
      set_tag("format", data_formats_[format]);
      set_tag("offset", (format==ASCII) ? -1 : offset );
    }
    //--------------------------------------------
    std::string quote(const std::string& str) const { return "=\""+str+"\" "; }
    //--------------------------------------------
    void set_tag(const std::string& tname, int val) { 
      tags_[tname] = (val>=0) ? tag_var_[tname]+quote(std::to_string(val)) : ""; 
    }
    //--------------------------------------------
    void set_tag(const std::string& tname, const std::string& val) { 
      tags_[tname] = (val.empty()) ? "" : tag_var_[tname]+quote(val); 
    }
    //--------------------------------------------
    const std::string& operator[](const std::string& tname) const { return tags_.at(tname); }
    //--------------------------------------------
    friend std::ostream& operator<<(std::ostream& out, const DataArray& da) {
      out << da["name"] << da["type"] << da["dim"] << da["format"] << da["offset"];
      return out;
    }
  };

  //============================================
  //               D A T A
  //============================================
  template <typename T>
  class Data {
    public:
    std::string name_ = "";
    const std::vector<T>& data_;
    int dim_ = 0;
    unsigned int offset_ = 0;
    unsigned int length_ = 0;
    DataArray dataarray_;
    unsigned int nbytes_ = 0;
    std::vector<int> index_;

    public:
    //--------------------------------------------
    Data() : data_(std::vector<T>()) {};
    //--------------------------------------------
    //--------------------------------------------
    Data(const std::string &name, const std::vector<T>& data, const int format, const int dim, const int offset=0, const int length=0) 
      : name_(name), data_(data), dim_(dim), offset_(offset), length_(length), dataarray_(name, dim, format)  
    //--------------------------------------------
    { 
      dataarray_.set_tag("type", datatype<T>::name());
      if (length == 0) {
        length_ = data.size();
      }
      if (format == BINARY) {
        nbytes_ = length_*sizeof(T);
      }
    }
    //--------------------------------------------
    Data(const std::string &name, const std::vector<T>& data, const int format, const int dim, std::vector<int>& index, const int offset=0, const int length=0) 
      : Data(name, data, format, dim, offset, length) 
    //--------------------------------------------
    { 
      index_ = index;
    }  

    // //--------------------------------------------
    // void set_data_index(std::vector<int>& ind) { 
    //   index_ = ind;
    //   length_ = ind.size();
    //   if (format == BINARY) {
    //     nbytes_ = length_*sizeof(T);
    //   }
    // } 
    //--------------------------------------------
    int size() const { return data_.size(); } 
    //--------------------------------------------
    //--------------------------------------------
    int dim() const { return dim_; } 
    //--------------------------------------------
    //--------------------------------------------
    bool is_binary() const {return (nbytes_ > 0) ? true : false; }
    //--------------------------------------------
    //--------------------------------------------
    void write_asciidata (std::ofstream& file, double min=-1) const {
    //--------------------------------------------
      if (!is_binary())
        file << vector_as_string_list(data_, offset_, offset_+length_, min);
    }
    //--------------------------------------------
    void write_binarydata(std::ofstream& file) const {
    //--------------------------------------------
      if (is_binary()) {
        file.write((char*)&nbytes_, sizeof(unsigned int));
        file.write((char*)&data_[offset_], nbytes_);
      }
    }
    //--------------------------------------------
    const std::string& dataarray(const std::string& tname) { return dataarray_[tname]; }
    //--------------------------------------------
    //--------------------------------------------
    int update_offset(int offset) { 
    //--------------------------------------------
      if (is_binary() && dataarray_["offset"].empty()) {
        dataarray_.set_tag("offset", offset); 
        return offset + nbytes_ + sizeof(unsigned int);
      } else {
        return offset;
      }
    }
    //--------------------------------------------
    void write_dataarray(std::ofstream& file, double min=-1) const {
    //--------------------------------------------
      file << "<DataArray " << dataarray_ << ">";
      write_asciidata(file, min);  
      file << "</DataArray>" << std::endl;
    }
  };

  //============================================
  //                G R I D
  //============================================
  template <typename T>
  class Grid {
    private:
    static const int point_dim_ = 3;
    Mesh<T> mesh_;
    Data<node_dtype> point_data_;
    std::vector<Data<int>> cell_data_;
    unsigned int offset_ = 0;

    public:
    //--------------------------------------------
    Grid() : mesh_(), point_data_(), cell_data_() { }
    //--------------------------------------------
    //--------------------------------------------
    Grid(std::vector<std::vector<node_dtype>>& nodes, const int format=BINARY) 
    : mesh_(nodes), point_data_("", mesh_.points(), format, point_dim_), cell_data_()
    //--------------------------------------------
    {
      cell_data_.emplace_back("connectivity", mesh_.connectivity(), format, 1);
      cell_data_.emplace_back("offsets",      mesh_.offsets()     , format, 1);
      cell_data_.emplace_back("types",        mesh_.types()       , format, 1);
      update_offset();
    }
    //--------------------------------------------
    void update_offset()
    //--------------------------------------------
    { 
      offset_ = point_data_.update_offset(offset_);
      for (auto& cell : cell_data_) 
        offset_ = cell.update_offset(offset_);
    }
    //--------------------------------------------
    int offset() { return offset_; }
    //--------------------------------------------
    //--------------------------------------------
    int num_points() { return int(point_data_.size()/point_data_.dim()); }
    //--------------------------------------------
    //--------------------------------------------
    int num_cells() { return cell_data_.back().size(); };
    //--------------------------------------------
    //--------------------------------------------
    Data<node_dtype>& point_data() { return point_data_; };
    //--------------------------------------------
    //--------------------------------------------
    std::vector<Data<int>>& cell_data() { return cell_data_; };
    //--------------------------------------------
  };

  //============================================
  //           V A R I A B L E S
  //============================================
  template <typename T>
  class Variables {
    public: 
    std::vector<Data<T>> datalist_;
    std::string cell_data_string_ = "";
    std::string scalar_string_="", vector_string_="";
    int num_scalar_=0, num_vector_=0;

    // Constructor
    //--------------------------------------------
    Variables() : datalist_() { }    
    //--------------------------------------------
    void add(const std::string& name, const std::vector<T>& data, const int format, const int dim, std::vector<int>& index, const int offset=0, const int length=0) 
    {
      datalist_.emplace_back(name, data, format, dim, index, offset, length);
      update_cell_data_string();
    }
    //--------------------------------------------
    std::vector<Data<T>>& data() { return datalist_; }
    //--------------------------------------------
    Data<T>& back() { return datalist_.back(); }

    private:
    //--------------------------------------------
    void update_cell_data_string() {
      // make lists of scalars and vectors
      std::string sep = "";
      const auto& var = datalist_.back();
      if (var.dim_>1) {
        // vector
        if (num_vector_>0) {
          sep = " ,";
        }
        ++num_vector_;
        vector_string_ += sep + var.name_;
      } else {
        // scalar
        if (num_scalar_>0) {
          sep = " ,";
        }
        ++num_scalar_;
        scalar_string_ += sep + var.name_;
      }
    }

  };

  //============================================
  //                F I L  E
  //                
  // Base class for the .vtu and .pvtu file classes
  //============================================
  class File {
    protected:
    std::string name_ = "";
    std::ofstream file_;
    std::string filename_ = "";
    std::vector<std::string> folders_;
    std::string path_ = "";
    std::string extension_ = "";
    const int precision_ = 5;

    public:
    int nwrite_ = 0;
    //--------------------------------------------
    File() : file_(), folders_() { }
    //--------------------------------------------
    //--------------------------------------------
    File(const std::string &name, const std::vector<std::string> &folders, const std::string &extension)
    : name_(name), file_(), folders_(folders), extension_(extension)
    //--------------------------------------------
    {
      for (const auto& f : folders_) {
        path_ += f;
        make_dir(path_);
      }
    }
    //--------------------------------------------
    void open() { file_.open(path_+filename_, std::ios::out); }
    //--------------------------------------------
    void close() {
      file_ << "</VTKFile>" << std::endl;
      file_.close();
    }
    //--------------------------------------------
    std::string& filename() { return filename_;}
    //--------------------------------------------
    std::ofstream& file(){ return file_; }

    private:
    //--------------------------------------------
    void make_dir(std::string &dir) {
      // check that dir has a '/' at the end
      if (dir.back() != '/') {
        dir += '/';
      }
      struct stat st = {0};
      if (stat(dir.c_str(), &st) == -1) {
    #ifdef _WIN32
          _mkdir(dir.c_str());
    #else
          mkdir(dir.c_str(), 0700);
    #endif
      }
    }
  };


  //============================================
  //               V T U F I L E
  //                
  // VTK Serial XML file for Unstructured Grid data
  //============================================
  template <typename S, typename T>  // S is cell-type, T is variable-type
  class VTU_file : public File {
    private:
    unsigned int offset_ = 0;
    Variables<T> variables_;
    Grid<S>& grid_; 
    //Data<node_dtype>& point_data_;
    //std::vector<Data<int>>& cell_data_;
    
    public:
    //--------------------------------------------
    //VTU_file(): File(), variables_(), point_data_(Data<node_dtype>()), cell_data_(std::vector<Data<int>>()) { } 
    //--------------------------------------------
    //--------------------------------------------
    // VTU_file(const std::string &_path, const std::string &_name, Data<node_dtype> &point_data, std::vector<Data<int>> &cell_data) 
    //   : File(_name, {_path, "vtu/"}, ".vtu"), variables_(), point_data_(point_data), cell_data_(cell_data) 
    VTU_file(const std::string &_path, const std::string &_name, Grid<S>& grid) 
      : File(_name, {_path, "vtu/"}, ".vtu"), variables_(), grid_(grid) { offset_ = grid.offset(); }
    //--------------------------------------------
    //{ 
      // offset_ = point_data_.update_offset(offset_);
      // for (auto& cell : cell_data_) 
      //   offset_ = cell.update_offset(offset_);
    //}
    //--------------------------------------------
    void write(const int rank)
    //--------------------------------------------
    {
      set_filename_and_open(rank);
      write_header();
      write_data();
      write_footer();
      write_appended_data(); 
      close();
      ++nwrite_;
    }
    // //--------------------------------------------
    // void add_variable() 
    // { 
      
    // }
    //--------------------------------------------
    void update_last_variable_offset() { 
    //--------------------------------------------
      offset_ = variables_.back().update_offset(offset_); 
    }
    //--------------------------------------------
    Variables<T>& variables() { return variables_; }
    //--------------------------------------------
    //--------------------------------------------
    Grid<S>& grid() { return grid_; }
    //--------------------------------------------
    //--------------------------------------------
    std::string endianess() {
    //--------------------------------------------
    // Function returning endianess for VTK-header 
    // Adapted from: https://stackoverflow.com/questions/4181951/how-to-check-whether-a-system-is-big-endian-or-little-endian/4181991
      int n = 1;
      if(*(char *)&n == 1) {
        return "LittleEndian";
      } else {
        return "BigEndian";
      }
    }
    //--------------------------------------------
    void write_header() {
    //--------------------------------------------
      file_ << "<?xml version=\"1.0\"?>" << std::endl;
      file_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianess() << "\">"   << std::endl;
      file_ << "  <UnstructuredGrid>" << std::endl;
      //file_ << "    <Piece NumberOfPoints=\"" << (int)(point_data_.size()/point_data_.dim_) << "\" NumberOfCells=\"" << cell_data_.back().size() << "\">" << std::endl;
      file_ << "    <Piece NumberOfPoints=\"" << grid_.num_points() << "\" NumberOfCells=\"" << grid_.num_cells() << "\">" << std::endl;
      file_ << "      <Points>" << std::endl;
      //point_data_.write_dataarray(file_);
      grid_.point_data().write_dataarray(file_);
      file_ << "      </Points>" << std::endl;
      file_ << "      <Cells>" << std::endl;
      //for (const auto& data : cell_data_)
      for (const auto& data : grid_.cell_data())
        data.write_dataarray(file_);
      file_ << "      </Cells>" << std::endl;
    };
    //--------------------------------------------
    void write_data() {
    //--------------------------------------------
      file_ << "      <CellData Scalars=\"" << variables_.scalar_string_ << "\" Vectors=\"" << variables_.vector_string_ << "\">" << std::endl;
      for (auto& data : variables_.data()) {
        data.write_dataarray(file_, 1e-20);
      }
      file_ << "      </CellData>" << std::endl;
    }
    //--------------------------------------------
    void write_appended_data() {
    //--------------------------------------------
      if (offset_ > 0) {
        file_ << "  <AppendedData encoding=\"raw\">" << std::endl;
        file_ << "_";
        grid_.point_data().write_binarydata(file_);
        for (const auto& cell : grid_.cell_data()) 
          cell.write_binarydata(file_);
        for (auto& data : variables_.data())
          data.write_binarydata(file_);
        file_ << "  </AppendedData>" << std::endl;
      }
    }
    //--------------------------------------------
    void write_footer() {
    //--------------------------------------------
      file_ << "    </Piece>" << std::endl;
      file_ << "  </UnstructuredGrid>" << std::endl;
    }
    //--------------------------------------------
    void set_filename_and_open(const int rank) {
    //--------------------------------------------
      std::ostringstream ss;
      ss << std::setfill('0') << std::setw(4) << rank << "_" << name_ << "_" << std::setw(7) << nwrite_ << extension_;
      filename_ = ss.str();
      open();
    }
    //--------------------------------------------
    const std::string piece_string() const {
    //--------------------------------------------
      std::ostringstream oss;
      oss << "    <Piece Source=\"" << folders_.back() + filename_ << "\" />";
      return oss.str();
    }
  };


  //============================================
  //             P V T U F I L  E
  //                
  // VTK Parallel XML file for Unstructured Grid
  //============================================
  template <typename S, typename T>
  class PVTU_file : public File {
    private:
    long file_position = 0;

    public:
    // constructor
    //--------------------------------------------
    PVTU_file(const std::string &_path, const std::string &_name) : File(_name, {_path}, ".pvtu") { }
    //--------------------------------------------
    //--------------------------------------------
    void write(const double time, const int rank, const int max_rank, VTU_file<S,T>& vtu)
    //--------------------------------------------
    {
      set_filename();
      if (rank == 0) {
        open();
        write_header(time, vtu);
        set_position_and_close();
      }
      write_piece(vtu.piece_string());
      //MPI_write_piece(vtu.piece_string(), rank);
      if (rank == max_rank) {
        write_footer();
        close();
      }
      ++nwrite_;
    }
    //--------------------------------------------
    std::string timestring() {
    //--------------------------------------------
      time_t t = time(NULL);
      tm tm = *localtime(&t);
      std::ostringstream ss;
      ss << std::setfill('0') << std::setw(2) << tm.tm_mday << "." << std::setw(2) << tm.tm_mon+1 << "."
          << std::setw(2) << tm.tm_year+1900 << " " << std::setw(2) << tm.tm_hour <<  ":" << std::setw(2)
      << tm.tm_min << ":" << std::setw(2) << tm.tm_sec;
      return ss.str();
    };
    //--------------------------------------------
    void write_piece(const std::string& piece_string) {
    //--------------------------------------------
      char piece[101] = {0}, piece_format[20] = {0};
      sprintf(piece_format, "%%-%ds\n", (int)(sizeof(piece))-1);  // -1 due to \n
      sprintf(piece, piece_format, piece_string.c_str());
      FILE *file = fopen((path_+filename_).c_str(), "a");
      if (file==nullptr) {
        std::cerr << "ERROR! Unable to open " << path_+filename_ << std::endl;
        exit(-1);
      }
      fseek(file, file_position, SEEK_SET);
      fwrite(piece, sizeof(char), sizeof(piece), file);
      fclose(file);
    }
    //--------------------------------------------
    void MPI_write_piece(const std::string& piece_string, const int rank) {
    //--------------------------------------------
      char piece[101] = {0}, piece_format[20] = {0};
      sprintf(piece_format, "%%-%ds\n", (int)(sizeof(piece))-2);
      sprintf(piece, piece_format, piece_string.c_str());
      MPI_File mpi_file;
      MPI_Status status;
      MPI_Bcast(&file_position, 1, MPI_LONG, 0, MPI_COMM_WORLD);
      MPI_Offset seek_position = (long long) (file_position + rank*(sizeof(piece)-1));
      int err = MPI_File_open(MPI_COMM_WORLD, (path_+filename_).c_str(), MPI_MODE_APPEND|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_file);
      if (err) {
        std::cerr << "ERROR! Unable to open " << path_+filename_ << std::endl;
      }
      MPI_File_seek(mpi_file, seek_position, MPI_SEEK_SET);
      MPI_File_write(mpi_file, piece, sizeof(piece)-1, MPI_CHAR, &status);
      MPI_File_close(&mpi_file);
    }

    //--------------------------------------------
    void write_header(double time, VTU_file<S,T>& vtu) {
    //--------------------------------------------
      //file_.precision(precision_);
      file_ << "<?xml version=\"1.0\"?>" << std::endl;
      file_ << "<!-- Created " << timestring() << " -->" << std::endl;
      file_ << "<!-- time = " << time << " s -->" << std::endl;
      file_ << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      file_ << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
      file_ << "    <PPoints>" << std::endl;
      file_ << "      <PDataArray type=\"Float32\" NumberOfComponents=\"" << vtu.grid().point_data().dim() << "\" />" << std::endl;
      file_ << "    </PPoints>" << std::endl;
      file_ << "    <PCellData>" << std::endl;
      file_ << "      <PDataArray type=\"Int32\" Name=\"connectivity\" />" << std::endl;
      file_ << "      <PDataArray type=\"Int32\" Name=\"offsets\" />" << std::endl;
      file_ << "      <PDataArray type=\"Int32\" Name=\"types\" />" << std::endl;
      file_ << "    </PCellData>" << std::endl;
      file_ << "    <PCellData>" << std::endl;
      for (auto& var : vtu.variables().data())
        file_ << "      <PDataArray " << var.dataarray("name") << var.dataarray("type") << var.dataarray("dim") << "/>" << std::endl;
      file_ << "    </PCellData>" << std::endl;
    }
    //--------------------------------------------
    void write_footer() {
    //--------------------------------------------
      file_.open(path_+filename_, std::ios::app);
      file_ << "  </PUnstructuredGrid>" << std::endl;
    }
    //--------------------------------------------
    void set_filename() {
    //--------------------------------------------
      std::ostringstream ss;
      ss << std::setfill('0') << name_ << "_" << std::setw(7) << nwrite_ << extension_;
      filename_ = ss.str();
    }
    //--------------------------------------------
    void set_position_and_close() {
    //--------------------------------------------
      file_position = file_.tellp();
      file_.close();
    }
  };

  //============================================
  //             O U T F I L E
  //============================================

  template <typename S, typename T>
  class Outfile {
    private:
    PVTU_file<S,T> pvtu_file_;
    VTU_file<S,T> vtu_file_;

    public:
    //--------------------------------------------
    //Outfile(std::string &_path, const std::string &_name, Data<node_dtype> &point_data, std::vector<Data<int>> &cell_data) 
    Outfile(std::string &_path, const std::string &_name, Grid<S>& grid) 
      : pvtu_file_(PVTU_file<S,T>(_path, _name)), vtu_file_(VTU_file<S,T>(_path, _name, grid)) { }
    //--------------------------------------------
    //--------------------------------------------
    void write(const double time, const int rank, const int max_rank)
    //--------------------------------------------
    {
      vtu_file_.write(rank);
      pvtu_file_.write(time, rank, max_rank, vtu_file_);
    }
    //--------------------------------------------
    VTU_file<S,T>& vtu_file() { return(vtu_file_); }
    //--------------------------------------------
    //--------------------------------------------
    PVTU_file<S,T>& pvtu_file() { return(pvtu_file_); }
    //--------------------------------------------

  };


  //============================================
  //               O U T P U T
  //============================================
  template <typename S, typename T>
  class Output {
    private:
    Grid<S> grid_;
    std::string path_;
    std::vector<Outfile<S,T>> outfiles_;
    std::unordered_map<std::string, int> get_index_;
    int nwrite_ = 0;
    int rank_ = 0;
    int max_rank_ = 0;
    int format_ = BINARY;

    public:
    // Constructors
    //--------------------------------------------
    // Output(const std::vector<VTK::point_dtype>& points, const std::vector<int>& connectivity, const std::string path, const int rank, const int num_procs, const int format) 
    //   : grid_(points, connectivity, format), path_(path), outfiles_(), get_index_(), rank_(rank), max_rank_(num_procs-1), format_(format) 
    // { }
    //--------------------------------------------
    Output(std::vector<std::vector<node_dtype>>& nodes, const std::string path, const int rank, const int num_procs, const int format) 
      : grid_(nodes, format), path_(path), outfiles_(), get_index_(), rank_(rank), max_rank_(num_procs-1), format_(format) { }
    //--------------------------------------------
    //--------------------------------------------
    Outfile<S,T>& operator[](const std::string& name) { return outfiles_[get_index_[name]]; } 
    //--------------------------------------------
    //--------------------------------------------
    Outfile<S,T>& operator[](const int index) { return outfiles_[index]; }
    //--------------------------------------------
    //--------------------------------------------
    Outfile<S,T>& file(const int index) { return outfiles_[index]; }
    //--------------------------------------------
    //--------------------------------------------
    inline void write(const std::string& var_name, const double time) { write(get_index_[var_name], time); }
    //--------------------------------------------
    //--------------------------------------------
    inline void write(const int index, const double time) {
    //--------------------------------------------
      write(outfiles_[index], time);
    }
    //--------------------------------------------
    const std::string& get_filename(const std::string& var_name) {
    //--------------------------------------------
      return(outfiles_[get_index_[var_name]].pvtu_file().get_filename());
    }
    //--------------------------------------------
    Outfile<S,T>& add_file(const std::string& name) {
    //--------------------------------------------
      //outfiles_.emplace_back(path_, name, grid_.point_data(), grid_.cell_data());
      outfiles_.emplace_back(path_, name, grid_);
      get_index_[name] = outfiles_.size()-1;
      return outfiles_.back();
    }    
    //--------------------------------------------
    Outfile<S,T>& last_outfile() { return outfiles_.back(); }
    //--------------------------------------------
    //--------------------------------------------
    void add_variable(const std::string &name, const int dim, std::vector<T>& data, std::vector<int>& index, int offset=0, const int length=0)
    //--------------------------------------------
    {
      last_outfile().vtu_file().variables().add(name, data, format_, dim, index, offset, length);
      last_outfile().vtu_file().update_last_variable_offset();
    }
    //--------------------------------------------
    void add_variable(const std::string &name, const int dim, std::vector<T>& data, int offset=0, const int length=0)
    //--------------------------------------------
    {
      std::vector<int> index;
      add_variable(name, dim, data, index, offset, length);
    }
    //--------------------------------------------
    void write(const double time) 
    //--------------------------------------------
    {
      for (auto& ofile : outfiles_) {
        write(ofile, time);    
      }
    }
    //--------------------------------------------
    void write(Outfile<S,T>& outfile, const double time)
    //--------------------------------------------
    { 
      outfile.write(time, rank_, max_rank_);
      ++nwrite_;
      // VTU_file<T>& vtu = outfile.vtu_file();
      // PVTU_file<T>& pvtu = outfile.pvtu_file();
      // vtu.set_filename_and_open(rank_);
      // vtu.write_header();
      // vtu.write_data();  
      // vtu.write_footer();
      // vtu.write_appended_data(); 
      // vtu.close();
      // pvtu.set_filename();
      // if (rank_ == 0) {
      //   pvtu.open();
      //   pvtu.write_header(time, vtu);
      //   pvtu.set_position_and_close();
      // }
      // pvtu.write_piece(vtu.piece_string());
      // //pvtu.MPI_write_piece(vtu.piece_string(), rank_);
      // if (rank_ == max_rank_) {
      //   pvtu.write_footer();
      //   pvtu.close();
      // }
      // ++(vtu.nwrite_);
      // ++(pvtu.nwrite_);
      // ++nwrite_;
    }
  };


}


#endif /* SRC_VTK_H_ */
