
#ifndef SRC_VTK_H_
#define SRC_VTK_H_
#include <unordered_map>
#include <map>
#include <deque>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <bitset>
#include <algorithm>
//#include <memory>
#include <mpi.h>
//#include <sys/types.h>
#include <sys/stat.h>
//#if defined (_WIN32)
//#include <direct.h>
//#else
#include <unistd.h>
//#endif
//#include "vector_func.h"

namespace util {

  template <typename T>
  //-----------------------------------------------------------------------------------  
  inline std::vector<T> linspace(T start, T stop, T inc=1) {
  //-----------------------------------------------------------------------------------  
    std::vector<T> vec((int)((stop-start)/inc), inc);
    vec[0] = start;
    std::partial_sum(vec.begin(), vec.end(), vec.begin(), std::plus<T>());
    return vec;
  }

  template <typename T>
  //-----------------------------------------------------------------------------------
  std::vector<T> add_coord(const std::vector<T>& vec, T pad_value, const int old_dim, const int new_dim) 
  //-----------------------------------------------------------------------------------
  {
    int N = int(vec.size()/old_dim);
    std::vector<T> vec_pad(new_dim * N, pad_value);
    for (int n=0; n<N; ++n) {
      for (int i=0; i<old_dim; ++i) {
        vec_pad[n*new_dim + i] = vec[n*old_dim + i];
      }
    }
    return vec_pad;
  }

  template <typename T>
  //-----------------------------------------------------------------------------------
  void print_vector(const std::vector<T>& vec, const int dim) 
  //-----------------------------------------------------------------------------------
  { 
    std::cout << "(";
    int n = vec.size();
    for (int i=0; i<n; ++i) {
        std::cout << vec[i] << (( i<n-1 && (i+1)%dim==0) ? "),(" : ",");
    }
    std::cout << ")" << std::endl;
  }

  template <typename T>
  //-----------------------------------------------------------------------------------
  const std::string to_string(const T value) { return std::to_string(value); }
  //-----------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------
  const std::string to_string(const std::string value) { return value; }
  //-----------------------------------------------------------------------------------

  template <typename T>
  //-----------------------------------------------------------------------------------
  const std::string vector_to_string(const std::vector<T>& vector, const std::string& separator = " ") 
  //-----------------------------------------------------------------------------------
  { 
    std::string str("");
    std::string sep("");
    for (const auto& elm : vector) {
        str += sep + to_string(elm);
        sep = separator;
    }
    return str;
  }

  //-----------------------------------------------------------------------------------
  void safe_exit(int exit_value) 
  //-----------------------------------------------------------------------------------
  {
    int mpi_running = 0;
    MPI_Initialized(&mpi_running);
    if (mpi_running) {
      MPI_Finalize();
    } 
    std::exit(exit_value);
  }

}

namespace VTK {
  
  using point_dtype = float;
  using nodes_vec = std::vector<std::vector<point_dtype>>;

  static constexpr int BINARY = 0;
  static constexpr int ASCII = 1;

  //-----------------------------------------------------------------------------------  
  //  Get datatype name
  //  Adapted from: https://stackoverflow.com/questions/4484982/how-to-convert-typename-t-to-string-in-c
  //-----------------------------------------------------------------------------------
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


  //--------------------------------------------------------------------------------------------------------------------------
  // 1: VTK_VERTEX
  //            o
  //-----------------------------------------------------------------------------------
  struct vertex {
    static constexpr char name[] = "Vertex";
    static constexpr int type = 1;
    static constexpr int dim = 1;
    static constexpr int n = 1;
    static constexpr std::array<std::array<int, dim>, n> points = { {{0}} };
  };
  // These three lines are necessary to avoid linker errors in c++11, but can be removed for c++17
  constexpr std::array<std::array<int, vertex::dim>, vertex::n> vertex::points;
  constexpr int vertex::type;
  constexpr char vertex::name[];
  
  //-----------------------------------------------------------------------------------
  // 3: VTK_LINE
  //             o
  //            /
  //           /
  //          / 
  //         o
  //----------------------------------------------------------------------------------- 
  struct line {
    static constexpr char name[] = "Line";
    static constexpr int type = 3;
    static constexpr int dim = 1;
    static constexpr int n = 2;
    static constexpr std::array<std::array<int, dim>, n> points = { {{0}, {1}} };
  };
  // These three lines are necessary to avoid linker errors in c++11, but can be removed in c++17
  constexpr std::array<std::array<int, line::dim>, line::n> line::points;
  constexpr int line::type;
  constexpr char line::name[];

  //----------------------------------------------------------------------------------- 
  //  8: VTK_PIXEL
  //
  //             2 o--------o 3
  //               |        | 
  //               |        |         y
  //               |        |         | 
  //             0 o--------o 1       |____ x
  //                      
  //----------------------------------------------------------------------------------- 
  struct pixel {
    static constexpr char name[] = "Pixel";
    static constexpr int type = 8;
    static constexpr int dim = 2;
    static constexpr int n = 4;
    static constexpr std::array<std::array<int, dim>, n> points = {{{0,0}, {1,0}, {0,1}, {1,1}}};
  };
  // These three lines are necessary to avoid linker errors in c++11, but can be removed in c++17
  constexpr std::array<std::array<int, pixel::dim>, pixel::n> pixel::points;
  constexpr int pixel::type;
  constexpr char pixel::name[];

  /*----------------------------------------------------------------------------------- 
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
  //----------------------------------------------------------------------------------- 
  */
  struct quad {
    static constexpr char name[] = "Quad";
    static constexpr int type = 9;
    static constexpr int dim = 2;
    static constexpr int n = 4;
    static constexpr std::array<std::array<int, dim>, n> points = {{{0,0}, {1,0}, {1,1}, {0,1}}};
  };
  // These three lines are necessary to avoid linker errors in c++11, but can be removed for c++17
  constexpr std::array<std::array<int, quad::dim>, quad::n> quad::points;
  constexpr int quad::type;
  constexpr char quad::name[];

  //-----------------------------------------------------------------------------------
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
  //-----------------------------------------------------------------------------------
  struct voxel {
    static constexpr char name[] = "Voxel";
    static constexpr int type = 11;
    static constexpr int dim = 3;
    static constexpr int n = 8;
    static constexpr std::array<std::array<int, dim>, n> points = {{{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}}};
  };
  // These three lines are necessary to avoid linker errors in c++11, but can be removed for c++17
  constexpr std::array<std::array<int, voxel::dim>, voxel::n> voxel::points;
  constexpr int voxel::type;
  constexpr char voxel::name[];

  //=====================================================================================
  //
  //                              D A T A _ W R A P P E R
  //
  //=====================================================================================
  template <typename T>
  class data_wrapper
  {
      public:
      data_wrapper() { } //{ std::cout << "data_wrapper constructor" << std::endl; }
      virtual ~data_wrapper() { };
      virtual const T at(const int pos) const = 0;
      virtual const T* ptr(const int pos) const = 0;
      virtual const T* begin() const = 0;
      virtual const T* end() const = 0;
      virtual const size_t size() const = 0;
  };

  template <typename T>
  class vec_wrapper : public data_wrapper<T>
  {
      private:
          const std::vector<T>& data_;
      public:
          vec_wrapper(const std::vector<T>& data) : data_(data) { } //{ std::cout << "vec_wrapper constructor" << std::endl; }
          const T at(const int pos) const { return data_[pos]; }
          const T* ptr(const int pos) const { return &data_[pos]; }
          const T* begin() const { return &data_[0]; }
          const T* end() const { return &data_[data_.size()]; }
          const size_t size() const { return data_.size(); }
  };

  template <typename T>
  class arr_wrapper : public data_wrapper<T>
  {
      private:
          const std::valarray<T>& data_;
      public:
          arr_wrapper(const std::valarray<T>& data) : data_(data) { } //{  std::cout << "arr_wrapper constructor" << std::endl;}
          const T at(const int pos) const { return data_[pos]; }
          const T* ptr(int pos) const { return &data_[pos]; }
          const T* begin() const { return &data_[0]; }
          const T* end() const { return &data_[data_.size()]; }
          const size_t size() const { return data_.size(); }
  };


  //=====================================================================================
  //
  //                                    M E S H
  //
  //=====================================================================================
  template <typename T>
  struct Mesh_data
  {
    std::vector<T> vec_;
    vec_wrapper<T> data_;
    Mesh_data() : vec_(), data_(vec_) { }
    Mesh_data(const std::vector<T>& v) : vec_(v), data_(vec_) { }
  };

  template <typename CELL>
  class Mesh
  {
  private:
    //std::vector<point_dtype> points_;
    Mesh_data<point_dtype> points_;
    //std::vector<int> conn_;
    Mesh_data<int> conn_;
    //std::vector<int> offsets_;
    Mesh_data<int> offsets_;
    //std::vector<int> types_;
    Mesh_data<int> types_;
    std::vector<int> size_;
    point_dtype unit_ = 1;
    static const int point_dim_ = 3;  // VTK seems to only accept 3D point data

  public:

    //                                   Mesh
    //-----------------------------------------------------------------------------------
    Mesh(std::vector<point_dtype> &points, const std::vector<int> &conn, point_dtype unit = 1) 
    //-----------------------------------------------------------------------------------
      : points_(points), conn_(conn), offsets_(), types_(), size_(), unit_(unit)
    {
      init();
    } 

    //                                   Mesh
    //-----------------------------------------------------------------------------------
    template <typename T>
    Mesh(const std::vector<std::vector<T>> &nodes, point_dtype unit = 1) 
    //-----------------------------------------------------------------------------------
      : points_(), conn_(), offsets_(), types_(), size_(), unit_(unit)
    {
      set_size(nodes);
      calc_points_and_conn(nodes);
      init();
    }
    //const std::vector<point_dtype>& points_vec() const { return points_.vec_; }
    const vec_wrapper<point_dtype>& points() const { return points_.data_; }
    //const std::vector<int>& connectivity_vec() const { return conn_.vec_; }
    const vec_wrapper<int>& connectivity() const { return conn_.data_; }
    //const std::vector<int>& offsets_vec() const { return offsets_.vec_; }
    const vec_wrapper<int>& offsets() const { return offsets_.data_; }
    //const std::vector<int>& types_vec() const { return types_.vec_; }
    const vec_wrapper<int>& types() const { return types_.data_; }
    int dim() const { return point_dim_; }
    int num() const { return int(conn_.vec_.size()/CELL::n); }

  private:
    //                                   Mesh
    //-----------------------------------------------------------------------------------
    void init() 
    //-----------------------------------------------------------------------------------
    {
      offsets_.vec_ = util::linspace<int>(CELL::n, conn_.vec_.size()+CELL::n, CELL::n);
      types_.vec_ = std::vector<int>(num(), CELL::type);
      // std::cout << "type, dim, n: " << CELL::type << ", " << CELL::dim << ", " << CELL::n << std::endl;
      // std::cout << "points: "; util::print_vector(points_, point_dim_);
      // std::cout << "connectivity: "; util::print_vector(conn_, CELL::n);
      // std::cout << "offsets: " << offsets_ << std::endl;
      // std::cout << "types: " << types_ << std::endl;
    }

    //                                   Mesh
    //-----------------------------------------------------------------------------------
    template <typename T>
    void set_size(const std::vector<std::vector<T>> &nodes)
    //-----------------------------------------------------------------------------------
    {
      std::vector<T> min_n(CELL::dim, 999999), max_n(CELL::dim, -1);
      for (const auto &node : nodes) {
        for (auto i=0; i < CELL::dim; i++) {
          if (node[i] > max_n[i])
            max_n[i] = node[i];
          if (node[i] < min_n[i])
            min_n[i] = node[i];
        }
      }
      std::vector<int> size(CELL::dim, 0); 
      for (auto i = 0; i < CELL::dim; ++i) {
        size[i] = int(max_n[i] - min_n[i]) + 1;
      }
      size_ = size;
      //std::cout << size_ << std::endl;
    }

    //                                   Mesh
    //-----------------------------------------------------------------------------------
    template <typename T>
    std::vector<int> get_stride(const std::vector<std::vector<T>> &nodes) const 
    //-----------------------------------------------------------------------------------
    {
      std::vector<int> stride_vec(CELL::dim, 1); 
      for (auto i = 1; i < CELL::dim; ++i) {
        stride_vec[i] = stride_vec[i-1] * (size_[i-1] + 1);
      }
      return stride_vec;
    }

    //                                   Mesh
    //-----------------------------------------------------------------------------------
    template <typename T>
    void calc_points_and_conn(const std::vector<std::vector<T>> &nodes) 
    //-----------------------------------------------------------------------------------
    {
      if (nodes[0].size() != CELL::dim) {
        std::cerr << "ERROR in VTK::Mesh: A " << CELL::name << " is " << CELL::dim << "-dimensional, but the given nodes are " << nodes[0].size() << "-dimensional" << std::endl;
        
        //std::exit(EXIT_FAILURE);
        util::safe_exit(EXIT_FAILURE);
      }
      // Loop over nodes and calculate corner points and unique indexes 
      auto stride = get_stride(nodes);
      std::vector<point_dtype> pts;
      std::vector<int> index;
      pts.reserve(CELL::n * nodes.size() * CELL::dim); 
      index.reserve(CELL::n * nodes.size());
      for (const auto &node : nodes) {
        for (const auto &cell_point : CELL::points) {
          int idx = 0;
          for (auto i=0; i < CELL::dim; i++) {
            auto p = node[i] + cell_point[i];
            idx += p*stride[i];
            pts.push_back(p);
          }
          index.push_back(idx);
        }
      }
      // Create unique index_list 
      auto minmax = std::minmax_element(index.begin(), index.end());
      auto min = *minmax.first;
      auto max = *minmax.second;
      std::vector<int> index_list(max - min + 1, -1);
      conn_.vec_.reserve(CELL::n * nodes.size());
      int n_pts = 0;
      auto pos = pts.begin(); 
      for (size_t i=0; i<index.size(); ++i) {
        int idx = index[i] - min;
        if (index_list[idx] < 0) {
          points_.vec_.insert(points_.vec_.end(), pos, pos + CELL::dim);
          index_list[idx] = n_pts;
          ++n_pts;
        }
        pos += CELL::dim;
        conn_.vec_.push_back(index_list[idx]);
      }
      // Shift to make nodes the cell centres 
      for (auto& p : points_.vec_)
        p -= 0.5*unit_;
      // Pad with 0 if dim is less than 3
      if (CELL::dim < point_dim_) {
        points_.vec_ = util::add_coord(points_.vec_,  (point_dtype)0, CELL::dim, point_dim_);
      }
    }
  };

  //=====================================================================================
  //
  //                                  D A T A A R R A Y
  //
  //=====================================================================================
  class DataArray {
    private:
    std::vector<std::string> data_formats_ = {"appended", "ascii"};
    std::map<std::string, std::string> tag_var_  = 
    { 
      {"name", "Name"}, {"type", "type"}, {"dim", "NumberOfComponents"}, {"format", "format"}, {"offset", "offset"}
    };
    std::map<std::string, std::string> tags_;

    public:
    //                                DataArray
    //-----------------------------------------------------------------------------------
    DataArray() : tags_() { }
    //-----------------------------------------------------------------------------------

    //                                DataArray
    //-----------------------------------------------------------------------------------
    DataArray(const std::string& name, int dim, const int format=0, const std::string& type="", int offset=-1)
    //-----------------------------------------------------------------------------------
     : tags_() 
    {
      set_tag("name", name);
      set_tag("type", type);
      set_tag("dim", dim);
      set_tag("format", data_formats_[format]);
      set_tag("offset", (format==ASCII) ? -1 : offset );
    }

    //                                DataArray
    //-----------------------------------------------------------------------------------
    std::string quote(const std::string& str) const { return "=\""+str+"\" "; }
    //-----------------------------------------------------------------------------------

    //                                DataArray
    //-----------------------------------------------------------------------------------
    void set_tag(const std::string& tname, int val)
    //-----------------------------------------------------------------------------------
    {
      tags_[tname] = (val>=0) ? tag_var_[tname]+quote(std::to_string(val)) : ""; 
    }

    //                                DataArray
    //-----------------------------------------------------------------------------------
    void set_tag(const std::string& tname, const std::string& val) 
    //-----------------------------------------------------------------------------------
    { 
      tags_[tname] = (val.empty()) ? "" : tag_var_[tname]+quote(val); 
    }

    //                                DataArray
    //-----------------------------------------------------------------------------------
    const std::string& operator[](const std::string& tname) const { return tags_.at(tname); }
    //-----------------------------------------------------------------------------------
 
    //                                DataArray
    //-----------------------------------------------------------------------------------
    friend std::ostream& operator<<(std::ostream& out, const DataArray& da)
    //-----------------------------------------------------------------------------------
    {
      out << da["name"] << da["type"] << da["dim"] << da["format"] << da["offset"];
      return out;
    }
  };

  //=====================================================================================
  //
  //                                     D A T A
  //
  //=====================================================================================
  template <typename T>
  class Data {
    public:
    std::string name_ = "";
    //const std::vector<T>& data_;
    const data_wrapper<T>& data_;
    int dim_ = 0;
    unsigned int offset_ = 0;
    unsigned int length_ = 0;
    DataArray dataarray_;
    unsigned int nbytes_ = 0;
    std::vector<int> index_;
    int contiguous_ = 1; 

    public:
    //                                   Data
    //-----------------------------------------------------------------------------------
    Data() : data_(std::vector<T>()), dataarray_(), index_() {};                                  
    //-----------------------------------------------------------------------------------

    //                                   Data
    //-----------------------------------------------------------------------------------
    //Data(const std::string &name, const std::vector<T>& data, const int format, const int dim, const std::vector<int>& index=std::vector<int>(), const int length=0, const int offset=0) 
    Data(const std::string &name, const data_wrapper<T>& data, const int format, const int dim, const std::vector<int>& index=std::vector<int>(), const int length=0, const int offset=0) 
    //-----------------------------------------------------------------------------------
      : name_(name), data_(data), dim_(dim), offset_(offset), length_(length), dataarray_(name, dim, format), index_(index)  
    { 
      dataarray_.set_tag("type", datatype<T>::name());
      if (index.empty()) {
        // Assume contiguous data, create default index-vector 
        contiguous_ = 1;
        set_length_nbytes(data_.size(), format);
        auto ind = util::linspace(offset_, offset_+length_);
        index_.insert(index_.begin(), ind.begin(), ind.end());
      } else {
        // Index-vector is provided, data is non-contiguous
        contiguous_ = 0;
        set_length_nbytes(index_.size(), format);
      }
      std::cout << dataarray_ << ", data-size: " << data_.size() << ", index-size: " << index_.size() << std::endl; 
    }  

    //                                   Data
    //-----------------------------------------------------------------------------------
    void set_length_nbytes(const int size, const int format)
    //-----------------------------------------------------------------------------------
    {
      if (length_ == 0) {
        length_ = size;
      }
      if (format == BINARY) {
        nbytes_ = length_*sizeof(T);
      }
    }
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    size_t size() const { return data_.size(); } 
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    const std::string& name() const { return name_; } 
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    int dim() const { return dim_; } 
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    bool is_binary() const {return (nbytes_ > 0) ? true : false; }
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    void write_asciidata (std::ofstream& file, double min=-1) const
    //-----------------------------------------------------------------------------------
    {
      if (!is_binary()) {
        for (const auto& ind : index_) {
          // T val = data_[ind];
          T val = data_.at(ind);
          if (min>0 && std::abs(val)<min)
            val = 0.0; 
          file << " " << val;
        }        
      }
    }

    //                                   Data
    //-----------------------------------------------------------------------------------
    void write_binarydata(std::ofstream& file) const {
    //-----------------------------------------------------------------------------------
      if (is_binary()) {
        file.write((char*)&nbytes_, sizeof(unsigned int));
        if (contiguous_) {
          // file.write((char*)&data_[offset_], nbytes_);
          file.write((char*)data_.ptr(offset_), nbytes_);
        } else {
          for (const auto& ind : index_) {
            // file.write((char*)&data_[ind], sizeof(T));            
            file.write((char*)data_.ptr(ind), sizeof(T));            
          }
        }
      }
    }
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    const std::string& dataarray(const std::string& tname) const { return dataarray_[tname]; }
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    int update_offset(int offset) { 
    //-----------------------------------------------------------------------------------
      if (is_binary() && dataarray_["offset"].empty()) {
        dataarray_.set_tag("offset", offset); 
        return offset + nbytes_ + sizeof(unsigned int);
      } else {
        return offset;
      }
    }
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    void write_dataarray(std::ofstream& file, double min=-1) const {
    //-----------------------------------------------------------------------------------
      file << "<DataArray " << dataarray_ << ">";
      write_asciidata(file, min);  
      file << "</DataArray>" << std::endl;
    }
  };

  //=====================================================================================
  //
  //                                     G R I D
  //
  //=====================================================================================
  template <typename CELL>
  class Grid {
    private:
    Mesh<CELL> mesh_;
    Data<point_dtype> point_data_;
    std::vector<Data<int>> cell_data_;
    unsigned int offset_ = 0;

    public:
    //                                   Grid
    //-----------------------------------------------------------------------------------
    Grid() : mesh_(), point_data_(), cell_data_() { }
    //-----------------------------------------------------------------------------------

    //                                   Grid
    //-----------------------------------------------------------------------------------
    template <typename S>
    Grid(const std::vector<std::vector<S>>& nodes, int format=BINARY) 
    : mesh_(nodes), point_data_("points", mesh_.points(), format, mesh_.dim()), cell_data_()
    //-----------------------------------------------------------------------------------
    {
      cell_data_.emplace_back("connectivity", mesh_.connectivity(), format, 1);
      cell_data_.emplace_back("offsets",      mesh_.offsets()     , format, 1);
      cell_data_.emplace_back("types",        mesh_.types()       , format, 1);
      update_offset();
    }

    //                                   Grid
    //-----------------------------------------------------------------------------------
    void update_offset()
    //-----------------------------------------------------------------------------------
    { 
      offset_ = point_data_.update_offset(offset_);
      for (auto& cell : cell_data_) 
        offset_ = cell.update_offset(offset_);
    }

    //                                   Grid
    //-----------------------------------------------------------------------------------
    int offset() const { return offset_; }
    //-----------------------------------------------------------------------------------

    //                                   Grid
    //-----------------------------------------------------------------------------------
    int num_points() const { return int(point_data_.size()/point_data_.dim()); }
    //-----------------------------------------------------------------------------------

    //                                   Grid
    //-----------------------------------------------------------------------------------
    size_t num_cells() const { return cell_data_.back().size(); };
    //-----------------------------------------------------------------------------------

    //                                   Grid
    //-----------------------------------------------------------------------------------
    const Data<point_dtype>& point_data() const { return point_data_; };
    //-----------------------------------------------------------------------------------

    //                                   Grid
    //-----------------------------------------------------------------------------------
    const std::vector<Data<int>>& cell_data() const { return cell_data_; };
    //-----------------------------------------------------------------------------------
  };


  //=====================================================================================
  // 
  //                                V A R I A B L E S
  //
  //=====================================================================================
  template <typename T>
  class Variables {
    public: 
    std::vector<Data<T>> datalist_;
    std::string cell_data_string_ = "";
    std::string scalar_names_= "", vector_names_= "";
    std::vector<std::string> scalars_, vectors_;
    
    //                                   Variables
    //-----------------------------------------------------------------------------------
    Variables() : datalist_(), scalars_(), vectors_() { }    
    //-----------------------------------------------------------------------------------

    //                                   Variables
    //-----------------------------------------------------------------------------------
    //void add(const std::string& name, const std::vector<T>& data, const int format, const int dim, const std::vector<int>& index, const int length=0, const int offset=0) 
    void add(const std::string& name, const data_wrapper<T>& data, const int format, const int dim, const std::vector<int>& index, const int length=0, const int offset=0) 
    //-----------------------------------------------------------------------------------
    {
      datalist_.emplace_back(name, data, format, dim, index, length, offset);      
      update_names(datalist_.back());
    }

    //                                   Variables
    //-----------------------------------------------------------------------------------
    const std::vector<Data<T>>& data() const { return datalist_; }
    //-----------------------------------------------------------------------------------

    //                                   Variables
    //-----------------------------------------------------------------------------------
    const Data<T>& back() const { return datalist_.back(); }
    Data<T>& back() { return datalist_.back(); }
    //-----------------------------------------------------------------------------------

    //                                   Variables
    //-----------------------------------------------------------------------------------
    const std::string& vector_names() const { return vector_names_; }
    //-----------------------------------------------------------------------------------

    //                                   Variables
    //-----------------------------------------------------------------------------------
    const std::string& scalar_names() const { return scalar_names_; }
    //-----------------------------------------------------------------------------------

    private:
    //                                   Variables
    //-----------------------------------------------------------------------------------
    void update_names(const Data<T>& var) {
    //-----------------------------------------------------------------------------------
      //const auto& var = datalist_.back();
      if (var.dim() > 1) {
        vectors_.push_back(var.name());
        vector_names_ = util::vector_to_string(vectors_, ", ");
      } else {
        scalars_.push_back(var.name());
        scalar_names_ = util::vector_to_string(scalars_, ", ");
      }
    }

  };

  //=====================================================================================
  //                                    F I L  E
  //                
  //                   Base class for the .vtu and .pvtu file classes
  //=====================================================================================
  class File {
    protected:
    std::string name_ = "";
    std::ofstream file_;
    std::string filename_ = "";
    std::vector<std::string> folders_;
    std::string path_ = "";
    std::string extension_ = "";
    const int precision_ = 5;
    int nwrite_ = 0;

    public:
    //                                    File
    //-----------------------------------------------------------------------------------
    File() : file_(), folders_() { }
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    File(const std::string &name, const std::vector<std::string> &folders, const std::string &extension)
    : name_(name), file_(), folders_(folders), extension_(extension)
    //-----------------------------------------------------------------------------------
    {
      for (const auto& f : folders_) {
        path_ += f;
        make_dir(path_);
      }
    }

    //                                    File
    //-----------------------------------------------------------------------------------
    void open() { file_.open(path_+filename_, std::ios::out); }
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    void close() {
    //-----------------------------------------------------------------------------------
      file_ << "</VTKFile>" << std::endl;
      file_.close();
    }

    //                                    File
    //-----------------------------------------------------------------------------------
    const std::string& filename() const { return filename_;}
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    const std::ofstream& file() const { return file_; }
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    const std::string path() const { return folders_.back()+filename_; }
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    void inc_nwrite() { ++nwrite_; }
    //-----------------------------------------------------------------------------------

    private:
    //                                    File
    //-----------------------------------------------------------------------------------
    void make_dir(std::string &dir) {
    //-----------------------------------------------------------------------------------
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


  //=====================================================================================
  //                                  V T U F I L E
  //                
  // VTK Serial XML file for Unstructured Grid data
  //=====================================================================================
  class VTU_file : public File {
    private:
    unsigned int offset_ = 0;
    
    public:
    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    VTU_file(const std::string &_path, const std::string &_name, unsigned int offset = 0) 
    //-----------------------------------------------------------------------------------
      : File(_name, {_path, "vtu/"}, ".vtu"), offset_(offset) { }

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    template <typename CELL, typename T>
    void write(const int rank, const Grid<CELL>& grid, Variables<T>& var) 
    //-----------------------------------------------------------------------------------
    {
      set_filename_and_open(rank);
      write_header(grid);
      write_data(var);
      write_footer();
      write_appended_data(grid, var); 
      close();
      inc_nwrite();
    }

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    void set_offset(unsigned int offset) { offset_ = offset; }
    //-----------------------------------------------------------------------------------

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    unsigned int offset() const { return offset_; }
    //-----------------------------------------------------------------------------------

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    std::string endianess() const 
    //-----------------------------------------------------------------------------------
    // Function returning endianess for VTK-header 
    // Adapted from: https://stackoverflow.com/questions/4181951/how-to-check-whether-a-system-is-big-endian-or-little-endian/4181991
    {
      int n = 1;
      if(*(char *)&n == 1) {
        return "LittleEndian";
      } else {
        return "BigEndian";
      }
    }

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    template <typename CELL>
    void write_header(const Grid<CELL>& grid) {
    //-----------------------------------------------------------------------------------
      file_ << "<?xml version=\"1.0\"?>" << std::endl;
      file_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianess() << "\">"   << std::endl;
      file_ << "  <UnstructuredGrid>" << std::endl;
      file_ << "    <Piece NumberOfPoints=\"" << grid.num_points() << "\" NumberOfCells=\"" << grid.num_cells() << "\">" << std::endl;
      file_ << "      <Points>" << std::endl;
      grid.point_data().write_dataarray(file_);
      file_ << "      </Points>" << std::endl;
      file_ << "      <Cells>" << std::endl;
      for (const auto& data : grid.cell_data())
        data.write_dataarray(file_);
      file_ << "      </Cells>" << std::endl;
    };

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    template <typename T>
    void write_data(const Variables<T>& var) {
    //-----------------------------------------------------------------------------------
      //file_ << "      <CellData Scalars=\"" << var.scalar_string_ << "\" Vectors=\"" << var.vector_string_ << "\">" << std::endl;
      file_ << "      <CellData Scalars=\"" << var.scalar_names() << "\" Vectors=\"" << var.vector_names() << "\">" << std::endl;
      for (auto& data : var.data()) {
        data.write_dataarray(file_, 1e-20);
      }
      file_ << "      </CellData>" << std::endl;
    }

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    template <typename CELL, typename T>
    void write_appended_data(const Grid<CELL>& grid, const Variables<T>& var) {
    //-----------------------------------------------------------------------------------
      if (offset_ > 0) {
        file_ << "  <AppendedData encoding=\"raw\">" << std::endl;
        file_ << "_";
        grid.point_data().write_binarydata(file_);
        for (const auto& cell : grid.cell_data()) {
          cell.write_binarydata(file_);
        }
        for (auto& data : var.data()) {
          data.write_binarydata(file_);
        }
        file_ << "  </AppendedData>" << std::endl;
      }
    }

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    void write_footer() {
    //-----------------------------------------------------------------------------------
      file_ << "    </Piece>" << std::endl;
      file_ << "  </UnstructuredGrid>" << std::endl;
    }

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    void set_filename_and_open(const int rank) {
    //-----------------------------------------------------------------------------------
      std::ostringstream ss;
      ss << std::setfill('0') << std::setw(4) << rank << "_" << name_ << "_" << std::setw(7) << nwrite_ << extension_;
      filename_ = ss.str();
      open();
    }
  };


  //=====================================================================================
  //
  //                                 P V T U F I L E
  //                
  //                  VTK Parallel XML file for Unstructured Grid
  //=====================================================================================
  class PVTU_file : public File {
    private:
    long file_position = 0;
    int mpi_running_ = 0;

    public:
    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    PVTU_file(const std::string &_path, const std::string &_name) : File(_name, {_path}, ".pvtu")
    //-----------------------------------------------------------------------------------
    { 
      // Check if this is a MPI-run
      MPI_Initialized(&mpi_running_); 
    }

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    template <typename CELL, typename T>
    void write(const double time, const int rank, const int max_rank, const Grid<CELL>& grid, const Variables<T>& var, const std::string& vtu_name)
    //-----------------------------------------------------------------------------------
    {
      set_filename();
      if (rank == 0) {
        open();
        write_header(time, grid, var);
        set_position_and_close();
      }
      if (mpi_running_)
        MPI_write_piece(vtu_name, rank);
      else
        write_piece(vtu_name);
      if (rank == max_rank) {
        write_footer();
        close();
      }
      inc_nwrite();
    }

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    std::string timestring() const {
    //-----------------------------------------------------------------------------------
      time_t t = time(NULL);
      tm tm = *localtime(&t);
      std::ostringstream ss;
      ss << std::setfill('0') << std::setw(2) << tm.tm_mday << "." << std::setw(2) << tm.tm_mon+1 << "."
          << std::setw(2) << tm.tm_year+1900 << " " << std::setw(2) << tm.tm_hour <<  ":" << std::setw(2)
      << tm.tm_min << ":" << std::setw(2) << tm.tm_sec;
      return ss.str();
    };

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    const std::string piece_string(const std::string& piece_file) const {
    //-----------------------------------------------------------------------------------
      std::ostringstream oss;
      oss << "    <Piece Source=\"" << piece_file << "\" />";
      return oss.str();
    }

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    void write_piece(const std::string& vtu_filename) const {
    //-----------------------------------------------------------------------------------
      char piece[101] = {0}, piece_format[20] = {0};
      sprintf(piece_format, "%%-%ds\n", (int)(sizeof(piece))-1);  // -1 due to \n
      sprintf(piece, piece_format, piece_string(vtu_filename).c_str());
      FILE *file = fopen((path_+filename_).c_str(), "a");
      if (file==nullptr) {
        std::cerr << "ERROR! Unable to open " << path_+filename_ << std::endl;
        exit(-1);
      }
      fseek(file, file_position, SEEK_SET);
      fwrite(piece, sizeof(char), sizeof(piece), file);
      fclose(file);
    }

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    void MPI_write_piece(const std::string& vtu_filename, const int rank) {
    //-----------------------------------------------------------------------------------
      char piece[101] = {0}, piece_format[20] = {0};
      sprintf(piece_format, "%%-%ds\n", (int)(sizeof(piece))-2);
      sprintf(piece, piece_format, piece_string(vtu_filename).c_str());
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

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    template <typename CELL, typename T>
    void write_header(double time, const Grid<CELL>& grid, const Variables<T>& var) {
    //-----------------------------------------------------------------------------------
      //file_.precision(precision_);
      file_ << "<?xml version=\"1.0\"?>" << std::endl;
      file_ << "<!-- Created " << timestring() << " -->" << std::endl;
      file_ << "<!-- time = " << time << " s -->" << std::endl;
      file_ << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      file_ << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
      file_ << "    <PPoints>" << std::endl;
      file_ << "      <PDataArray type=\"Float32\" NumberOfComponents=\"" << grid.point_data().dim() << "\" />" << std::endl;
      file_ << "    </PPoints>" << std::endl;
      file_ << "    <PCellData>" << std::endl;
      file_ << "      <PDataArray type=\"Int32\" Name=\"connectivity\" />" << std::endl;
      file_ << "      <PDataArray type=\"Int32\" Name=\"offsets\" />" << std::endl;
      file_ << "      <PDataArray type=\"Int32\" Name=\"types\" />" << std::endl;
      file_ << "    </PCellData>" << std::endl;
      file_ << "    <PCellData>" << std::endl;
      for (auto& v : var.data())
        file_ << "      <PDataArray " << v.dataarray("name") << v.dataarray("type") << v.dataarray("dim") << "/>" << std::endl;
      file_ << "    </PCellData>" << std::endl;
    }

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    void write_footer() {
    //-----------------------------------------------------------------------------------
      file_.open(path_+filename_, std::ios::app);
      file_ << "  </PUnstructuredGrid>" << std::endl;
    }

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    void set_filename() {
    //-----------------------------------------------------------------------------------
      std::ostringstream ss;
      ss << std::setfill('0') << name_ << "_" << std::setw(7) << nwrite_ << extension_;
      filename_ = ss.str();
    }

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    void set_position_and_close() {
    //-----------------------------------------------------------------------------------
      file_position = file_.tellp();
      file_.close();
    }
  };

  //=====================================================================================
  //  
  //                                  O U T F I L E
  //
  //=====================================================================================

  template <typename T>
  class Outfile {
    private:
    Variables<T> variables_;
    PVTU_file pvtu_file_;
    VTU_file vtu_file_;

    public:
    //                                  Outfile
    //-----------------------------------------------------------------------------------
    Outfile(std::string &_path, const std::string &_name, const unsigned int _offset) 
      : variables_(), pvtu_file_(_path, _name), vtu_file_(_path, _name, _offset) { }
    //-----------------------------------------------------------------------------------

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    template <typename CELL>
    void write(const Grid<CELL>& grid, double time, int rank, int max_rank)
    //-----------------------------------------------------------------------------------
    {
      vtu_file_.write(rank, grid, variables_);
      pvtu_file_.write(time, rank, max_rank, grid, variables_, vtu_file_.path());
    }

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    void update_offset() { 
    //-----------------------------------------------------------------------------------
      vtu_file_.set_offset( variables_.back().update_offset( vtu_file_.offset() ) ); 
    }

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    const VTU_file& vtu_file() const { return vtu_file_; }
    //-----------------------------------------------------------------------------------

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    const PVTU_file& pvtu_file() const { return pvtu_file_; }
    //-----------------------------------------------------------------------------------

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    const Variables<T>& variables() const { return variables_; }
    //-----------------------------------------------------------------------------------

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    Variables<T>& variables() { return variables_; }
    //-----------------------------------------------------------------------------------

  };

  //=====================================================================================
  //
  //                                    B U F F E R 
  //
  //=====================================================================================

  // template <typename T>
  // class base_data 
  // {
  //   // private:
  //   //   T val;
  //   public:
  //     base_data() {  std::cout << "base_data constructor" << std::endl; };
  //     //virtual ~base_data() { };
  //     virtual const T* start() const = 0; 
  //     virtual const T value(const int pos) const = 0;
  // };

  // template <typename T>
  // class vec_data : public base_data<T> 
  // {
  //   private:
  //     const std::vector<T>& data_;
  //   public:
  //   vec_data(const std::vector<T>& data) : data_(data) { std::cout << "vec_data constructor" << std::endl; }
  //   const T* start() const {  std::cout << "vec_data::start()" << std::endl; return &data_[0]; };
  //   const T value(const int pos) const = 0;
  // };

  // template <typename T>
  // class arr_data : public base_data<T> 
  // {
  //   private:
  //     const std::valarray<T>& data_;
  //   public:
  //   arr_data(const std::valarray<T>& data) : data_(data) {  std::cout << "arr_data constructor" << std::endl;}
  //   const T* start() const { std::cout << "arr_data::start()" << std::endl; return &data_[0]; };
  // };

  // template <typename T>
  // class Buffer 
  // {
  //   private:
  //   std::string name_;
  //   std::vector<T> buffer_;
  //   //const std::valarray<T>& data_;
  //   const std::valarray<T>& arr_data_;
  //   const std::vector<T>& vec_data_;
  //   int vec_buffer_;
  //   // const base_data<T>* data_;


  //   public:
  //   void info() {
  //     // std::cout << name_ << ": buf = " << buffer_.size() << ", arr = " << arr_data_.size() << ", vec = " << vec_data_.size() << ", base_data[0] = " << *(data_->start()) << std::endl;
  //     std::cout << name_ << ": buf = " << buffer_.size() << ", arr = " << arr_data_.size() << ", vec = " << vec_data_.size() << std::endl;
  //   }

  //   // Valarray constructor               Buffer
  //   //-----------------------------------------------------------------------------------
  //   Buffer(const std::string& name, const std::valarray<T>& data, const size_t size=0, const std::vector<T>& vec=std::vector<T>()) 
  //     : name_(name), buffer_(std::max(data.size(), size), 0), arr_data_(data), vec_data_(vec), vec_buffer_(0) //, data_(new arr_data<T>(data)) 
  //     { std::cout << "valarray: " << name_ << ", buf = " << buffer_.size() << ", arr = " << arr_data_.size() << ", vec = " << vec_data_.size() << std::endl;}
  //   //-----------------------------------------------------------------------------------

  //   // Vector constructor                 Buffer
  //   //-----------------------------------------------------------------------------------
  //   Buffer(const std::string& name, const std::vector<T>& data, const size_t size=0, const std::valarray<T>& arr=std::valarray<T>()) 
  //     : name_(name), buffer_(std::max(data.size(), size), 0), arr_data_(arr), vec_data_(data), vec_buffer_(1) //, data_(new vec_data<T>(data))
  //     { std::cout << "vector: " << name_ << ", buf = " << buffer_.size() << ", arr = " << arr_data_.size() << ", vec = " << vec_data_.size() << std::endl;}
  //   //-----------------------------------------------------------------------------------

  //   // ~Buffer() { delete data_; }
  //   // //                                    Buffer
  //   // //-----------------------------------------------------------------------------------
  //   // Buffer(const std::valarray<T>& data, const size_t size) : buffer_(size, 0), data_(data) { }
  //   // //-----------------------------------------------------------------------------------
    
  //   //                                    Buffer
  //   //-----------------------------------------------------------------------------------
  //   void update()  
  //   //-----------------------------------------------------------------------------------
  //   {
  //     if (vec_buffer_) {
  //       std::copy(std::begin(vec_data_), std::end(vec_data_), buffer_.begin()); 
  //     } else {
  //       std::copy(std::begin(arr_data_), std::end(arr_data_), buffer_.begin()); 
  //     }
  //   }

  //   //                                    Buffer
  //   //-----------------------------------------------------------------------------------
  //   const std::vector<T>& buffer() const { return buffer_; }
  //   //-----------------------------------------------------------------------------------

  //   // //                                    Buffer
  //   // //-----------------------------------------------------------------------------------
  //   // void push_back(const T& value) { buffer_.push_back(value); }
  //   // //-----------------------------------------------------------------------------------

  //   //                                    Buffer
  //   //-----------------------------------------------------------------------------------
  //   int buffer_size() const { return static_cast<int>(buffer_.size()); }
  //   //-----------------------------------------------------------------------------------

  //   //                                    Buffer
  //   //-----------------------------------------------------------------------------------
  //   const T* data_ptr() const 
  //   //-----------------------------------------------------------------------------------
  //   { 
  //     if (vec_buffer_) 
  //       return &vec_data_[0]; 
  //     else
  //       return &arr_data_[0]; 
  //   }
  // };

  // class Buffer_vector : public Buffer
  // {
  //   private:
  //   const std::vector<T>& data_;

  //   public:
  //   //                                    Buffer_vector
  //   //-----------------------------------------------------------------------------------
  //   Buffer(const std::valarray<T>& data) : buffer_(data.size(), 0), data_(data) { }
  //   //-----------------------------------------------------------------------------------

  //   //                                    Buffer_vector
  //   //-----------------------------------------------------------------------------------
  //   Buffer(const std::valarray<T>& data, const size_t size) : buffer_(size, 0), data_(data) { }
  //   //-----------------------------------------------------------------------------------

  // }


  //=====================================================================================
  //
  //                                    O U T P U T
  //
  //=====================================================================================

  template <typename CELL, typename T>
  class Output 
  {
    private:
    int format_ = BINARY;
    Grid<CELL> grid_;
    std::string path_;
    std::vector<Outfile<T>> outfiles_;
    std::unordered_map<std::string, int> get_index_;
    int nwrite_ = 0;
    int rank_ = 0;
    int max_rank_ = 0;
    // std::deque<Buffer<T>> buffers_;  // If vector is used instead of deque, we need to call reserve (to avoid reallocation) and specify a maximum number of buffers 
    //std::vector<T> empty_vec = std::vector<T>();
    //std::valarray<T> empty_arr = std::valarray<T>();
    std::vector< std::unique_ptr<data_wrapper<T>> > wrappers_;
    std::deque<std::vector<T>> buffer_2d_; // Buffers for 2D-data
    std::vector<int> buffer_index_;


    public:
    //                                     Output
    //-----------------------------------------------------------------------------------
    template <typename S>
    Output(int format, const std::vector<std::vector<S>>& nodes, const std::string path="out", int rank=0, int num_procs=1) 
    //-----------------------------------------------------------------------------------
      : format_(format), grid_(nodes, format), path_(path), outfiles_(), get_index_(), rank_(rank), max_rank_(num_procs-1), wrappers_(), buffer_2d_(), buffer_index_() { }

    // //                                     Output
    // //-----------------------------------------------------------------------------------
    // Outfile<T>& operator[](const std::string& name) { return outfiles_[get_index_[name]]; } 
    // //-----------------------------------------------------------------------------------

    // //                                     Output
    // //-----------------------------------------------------------------------------------
    // Outfile<T>& operator[](const int index) { return outfiles_[index]; }
    // //-----------------------------------------------------------------------------------

    // //                                     Output
    // //-----------------------------------------------------------------------------------
    // Outfile<T>& file(const int index) { return outfiles_[index]; }
    // //-----------------------------------------------------------------------------------

    //                                     Output
    //-----------------------------------------------------------------------------------
    Outfile<T>& last_file() const { return outfiles_.back(); }
    //-----------------------------------------------------------------------------------

    //                                     Output
    //-----------------------------------------------------------------------------------
    Outfile<T>& add_file(const std::string& name)
    //-----------------------------------------------------------------------------------
    {
      outfiles_.emplace_back(path_, name, grid_.offset());
      get_index_[name] = outfiles_.size()-1;
      return outfiles_.back();
    }

    // //-----------------------------------------------------------------------------------
    // const std::vector<Buffer<T>>& buffers() const { return buffers_; }
    // //-----------------------------------------------------------------------------------


    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variable(const std::string& name, int dim, const std::vector<T>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      wrappers_.emplace_back(new vec_wrapper<T>(data));
      add_variable_(name, dim, index, length, offset);
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variable(const std::string& name, int dim, const std::valarray<T>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      wrappers_.emplace_back(new arr_wrapper<T>(data));
      add_variable_(name, dim, index, length, offset);
    }

    // //                                     Output
    // //-----------------------------------------------------------------------------------
    // void add_variable(const std::string& name, int dim, const std::valarray<T>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    // //-----------------------------------------------------------------------------------
    // // Valarray version
    // // For valarray data a vector-copy (buffer) of the original data is made 
    // {
    //   // Check if data is previously buffered
    //   const Buffer<T>* buffer = nullptr;
    //   const T* data_ptr = &data[0];
    //   for (const auto& buf : buffers_) {
    //     if (buf.data_ptr() == data_ptr) {
    //       buffer = &buf;
    //       break;
    //     }
    //   }
    //   if (buffer == nullptr) {
    //     // Buffer does not exist, create it
    //     auto size = data.size();
    //     if (dim == 2)
    //       size += 1; 
    //     buffers_.emplace_back(name, data, size);
    //     buffers_.back().info();
    //     buffer = &(buffers_.back());
    //   } 

    //   if (dim == 2) {
    //     dim = 3;          
    //     add_variable(name, dim, buffer->buffer(), util::add_coord(index, buffer->buffer_size()-1, 2, 3), length, offset);
    //   } else {
    //     add_variable(name, dim, buffer->buffer(), index, length, offset);
    //   }
    // }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void write(double time)
    //-----------------------------------------------------------------------------------
    {
      //std::cout << "A" << std::endl;
      for (size_t i=0; i<buffer_2d_.size(); ++i ) {
        std::cout << "buffer: " << i << ", " << buffer_index_[i] << std::endl;
        const data_wrapper<T>& src_data = *(wrappers_[ buffer_index_[i] ]);
        const data_wrapper<T>& wrap_vec = *(wrappers_[ buffer_index_[i+1] ]);
        for (size_t i=0; i<src_data.size(); ++i)
          std::cout << src_data.at(i) << ", ";
        for (size_t i=0; i<src_data.size(); ++i)
          std::cout << wrap_vec.at(i) << ", ";
        std::cout << std::endl;
        std::vector<T>& dst_data = buffer_2d_[i];
        util::print_vector(dst_data, 2);
        std::cout << "A: " << dst_data.size() << std::endl;
        // const T* a = src_data.begin();
        // const T* b = src_data.end();
        // std::cout << *(a+2) << ", " << *(b-1) << std::endl;
        std::cout << &dst_data[0] << ", " << &dst_data[dst_data.size()-1] << std::endl;
        std::copy(src_data.begin(), src_data.end(), dst_data.begin());
        util::print_vector(dst_data, 2);
        std::cout << "B: " << dst_data.size() << std::endl;
        std::cout << &dst_data[0] << ", " << &dst_data[dst_data.size()-1] << std::endl;
        for (size_t i=0; i<src_data.size(); ++i)
          std::cout << wrap_vec.at(i) << ", ";
        //buffer.info();
        //buffer.update();
      }
      //std::cout << "B" << std::endl;
      for (auto& outfile : outfiles_) {
        outfile.write(grid_, time, rank_, max_rank_);    
      }
      ++nwrite_;
    }
  
  private:
    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variable_(const std::string& name, int dim, const std::vector<int>& index, int length, int offset)
    //-----------------------------------------------------------------------------------
    {
      const data_wrapper<T>& data = *(wrappers_.back());
      for (size_t i=0; i<data.size(); ++i)
        std::cout << data.at(i) << ", ";
      std::cout << std::endl;
      if (index.size() > 0 && index.size() != grid_.num_cells()*dim) {
        std::cerr << "  ERROR in VTK::Output::add_variable(" << name << ", " << dim << "):" << std::endl;
        std::cerr << "  Size of index-vector (" << index.size() <<  ") does not match number of grid cells X dim (" << grid_.num_cells()*dim << ")" << std::endl;
        util::safe_exit(EXIT_FAILURE);
      }
      if (dim == 2) {
        dim = 3;
        int size = data.size();
        //buffers_.emplace_back(name, data, size+1);
        std::vector<T> data2d(size+1);
        std::copy(data.begin(), data.end(), data2d.begin());
        buffer_2d_.push_back(data2d);
        buffer_index_.push_back(wrappers_.size()-1); // Store index of original data for later copy
        wrappers_.emplace_back(new vec_wrapper<T>(data2d));
        //auto data = wrappers_.back();
        // Fix index-vector
        std::vector<int> index_3d(index);
        if (index_3d.empty()) {
          // Create default contiguous index-vector
          auto ind = util::linspace(offset, offset+size);
          index_3d.insert(index_3d.begin(), ind.begin(), ind.end());
        }
        index_3d = util::add_coord(index_3d, size, 2, 3);
        util::print_vector(data2d, 2);
        //util::print_vector(index_3d, 3);
        outfiles_.back().variables().add(name, *(wrappers_.back()), format_, dim, index_3d, length, offset);
      } else {
        outfiles_.back().variables().add(name, data, format_, dim, index, length, offset);
      }

      //outfiles_.back().variables().add(name, data, format_, dim, index, length, offset);
      outfiles_.back().update_offset();
      //std::cout << "A" << std::endl;
    }

  };


}


#endif /* SRC_VTK_H_ */
