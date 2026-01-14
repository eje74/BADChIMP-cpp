#ifndef SRC_VTK_COMMON_H_
#define SRC_VTK_COMMON_H_
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
#include <memory>
#include <array>
#include <type_traits>
#include <limits>
#include <mpi.h>
//#include <sys/types.h>
#include <sys/stat.h>
//#if defined (_WIN32)
//#include <direct.h>
//#else
#include <unistd.h>
//#endif
//#include "vector_func.h"
//#define TIMER
#ifdef TIMER
#include <chrono>
#endif

namespace util {

  template <typename T>
  //-----------------------------------------------------------------------------------  
  // Generate a linearly spaced vector.
  std::vector<T> linspace(T start, T stop, T inc=1) {
  //-----------------------------------------------------------------------------------  
    std::vector<T> vec((int)((stop-start)/inc), inc);
    vec[0] = start;
    std::partial_sum(vec.begin(), vec.end(), vec.begin(), std::plus<T>());
    return vec;
  }

  template <typename T>
  //-----------------------------------------------------------------------------------
  // Pad a vector with an extra coordinate.
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
  // Print a vector grouped by dimension.
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
  // Convert a value to string.
  const std::string to_string(const T value) { return std::to_string(value); }
  //-----------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------
  // Passthrough for strings.
  const std::string to_string(const std::string value) { return value; }
  //-----------------------------------------------------------------------------------

  template <typename T>
  //-----------------------------------------------------------------------------------
  // Join vector elements into a string.
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
  // Abort with MPI finalize if needed.
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

  template <int DIM>
  //-----------------------------------------------------------------------------------
  // Find per-dimension min/max indices.
  void find_min_max(const std::vector<int>& pos, std::array<int, DIM>& min_pos, std::array<int, DIM>& max_pos)
  //-----------------------------------------------------------------------------------
  {
    min_pos.fill(std::numeric_limits<int>::max());
    max_pos.fill(std::numeric_limits<int>::min());
    const size_t node_count = pos.size()/DIM;
    for (size_t n=0; n<node_count; ++n) {
      for (int d=0; d<DIM; ++d) {
        const int val = pos[n*DIM + d];
        if (val < min_pos[d]) min_pos[d] = val;
        if (val > max_pos[d]) max_pos[d] = val;
      }
    }
  }

  template <int DIM>
  //-----------------------------------------------------------------------------------
  // Compute block size and expected node count.
  long long block_size(const std::array<int, DIM>& min_pos, const std::array<int, DIM>& max_pos, std::array<int, DIM>& size)
  //-----------------------------------------------------------------------------------
  {
    size.fill(1);
    long long expected = 1;
    for (int d=0; d<DIM; ++d) {
      size[d] = max_pos[d] - min_pos[d] + 1;
      expected *= size[d];
    }
    return expected;
  }

  template <int DIM>
  //-----------------------------------------------------------------------------------
  // Build a 6-value extent from min/max indices.
  void build_piece_extent(const std::array<int, DIM>& min_pos, const std::array<int, DIM>& max_pos, std::array<int, 6>& piece_extent)
  //-----------------------------------------------------------------------------------
  {
    piece_extent = {0, 0, 0, 0, 0, 0};
    piece_extent[0] = min_pos[0];
    piece_extent[1] = max_pos[0] + 1;
    piece_extent[2] = (DIM > 1) ? min_pos[1] : 0;
    piece_extent[3] = (DIM > 1) ? max_pos[1] + 1 : 0;
    piece_extent[4] = (DIM > 2) ? min_pos[2] : 0;
    piece_extent[5] = (DIM > 2) ? max_pos[2] + 1 : 0;
  }

}

namespace VTK {
  
  using point_dtype = float;
  //using nodes_vec = std::vector<std::vector<point_dtype>>;

  static constexpr int BINARY = 0;
  static constexpr int ASCII = 1;

  //-----------------------------------------------------------------------------------  
  //  Get datatype name
  //  Adapted from: https://stackoverflow.com/questions/4484982/how-to-convert-typename-t-to-string-in-c
  //-----------------------------------------------------------------------------------
  template <typename T>
  struct datatype {
    // Return VTK type name for T.
    static const char* name() {
      return typeid(T).name();
    } 
  };
  template <>
  struct datatype<int> {
    // Return VTK type name for int.
    static const char* name() {
      return "Int32";
    } 
  };
  template <>
  struct datatype<unsigned int> {
    // Return VTK type name for unsigned int.
    static const char* name() {
      return "UInt32";
    } 
  };
  template <>
  struct datatype<float> {
    // Return VTK type name for float.
    static const char* name() {
      return "Float32";
    } 
  };
  template <>
  struct datatype<double> {
    // Return VTK type name for double.
    static const char* name() {
      return "Float64";
    } 
  };


  //=====================================================================================
  //
  //                              D A T A _ W R A P P E R
  //
  //=====================================================================================
  // data_wrapper provides a uniform view of underlying data containers.
  // Cast wrappers keep a float buffer in sync when output type differs
  // from the input container element type.
  template <typename T>
  class data_wrapper
  {
      public:
      data_wrapper() { } 
      virtual ~data_wrapper() { };
      // Return pointer to the requested element.
      virtual const T* ptr(const int pos) const = 0;      
      // Return number of elements.
      virtual const size_t size() const = 0;
      // Sync backing storage if needed (no-op by default).
      virtual void sync() { }
      // Read element by index.
      const T at(const int pos) const { return *ptr(pos); }
      // Begin pointer for contiguous access.
      const T* begin() const { return ptr(0); }
      // End pointer for contiguous access.
      const T* end() const { return ptr(size()); }
  };

  template <typename T>
  class vec_wrapper : public data_wrapper<T>
  {
      private:
          const std::vector<T>& data_;
      public:
          vec_wrapper(const std::vector<T>& data) : data_(data) { }
          // Pointer to vector data.
          const T* ptr(const int pos) const { return &data_[pos]; }
          // Vector size.
          const size_t size() const { return data_.size(); }
  };

  template <typename T>
  class arr_wrapper : public data_wrapper<T>
  {
      private:
          const std::valarray<T>& data_;
      public:
          arr_wrapper(const std::valarray<T>& data) : data_(data) { }
          // Pointer to valarray data.
          const T* ptr(int pos) const { return &data_[pos]; }
          // Valarray size.
          const size_t size() const { return data_.size(); }
  };

  template <typename OutT, typename InT>
  class vec_cast_wrapper : public data_wrapper<OutT>
  {
      private:
          const std::vector<InT>& data_;
          std::vector<OutT> buffer_;
      public:
          vec_cast_wrapper(const std::vector<InT>& data) : data_(data), buffer_(data.size()) { sync(); }
          // Pointer to casted buffer.
          const OutT* ptr(const int pos) const { return &buffer_[pos]; }
          // Buffer size.
          const size_t size() const { return buffer_.size(); }
          // Update buffer with casted values.
          void sync() {
            if (buffer_.size() != data_.size()) {
              buffer_.resize(data_.size());
            }
            for (size_t i = 0; i < data_.size(); ++i) {
              buffer_[i] = static_cast<OutT>(data_[i]);
            }
          }
  };

  template <typename OutT, typename InT>
  class arr_cast_wrapper : public data_wrapper<OutT>
  {
      private:
          const std::valarray<InT>& data_;
          std::vector<OutT> buffer_;
      public:
          arr_cast_wrapper(const std::valarray<InT>& data) : data_(data), buffer_(data.size()) { sync(); }
          // Pointer to casted buffer.
          const OutT* ptr(const int pos) const { return &buffer_[pos]; }
          // Buffer size.
          const size_t size() const { return buffer_.size(); }
          // Update buffer with casted values.
          void sync() {
            if (buffer_.size() != data_.size()) {
              buffer_.resize(data_.size());
            }
            for (size_t i = 0; i < data_.size(); ++i) {
              buffer_[i] = static_cast<OutT>(data_[i]);
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
    // Default constructor.
    DataArray() : tags_() { }
    //-----------------------------------------------------------------------------------

    //                                DataArray
    //-----------------------------------------------------------------------------------
    // Construct DataArray with name, dimension, and format metadata.
    DataArray(const std::string& name, int dim, const int format=0, const std::string& type="", int offset=-1)
    //-----------------------------------------------------------------------------------
     : tags_() 
    {
      //std::cout << "DataArray" << std::endl;
      set_tag("name", name);
      set_tag("type", type);
      set_tag("dim", dim);
      set_tag("format", data_formats_[format]);
      set_tag("offset", (format==ASCII) ? -1 : offset );
    }

    //                                DataArray
    //-----------------------------------------------------------------------------------
    // Quote and format tag values.
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
    const data_wrapper<T>& data_;
    int dim_ = 0;
    unsigned int offset_ = 0;
    unsigned int length_ = 0;
    DataArray dataarray_;
    unsigned int nbytes_ = 0;
    std::vector<int> index_;
    int contiguous_ = 1; 
    double min_value_ = 1e-20;
    std::string indent_;

    public:
    //                                   Data
    //-----------------------------------------------------------------------------------
    // Default constructor.
    Data() : data_(std::vector<T>()), dataarray_(), index_() {};                                  
    //-----------------------------------------------------------------------------------

    //                                   Data
    //-----------------------------------------------------------------------------------
    // Construct Data with data wrapper and indexing metadata.
    Data(const std::string &name, const data_wrapper<T>& data, const int format, const int dim, const int level=0, const std::vector<int>& index=std::vector<int>(), const int length=0, const int offset=0) 
    //-----------------------------------------------------------------------------------
      : name_(name), data_(data), dim_(dim), offset_(offset), length_(length), dataarray_(name, dim, format, datatype<T>::name()), index_(index), indent_(level*2, ' ')  
    { 
      // std::cout << "Data" << std::endl;
      if (index.empty()) {
        // Assume contiguous data, create default index-vector 
        contiguous_ = 1;
        index_ = util::linspace<int>(offset_, offset_+data_.size());
      } else {
        // Index-vector is provided, data is non-contiguous
        contiguous_ = 0;
      }
      if (dim == 2) {
        // Paraview gives an error message for 2D-data. A third coordinate (0)
        // is added to fix this. The data-vector is left unchanged, but the index-vector
        // is given a negative index for the third coordinate. The write-functions write  
        // zero if the index is negative.    
        dim_ = 3;
        dataarray_.set_tag("dim", dim_);  // Update the data-array
        contiguous_ = 0;
        index_ = util::add_coord(index_, -1, 2, 3); //ind3d;
      }
      if (length_ == 0) {
        length_ = index_.size();
      }
      if (format == BINARY) {
        nbytes_ = length_*sizeof(T);
      }
      // std::cout << dataarray_ << ", data-size: " << data_.size() << ", index-size: " << index_.size() << std::endl; 
    }  

    //                                   Data
    //-----------------------------------------------------------------------------------
    // Return number of elements.
    size_t size() const { return data_.size(); } 
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    // Return variable name.
    const std::string& name() const { return name_; } 
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    // Return variable dimension (components).
    int dim() const { return dim_; } 
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    // Report if the data is written as binary.
    bool is_binary() const {return (nbytes_ > 0) ? true : false; }
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    // Write ASCII payload.
    void write_asciidata (std::ofstream& file) const
    //-----------------------------------------------------------------------------------
    {
      if (not is_binary()) {
        for (const auto& ind : index_) {
          if (ind<0) {
            file << "0 ";
          } else {
            T val = data_.at(ind);
            if (std::abs(val)<min_value_) {
              file << "0 ";
            } else {
              file << val << " ";
            }
          }
        }        
      }
    }

    //                                   Data
    //-----------------------------------------------------------------------------------
    // Write appended binary payload.
    void write_binarydata(std::ofstream& file) const {
    //-----------------------------------------------------------------------------------
      if (is_binary()) {
        file.write((char*)&nbytes_, sizeof(unsigned int));
        if (contiguous_) {
          // file.write((char*)&data_[offset_], nbytes_);
          file.write((char*)data_.ptr(offset_), nbytes_);
        } else {          
          T zero = 0;
          for (const auto& ind : index_) {
            if (ind<0) 
              file.write((char*)&zero, sizeof(T));            
            else             
              file.write((char*)data_.ptr(ind), sizeof(T));            
          }
        }
      }
    }
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    // Access DataArray tag values.
    const std::string& dataarray(const std::string& tname) const { return dataarray_[tname]; }
    //-----------------------------------------------------------------------------------
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    // Update and return the next appended offset.
    int update_offset(int offset)
    //-----------------------------------------------------------------------------------
    {
      if (is_binary() && dataarray_["offset"].empty()) {
        dataarray_.set_tag("offset", offset); 
        return offset + nbytes_ + sizeof(unsigned int);
      } else {
        return offset;
      }
    }
    
    //                                   Data
    //-----------------------------------------------------------------------------------
    // Write the XML DataArray element.
    void write_dataarray(std::ofstream& file, int indent=0) const
    //-----------------------------------------------------------------------------------
    { 
      file << indent_ << "<DataArray " << dataarray_ << "> ";
      write_asciidata(file);  
      file << "</DataArray>" << std::endl;
    }

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
    static constexpr int indent_ = 4;

    //                                   Variables
    //-----------------------------------------------------------------------------------
    // Default constructor.
    Variables() : datalist_(), scalars_(), vectors_() { }    
    //-----------------------------------------------------------------------------------

    //                                   Variables
    //-----------------------------------------------------------------------------------
    // Register a variable and its data wrapper.
    void add(const std::string& name, const data_wrapper<T>& data, const int format, const int dim, const std::vector<int>& index, const int length=0, const int offset=0) 
    //-----------------------------------------------------------------------------------
    {
      // std::cout << "Variables" << std::endl;
      datalist_.emplace_back(name, data, format, dim, indent_, index, length, offset);      
      update_names(datalist_.back());
    }

    //                                   Variables
    //-----------------------------------------------------------------------------------
    // Access all variables.
    const std::vector<Data<T>>& data() const { return datalist_; }
    //-----------------------------------------------------------------------------------

    //                                   Variables
    //-----------------------------------------------------------------------------------
    // Access the most recently added variable.
    const Data<T>& back() const { return datalist_.back(); }
    Data<T>& back() { return datalist_.back(); }
    //-----------------------------------------------------------------------------------

    //                                   Variables
    //-----------------------------------------------------------------------------------
    // Return vector variable names.
    const std::string& vector_names() const { return vector_names_; }
    //-----------------------------------------------------------------------------------

    //                                   Variables
    //-----------------------------------------------------------------------------------
    // Return scalar variable names.
    const std::string& scalar_names() const { return scalar_names_; }
    //-----------------------------------------------------------------------------------

    private:
    //                                   Variables
    //-----------------------------------------------------------------------------------
    // Track scalar/vector names for XML tags.
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
    // Default constructor.
    File() : file_(), folders_() { }
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    // Construct and create output directories.
    File(const std::string &name, const std::vector<std::string> &folders, const std::string &extension)
    : name_(name), file_(), folders_(folders), extension_(extension)
    //-----------------------------------------------------------------------------------
    {
      // std::cout << "File" << std::endl;
      for (const auto& f : folders_) {
        path_ += f;
        make_dir(path_);
      }
    }

    //                                    File
    //-----------------------------------------------------------------------------------
    // Open output file.
    void open() { file_.open(path_+filename_, std::ios::out); }
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    // Close output file with final XML tag.
    void close() {
    //-----------------------------------------------------------------------------------
      file_ << "</VTKFile>" << std::endl;
      file_.close();
    }

    //                                    File
    //-----------------------------------------------------------------------------------
    // Return filename (no path).
    const std::string& filename() const { return filename_;}
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    // Access file stream (read-only).
    const std::ofstream& file() const { return file_; }
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    // Return full path to file.
    const std::string path() const { return folders_.back()+filename_; }
    //-----------------------------------------------------------------------------------

    //                                    File
    //-----------------------------------------------------------------------------------
    // Increment write counter.
    void inc_nwrite() { ++nwrite_; }
    //-----------------------------------------------------------------------------------

    private:
    //                                    File
    //-----------------------------------------------------------------------------------
    // Ensure directory exists.
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


}

#endif /* SRC_VTK_COMMON_H_ */
