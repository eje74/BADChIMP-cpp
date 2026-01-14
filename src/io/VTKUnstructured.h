#ifndef SRC_VTK_UNSTRUCTURED_H_
#define SRC_VTK_UNSTRUCTURED_H_
#include "VTKCommon.h"

namespace VTK {

  //-----------------------------------------------------------------------------------
  // 1: VTK_VERTEX
  //            o
  //-----------------------------------------------------------------------------------
  struct vertex {
    static constexpr char name[] = "Vertex";
    static constexpr int type = 1;
    static constexpr bool center = false;
    static constexpr int n = 1;
    static constexpr std::array<std::array<int,3>, n> points = { {{0,0,0}} };
  };
  // These three lines are necessary to avoid linker errors in c++11, but can be removed for c++17
  // constexpr std::array<std::array<int, vertex::dim>, vertex::n> vertex::points;
  // constexpr int vertex::type;
  // constexpr char vertex::name[];
  
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
    static constexpr bool center = false;
    static constexpr int n = 2;
    static constexpr std::array<std::array<int,3>, n> points = { {{0,0,0}, {1,0,0}} };
  };

  // These three lines are necessary to avoid linker errors in c++11, but can be removed in c++17
  // constexpr std::array<std::array<int, line::dim>, line::n> line::points;
  // constexpr int line::type;
  // constexpr char line::name[];

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
    static constexpr bool center = true;
    static constexpr int n = 4;
    static constexpr std::array<std::array<int,3>, n> points = {{{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}}};
  };
  // These three lines are necessary to avoid linker errors in c++11, but can be removed in c++17
  // constexpr std::array<std::array<int, pixel::dim>, pixel::n> pixel::points;
  // constexpr int pixel::type;
  // constexpr char pixel::name[];

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
  //                  o
  //                  0
  //----------------------------------------------------------------------------------- 
  */
  struct quad {
    static constexpr char name[] = "Quad";
    static constexpr bool center = true;
    static constexpr int type = 9;
    static constexpr int n = 4;
    static constexpr std::array<std::array<int,3>, n> points = {{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}}};
  };
  // These three lines are necessary to avoid linker errors in c++11, but can be removed for c++17
  // constexpr std::array<std::array<int, quad::dim>, quad::n> quad::points;
  // constexpr int quad::type;
  // constexpr char quad::name[];

  //----------------------------------------------------------------------------------- 
  // 11: VTK_VOXEL
  //
  //            6 o--------o 7
  //             /|       /|
  //            / |      / |
  //         4 o--------o 5 |
  //           |  |     |  |
  //           | 2 o----|---o 3
  //           | /      | /
  //           |/       |/
  //         0 o--------o 1
  //                      
  //----------------------------------------------------------------------------------- 
  struct voxel {
    static constexpr char name[] = "Voxel";
    static constexpr bool center = true;
    static constexpr int type = 11;
    static constexpr int dim = 3;
    static constexpr int n = 8;
    static constexpr std::array<std::array<int, dim>, n> points = {{{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}}};
  };
  // These three lines are necessary to avoid linker errors in c++11, but can be removed for c++17
  // constexpr std::array<std::array<int, voxel::dim>, voxel::n> voxel::points;
  // constexpr int voxel::type;
  // constexpr char voxel::name[];

  //-----------------------------------------------------------------------------------
  // Lattice-Boltzmann 3D lattice
  //-----------------------------------------------------------------------------------
    struct D3Q19 {
    static constexpr char name[] = "D3Q19";
    static constexpr bool center = false;
    static constexpr int type = vertex::type;
    static constexpr int n = 19;
    static constexpr std::array<std::array<int,3>, n> points = { {{1,0,0}, {0,1,0}, {0,0,1}, {1,1,0}, {1,-1,0}, {1,0,1}, {1,0,-1}, 
                                                                  {0,1,1}, {0,1,-1}, {-1,0,0}, {0,-1,0}, {0,0,-1}, {-1,-1,0}, {-1,1,0}, 
                                                                  {-1,0,-1}, {-1,0,1}, {0,-1,-1}, {0,-1,1}, {0,0,0}} };
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

  // Use template specialization to choose a 2D pixel if a 2D voxel is requested
  template <typename CELL, int DIM>
  struct Cell 
  { 
    typedef CELL cell; 
    static constexpr int dim = DIM; 
  };
  template<>
  struct Cell<voxel,2> 
  { 
    typedef pixel cell; 
    static constexpr int dim = 2; 
  };
  
  
  // Use template CELL instead of C and D to get effect of 
  // the template specialization above
  template <typename C, int D, typename CELL=Cell<C,D>>
  class Mesh
  {
  private:
    Mesh_data<point_dtype> points_;
    Mesh_data<int> conn_;
    Mesh_data<int> offsets_;
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
    Mesh(const std::vector<T> &nodes, point_dtype unit = 1) 
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
    // int num() const { return int(conn_.vec_.size()/CELL::n); }
    int num() const { return int(conn_.vec_.size()/CELL::cell::n); }

  private:
    //                                   Mesh
    //-----------------------------------------------------------------------------------
    void init() 
    //-----------------------------------------------------------------------------------
    {
      offsets_.vec_ = util::linspace<int>(CELL::cell::n, conn_.vec_.size()+CELL::cell::n, CELL::cell::n);
      types_.vec_ = std::vector<int>(num(), CELL::cell::type);
      // std::cout << "type, dim, n: " << CELL::type << ", " << CELL::dim << ", " << CELL::n << std::endl;
      // std::cout << "points: "; util::print_vector(points_, point_dim_);
      // std::cout << "connectivity: "; util::print_vector(conn_, CELL::n);
      // std::cout << "offsets: " << offsets_ << std::endl;
      // std::cout << "types: " << types_ << std::endl;
    }

    //                                   Mesh
    //-----------------------------------------------------------------------------------
    template <typename T>
    void set_size(const std::vector<T> &nodes)
    //-----------------------------------------------------------------------------------
    {
      std::vector<T> min_n(CELL::dim, 999999), max_n(CELL::dim, -1);
      int n = 0;
      for (const auto node : nodes) {
        int i = n%CELL::dim;
        if (node > max_n[i])
          max_n[i] = node;
        if (node < min_n[i])
          min_n[i] = node;
        ++n;
      }
      std::vector<int> size(CELL::dim, 0); 
      for (auto i = 0; i < CELL::dim; ++i) {
        size[i] = int(max_n[i] - min_n[i]) + 1;
      }
      size_ = size;
      // util::print_vector(size_, CELL::dim);
    }

    //                                   Mesh
    //-----------------------------------------------------------------------------------
    std::vector<int> get_stride() const 
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
    void calc_points_and_conn(const std::vector<T> &nodes) 
    //-----------------------------------------------------------------------------------
    {
      // if (nodes[0].size() != CELL::dim*CELL::cell::n) {
      //   std::cerr << "ERROR in VTK::Mesh: A " << CELL::name << " is " << CELL::dim << "-dimensional, but the given nodes are " << nodes[0].size() << "-dimensional" << std::endl;        
      //   util::safe_exit(EXIT_FAILURE);
      // }
      // Loop over nodes and calculate corner points and unique indexes 
      int num_nodes = static_cast<int>(nodes.size()/CELL::dim);
      auto stride = get_stride();
      std::vector<point_dtype> pts; 
      std::vector<int> index;
      pts.reserve(CELL::cell::n * num_nodes * CELL::dim);
      index.reserve(CELL::cell::n * num_nodes);
      for (int n=0; n<num_nodes; ++n) {
        for (size_t p=0; p<CELL::cell::n; ++p) {
          int idx = 0;
          for (auto i=0; i < CELL::dim; ++i) {
            auto coord = nodes[n*CELL::dim +i] + CELL::cell::points[p][i];
            idx += coord*stride[i];
            pts.push_back(coord);
          }
          index.push_back(idx);
        }
      }
      // Create unique index_list 
      auto minmax = std::minmax_element(index.begin(), index.end());
      auto min = *minmax.first;
      auto max = *minmax.second;
      // std::cout << "min, max = " << min << ", " << max << std::endl;
      std::vector<int> index_list(max - min + 1, -1);
      conn_.vec_.reserve(CELL::cell::n * num_nodes);
      std::vector<int> pos_index;
      pos_index.reserve(CELL::cell::n * num_nodes);
      int n_pts = 0;
      int pos = 0;
      for (size_t i=0; i<index.size(); ++i) {
        int idx = index[i] - min;
        if (index_list[idx] < 0) {
          pos_index.push_back(pos);
          index_list[idx] = n_pts++;
        }
        pos += CELL::dim;
        conn_.vec_.push_back(index_list[idx]);
      }

      points_.vec_.reserve(n_pts*point_dim_);
      double shift = 0;
      // Shift points to center nodes inside the cell
      if (CELL::cell::center)
        shift = 0.5*unit_;
      for (const auto p : pos_index) {
        int i, ii;
        for (i=0; i<CELL::dim; ++i)
          points_.vec_.push_back(pts[p+i]-shift);
        // Pad with 0 if dim < 3
        for (ii=i; ii<point_dim_; ++ii)
          points_.vec_.push_back(0);
      }

    }
  };


  //=====================================================================================
  //
  //                                     G R I D
  //
  //=====================================================================================
  template <typename CELL, int DIM>
  class Grid {
    private:
    Mesh<CELL, DIM> mesh_;
    Data<point_dtype> point_data_;
    std::vector<Data<int>> cell_data_;
    unsigned int offset_ = 0;
    static constexpr int indent_ = 4;

    public:
    //                                   Grid
    //-----------------------------------------------------------------------------------
    Grid() : mesh_(), point_data_(), cell_data_() { }
    //-----------------------------------------------------------------------------------

    //                                   Grid
    //-----------------------------------------------------------------------------------
    template <typename S>
    Grid(const std::vector<S>& nodes, int format=BINARY) 
    : mesh_(nodes), point_data_("points", mesh_.points(), format, mesh_.dim(), indent_), cell_data_()
    //-----------------------------------------------------------------------------------
    {
      // std::cout << "Grid" << std::endl;
      cell_data_.emplace_back("connectivity", mesh_.connectivity(), format, 1, indent_);
      cell_data_.emplace_back("offsets",      mesh_.offsets()     , format, 1, indent_);
      cell_data_.emplace_back("types",        mesh_.types()       , format, 1, indent_);
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
    template <typename CELL, typename T, int DIM>
    void write(const int rank, const Grid<CELL,DIM>& grid, Variables<T>& var) 
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
    template <typename CELL, int DIM>
    void write_header(const Grid<CELL,DIM>& grid) {
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
      file_ << "      <CellData Scalars=\"" << var.scalar_names() << "\" Vectors=\"" << var.vector_names() << "\">" << std::endl;
      for (auto& data : var.data()) {
        data.write_dataarray(file_);
      }
      file_ << "      </CellData>" << std::endl;
    }

    //                                   VTU_file
    //-----------------------------------------------------------------------------------
    template <typename CELL, typename T, int DIM>
    void write_appended_data(const Grid<CELL,DIM>& grid, const Variables<T>& var) {
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
      //std::cout << "PVTU_File" << std::endl;
      // Check if this is a MPI-run
      MPI_Initialized(&mpi_running_); 
    }

    //                                   PVTU_file
    //-----------------------------------------------------------------------------------
    template <typename CELL, typename T, int DIM>
    void write(const double time, const int rank, const int max_rank, const Grid<CELL,DIM>& grid, const Variables<T>& var, const std::string& vtu_name)
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
        util::safe_exit(-1);
        // exit(-1);
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
    template <typename CELL, typename T, int DIM>
    void write_header(double time, const Grid<CELL,DIM>& grid, const Variables<T>& var) {
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
  class OutfileUnstructured {
    private:
    Variables<T> variables_;
    PVTU_file pvtu_file_;
    VTU_file vtu_file_;

    public:
    //                                  Outfile
    //-----------------------------------------------------------------------------------
    OutfileUnstructured(std::string &_path, const std::string &_name, const unsigned int _offset) 
      : variables_(), pvtu_file_(_path, _name), vtu_file_(_path, _name, _offset) { }
    //-----------------------------------------------------------------------------------

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    template <typename CELL, int DIM>
    void write(const Grid<CELL,DIM>& grid, double time, int rank, int max_rank)
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
    Variables<T>& variables() { return variables_; }
    //-----------------------------------------------------------------------------------

  };


  //=====================================================================================
  //
  //                                    O U T P U T
  //
  //=====================================================================================
  // VTK::OutputUnstructured manages one or more VTK files and the variables attached to them.
  // The template parameter T determines the on-disk data type (Float32/Float64).
  template <typename CELL, int DIM=3, typename T=double>
  class OutputUnstructured 
  {
    private:
    int format_ = BINARY;
    Grid<CELL,DIM> grid_;
    std::string path_;
    std::vector<OutfileUnstructured<T>> outfiles_;
    std::unordered_map<std::string, int> get_index_;
    int nwrite_ = 0;
    int rank_ = 0;
    int max_rank_ = 0;
    std::vector< std::unique_ptr<data_wrapper<T>> > wrappers_;
    //int dim_ = 3;

    public:
    //                                     Output
    //-----------------------------------------------------------------------------------
    template <typename S>
    OutputUnstructured(int format, const std::vector<S>& nodes, const std::string path="out", int rank=0, int num_procs=1) 
    //-----------------------------------------------------------------------------------
      : format_(format), grid_(nodes, format), path_(path), outfiles_(), get_index_(), rank_(rank), max_rank_(num_procs-1), wrappers_() { }


    //                                     Output
    //-----------------------------------------------------------------------------------
    OutfileUnstructured<T>& add_file(const std::string& name)
    //-----------------------------------------------------------------------------------
    {
      outfiles_.emplace_back(path_, name, grid_.offset());
      get_index_[name] = outfiles_.size()-1;
      return outfiles_.back();
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variable(const std::string& name, const std::vector<T>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      // wrappers_.emplace_back(new vec_wrapper<T>(data)); // c++11 version
      wrappers_.emplace_back(std::make_unique< vec_wrapper<T> >(data));
      add_variable_(name, index, length, offset);
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    // Enabled only when InT != T so mismatched input types go through the cast wrapper.
    template <typename InT, typename std::enable_if<!std::is_same<InT, T>::value, int>::type = 0>
    void add_variable(const std::string& name, const std::vector<InT>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      wrappers_.emplace_back(std::make_unique< vec_cast_wrapper<T, InT> >(data));
      add_variable_(name, index, length, offset);
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variable(const std::string& name, const std::valarray<T>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      // wrappers_.emplace_back(new arr_wrapper<T>(data));
      wrappers_.emplace_back(std::make_unique< arr_wrapper<T> >(data));
      add_variable_(name, index, length, offset);
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    // Enabled only when InT != T so mismatched input types go through the cast wrapper.
    template <typename InT, typename std::enable_if<!std::is_same<InT, T>::value, int>::type = 0>
    void add_variable(const std::string& name, const std::valarray<InT>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      wrappers_.emplace_back(std::make_unique< arr_cast_wrapper<T, InT> >(data));
      add_variable_(name, index, length, offset);
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void write(double time=0.0)
    //-----------------------------------------------------------------------------------
    {
#ifdef TIMER
      std::chrono::steady_clock::time_point begin =  std::chrono::steady_clock::now();
#endif
      // Ensure casted buffers are updated before writing.
      for (auto& wrapper : wrappers_) {
        wrapper->sync();
      }
      for (auto& outfile : outfiles_) {
        outfile.write(grid_, time, rank_, max_rank_);    
      }
      ++nwrite_;
#ifdef TIMER
      std::chrono::steady_clock::time_point end =  std::chrono::steady_clock::now();
      std::cout << "VTK::OutputUnstructured::write() duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << std::endl;
#endif
    }
  
  private:
    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variable_(const std::string& name, const std::vector<int>& index, int length, int offset)
    //-----------------------------------------------------------------------------------
    {
      int size = 1;
      if (index.size() > 0) {
        size = index.size();
      } else {
        size = (*wrappers_.back()).size();
      }
      int dim = size/grid_.num_cells();
      //std::cout << name << " : " << size << ", " << grid_.num_cells() << ", " << grid_.num_points() << std::endl;
      // if (index.size() > 0 && index.size() != grid_.num_cells()*dim) {
      if (dim != 1 and dim != DIM) {
        std::cerr << "  ERROR in VTK::OutputUnstructured::add_variable(" << name << ", " << dim << "):" << std::endl;
        std::cerr << "  Wrong dimension: Expected " << DIM << " or 1, got " << dim << std::endl;
        util::safe_exit(EXIT_FAILURE);
      }
      if (outfiles_.empty()) {
        std::cerr << "  ERROR in VTK::OutputUnstructured::add_variable(" << name << ", " << dim << "):" << std::endl;
        std::cerr << "  Add an output file using 'add_file(name)' before adding variables" << std::endl;
        util::safe_exit(EXIT_FAILURE);
      }
      outfiles_.back().variables().add(name, *wrappers_.back(), format_, dim, index, length, offset);
      outfiles_.back().update_offset();
    }

  };

}

#endif /* SRC_VTK_UNSTRUCTURED_H_ */
