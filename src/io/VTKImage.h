#ifndef SRC_VTK_IMAGE_H_
#define SRC_VTK_IMAGE_H_
#include "VTKCommon.h"
#include <unordered_set>
#include <limits>

namespace VTK {

  //=====================================================================================
  //
  //                                I M A G E  G R I D
  //
  //=====================================================================================
  // ImageGrid uses cell counts per dimension. Cell data length must match
  // cells_[0]*cells_[1]*(cells_[2] if DIM==3).
  template <int DIM>
  class ImageGrid {
    private:
    std::array<int, DIM> cells_;
    std::array<double, 3> origin_;
    std::array<double, 3> spacing_;
    std::array<int, 6> whole_extent_;
    std::array<int, 6> piece_extent_;
    int num_cells_ = 0;
    int num_points_ = 0;

    private:
    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Pad a DIM-sized vector to 3D for VTK metadata.
    static std::array<double, 3> pad_vec(const std::array<double, DIM>& vec, double def)
    //-----------------------------------------------------------------------------------
    {
      std::array<double, 3> out = {def, def, def};
      for (int i=0; i<DIM; ++i)
        out[i] = vec[i];
      return out;
    }

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Build a VTK extent array from cell counts.
    static std::array<int, 6> make_extent(const std::array<int, DIM>& cells)
    //-----------------------------------------------------------------------------------
    {
      std::array<int, 6> ext = {0, 0, 0, 0, 0, 0};
      ext[0] = 0;
      ext[1] = cells[0];
      ext[2] = 0;
      ext[3] = (DIM > 1) ? cells[1] : 0;
      ext[4] = 0;
      ext[5] = (DIM > 2) ? cells[2] : 0;
      return ext;
    }

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Update cached cell/point counts from cells_.
    void update_counts()
    //-----------------------------------------------------------------------------------
    {
      num_cells_ = 1;
      num_points_ = 1;
      for (int i=0; i<DIM; ++i) {
        num_cells_ *= cells_[i];
        num_points_ *= (cells_[i] + 1);
      }
    }

    public:
    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Construct an image grid with optional piece extent (defaults to whole extent).
    ImageGrid(const std::array<int, DIM>& cells, const std::array<double, DIM>& origin, const std::array<double, DIM>& spacing, const std::array<int, 6>& piece_extent)
    //-----------------------------------------------------------------------------------
      : cells_(cells), origin_(pad_vec(origin, 0.0)), spacing_(pad_vec(spacing, 1.0)), whole_extent_(make_extent(cells_)), piece_extent_(piece_extent)
    {
      update_counts();
      if (piece_extent_ == std::array<int, 6>{0, 0, 0, 0, 0, 0}) {
        piece_extent_ = whole_extent_;
      }
    }

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Return total number of cells.
    int num_cells() const { return num_cells_; }
    //-----------------------------------------------------------------------------------

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Return total number of points.
    int num_points() const { return num_points_; }

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Return number of cells in the current piece.
    int num_piece_cells() const {
    //-----------------------------------------------------------------------------------
      int nx = piece_extent_[1] - piece_extent_[0];
      int ny = (DIM > 1) ? (piece_extent_[3] - piece_extent_[2]) : 1;
      int nz = (DIM > 2) ? (piece_extent_[5] - piece_extent_[4]) : 1;
      return nx*ny*nz;
    }
    //-----------------------------------------------------------------------------------

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Return the global (whole) extent.
    const std::array<int, 6>& whole_extent() const { return whole_extent_; }
    //-----------------------------------------------------------------------------------

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Return the per-piece extent.
    const std::array<int, 6>& piece_extent() const { return piece_extent_; }
    //-----------------------------------------------------------------------------------

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Return origin padded to 3D.
    const std::array<double, 3>& origin() const { return origin_; }
    //-----------------------------------------------------------------------------------

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Return spacing padded to 3D.
    const std::array<double, 3>& spacing() const { return spacing_; }
    //-----------------------------------------------------------------------------------

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Format an extent array for XML attributes.
    static std::string extent_string(const std::array<int, 6>& ext)
    //-----------------------------------------------------------------------------------
    {
      std::ostringstream ss;
      ss << ext[0] << " " << ext[1] << " " << ext[2] << " " << ext[3] << " " << ext[4] << " " << ext[5];
      return ss.str();
    }

    //                                   ImageGrid
    //-----------------------------------------------------------------------------------
    // Format a vector (origin/spacing) for XML attributes.
    static std::string vec_string(const std::array<double, 3>& vec)
    //-----------------------------------------------------------------------------------
    {
      std::ostringstream ss;
      ss << vec[0] << " " << vec[1] << " " << vec[2];
      return ss.str();
    }
  };


  //=====================================================================================
  //                                  V T I  F I L E
  //                
  // VTK Serial XML file for Image Data
  //=====================================================================================
  class VTI_file : public File {
    private:
    unsigned int offset_ = 0;
    
    public:
    //                                   VTI_file
    //-----------------------------------------------------------------------------------
    // Initialize .vti writer for image data.
    VTI_file(const std::string &_path, const std::string &_name, unsigned int offset = 0) 
    //-----------------------------------------------------------------------------------
      : File(_name, {_path, "vti/"}, ".vti"), offset_(offset) { }

    //                                   VTI_file
    //-----------------------------------------------------------------------------------
    // Write a single .vti file and advance counters.
    template <int DIM, typename T>
    void write(const int rank, const ImageGrid<DIM>& grid, Variables<T>& var) 
    //-----------------------------------------------------------------------------------
    {
      set_filename_and_open(rank);
      write_header(grid);
      write_data(var);
      write_footer();
      write_appended_data(var); 
      close();
      inc_nwrite();
    }

    //                                   VTI_file
    //-----------------------------------------------------------------------------------
    // Set the appended-data offset.
    void set_offset(unsigned int offset) { offset_ = offset; }
    //-----------------------------------------------------------------------------------

    //                                   VTI_file
    //-----------------------------------------------------------------------------------
    // Get the appended-data offset.
    unsigned int offset() const { return offset_; }
    //-----------------------------------------------------------------------------------

    //                                   VTI_file
    //-----------------------------------------------------------------------------------
    // Write the VTI header and piece tags.
    template <int DIM>
    void write_header(const ImageGrid<DIM>& grid) {
    //-----------------------------------------------------------------------------------
      file_ << "<?xml version=\"1.0\"?>" << std::endl;
      file_ << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      file_ << "  <ImageData WholeExtent=\"" << grid.extent_string(grid.whole_extent()) << "\" Origin=\"" << grid.vec_string(grid.origin()) << "\" Spacing=\"" << grid.vec_string(grid.spacing()) << "\">" << std::endl;
      file_ << "    <Piece Extent=\"" << grid.extent_string(grid.piece_extent()) << "\">" << std::endl;
    }

    //                                   VTI_file
    //-----------------------------------------------------------------------------------
    // Write CellData tags for variables.
    template <typename T>
    void write_data(const Variables<T>& var) {
    //-----------------------------------------------------------------------------------
      file_ << "      <CellData Scalars=\"" << var.scalar_names() << "\" Vectors=\"" << var.vector_names() << "\">" << std::endl;
      for (auto& data : var.data()) {
        data.write_dataarray(file_);
      }
      file_ << "      </CellData>" << std::endl;
    }

    //                                   VTI_file
    //-----------------------------------------------------------------------------------
    // Write appended binary payload for variables.
    template <typename T>
    void write_appended_data(const Variables<T>& var) {
    //-----------------------------------------------------------------------------------
      if (offset_ > 0) {
        file_ << "  <AppendedData encoding=\"raw\">" << std::endl;
        file_ << "_";
        for (auto& data : var.data()) {
          data.write_binarydata(file_);
        }
        file_ << "  </AppendedData>" << std::endl;
      }
    }

    //                                   VTI_file
    //-----------------------------------------------------------------------------------
    // Close ImageData tags.
    void write_footer() {
    //-----------------------------------------------------------------------------------
      file_ << "    </Piece>" << std::endl;
      file_ << "  </ImageData>" << std::endl;
    }

    //                                   VTI_file
    //-----------------------------------------------------------------------------------
    // Build filename and open the file.
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
  //                                P V T I  F I L E
  //                
  //                    VTK Parallel XML file for Image Data
  //=====================================================================================
  class PVTI_file : public File {
    private:
    long file_position = 0;
    int mpi_running_ = 0;

    public:
    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Initialize .pvti writer for parallel image data.
    PVTI_file(const std::string &_path, const std::string &_name) : File(_name, {_path}, ".pvti")
    //-----------------------------------------------------------------------------------
    { 
      //std::cout << "PVTI_File" << std::endl;
      // Check if this is a MPI-run
      MPI_Initialized(&mpi_running_); 
    }

    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Write the .pvti file and register the piece entry.
    template <int DIM, typename T>
    void write(const double time, const int rank, const int max_rank, const ImageGrid<DIM>& grid, const Variables<T>& var, const std::string& vti_name)
    //-----------------------------------------------------------------------------------
    {
      set_filename();
      if (rank == 0) {
        open();
        write_header(time, grid, var);
        set_position_and_close();
      }
      if (mpi_running_)
        MPI_write_piece(vti_name, grid, rank);
      else
        write_piece(vti_name, grid);
      if (rank == max_rank) {
        write_footer();
        close();
      }
      inc_nwrite();
    }

    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Create a timestamp string for the header.
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

    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Build a Piece entry string for the given extent.
    template <int DIM>
    const std::string piece_string(const std::string& piece_file, const ImageGrid<DIM>& grid) const {
    //-----------------------------------------------------------------------------------
      std::ostringstream oss;
      oss << "    <Piece Extent=\"" << grid.extent_string(grid.piece_extent()) << "\" Source=\"" << piece_file << "\" />";
      return oss.str();
    }

    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Append a Piece entry (non-MPI).
    template <int DIM>
    void write_piece(const std::string& vti_filename, const ImageGrid<DIM>& grid) const {
    //-----------------------------------------------------------------------------------
      char piece[121] = {0}, piece_format[20] = {0};
      sprintf(piece_format, "%%-%ds\n", (int)(sizeof(piece))-1);  // -1 due to \n
      sprintf(piece, piece_format, piece_string(vti_filename, grid).c_str());
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

    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Append a Piece entry using MPI-IO.
    template <int DIM>
    void MPI_write_piece(const std::string& vti_filename, const ImageGrid<DIM>& grid, const int rank) {
    //-----------------------------------------------------------------------------------
      char piece[121] = {0}, piece_format[20] = {0};
      sprintf(piece_format, "%%-%ds\n", (int)(sizeof(piece))-2);
      sprintf(piece, piece_format, piece_string(vti_filename, grid).c_str());
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

    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Write the PVTI header and variable declarations.
    template <int DIM, typename T>
    void write_header(double time, const ImageGrid<DIM>& grid, const Variables<T>& var) {
    //-----------------------------------------------------------------------------------
      file_ << "<?xml version=\"1.0\"?>" << std::endl;
      file_ << "<!-- Created " << timestring() << " -->" << std::endl;
      file_ << "<!-- time = " << time << " s -->" << std::endl;
      file_ << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      file_ << "  <PImageData WholeExtent=\"" << grid.extent_string(grid.whole_extent()) << "\" Origin=\"" << grid.vec_string(grid.origin()) << "\" Spacing=\"" << grid.vec_string(grid.spacing()) << "\">" << std::endl;
      file_ << "    <PCellData>" << std::endl;
      for (auto& v : var.data())
        file_ << "      <PDataArray " << v.dataarray("name") << v.dataarray("type") << v.dataarray("dim") << "/>" << std::endl;
      file_ << "    </PCellData>" << std::endl;
    }

    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Close PImageData tag.
    void write_footer() {
    //-----------------------------------------------------------------------------------
      file_.open(path_+filename_, std::ios::app);
      file_ << "  </PImageData>" << std::endl;
    }

    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Create the pvti filename.
    void set_filename() {
    //-----------------------------------------------------------------------------------
      std::ostringstream ss;
      ss << std::setfill('0') << name_ << "_" << std::setw(7) << nwrite_ << extension_;
      filename_ = ss.str();
    }

    //                                   PVTI_file
    //-----------------------------------------------------------------------------------
    // Store insertion position for pieces and close the file.
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
  class OutfileImage {
    private:
    Variables<T> variables_;
    PVTI_file pvti_file_;
    VTI_file vti_file_;

    public:
    //                                  Outfile
    //-----------------------------------------------------------------------------------
    // Initialize per-file writers and variable storage.
    OutfileImage(std::string &_path, const std::string &_name, const unsigned int _offset) 
      : variables_(), pvti_file_(_path, _name), vti_file_(_path, _name, _offset) { }
    //-----------------------------------------------------------------------------------

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    // Write VTI/PVTI pair for the current step.
    template <int DIM>
    void write(const ImageGrid<DIM>& grid, double time, int rank, int max_rank)
    //-----------------------------------------------------------------------------------
    {
      vti_file_.write(rank, grid, variables_);
      pvti_file_.write(time, rank, max_rank, grid, variables_, vti_file_.path());
    }

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    // Update appended-data offsets after adding a variable.
    void update_offset() { 
    //-----------------------------------------------------------------------------------
      vti_file_.set_offset( variables_.back().update_offset( vti_file_.offset() ) ); 
    }

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    // Access the VTI writer.
    const VTI_file& vti_file() const { return vti_file_; }
    //-----------------------------------------------------------------------------------

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    // Access the PVTI writer.
    const PVTI_file& pvti_file() const { return pvti_file_; }
    //-----------------------------------------------------------------------------------

    //                                  Outfile
    //-----------------------------------------------------------------------------------
    // Access variable list.
    const Variables<T>& variables() const { return variables_; }
    //-----------------------------------------------------------------------------------
    // Access variable list (mutable).
    Variables<T>& variables() { return variables_; }
    //-----------------------------------------------------------------------------------

  };


  //=====================================================================================
  //
  //                               O U T P U T  I M A G E
  //
  //=====================================================================================
  // VTK::OutputImage manages one or more ImageData files and the variables attached to them.
  // The template parameter T determines the on-disk data type (Float32/Float64).
  template <int DIM=3, typename T=double>
  class OutputImage 
  {
    private:
    int format_ = BINARY;
    ImageGrid<DIM> grid_;
    std::string path_;
    std::vector<OutfileImage<T>> outfiles_;
    std::unordered_map<std::string, int> get_index_;
    int nwrite_ = 0;
    int rank_ = 0;
    int max_rank_ = 0;
    std::vector< std::unique_ptr<data_wrapper<T>> > wrappers_;
    //int dim_ = 3;

    private:
    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Verify that all points inside the bounding box are present (0-based index space).
    // This ensures the node list forms a full rectangular block required by ImageData.
    static void validate_full_block(const std::vector<int>& pos, const std::array<int, DIM>& min_pos, const std::array<int, DIM>& size, long long expected)
    //-----------------------------------------------------------------------------------
    {
      const size_t node_count = pos.size()/DIM;
      std::unordered_set<long long> keys;
      keys.reserve(node_count*2);
      long long stride = 1;
      std::array<long long, DIM> strides;
      for (int d=0; d<DIM; ++d) {
        strides[d] = stride;
        stride *= size[d];
      }

      // Shift into 0-based coordinates so the key space is dense.
      for (size_t n=0; n<node_count; ++n) {
        long long key = 0;
        for (int d=0; d<DIM; ++d) {
          key += (pos[n*DIM + d] - min_pos[d]) * strides[d];
        }
        keys.insert(key);
      }

      if (keys.size() != static_cast<size_t>(expected)) {
        std::cerr << "  ERROR in VTK::OutputImage: Node set does not form a full rectangular block" << std::endl;
        std::cerr << "  Expected " << expected << " nodes in the bounding box, got " << keys.size() << std::endl;
        std::cerr << "  Use unstructured output for irregular domains" << std::endl;
        util::safe_exit(EXIT_FAILURE);
      }
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Build global min/max using MPI.
    // This allows ImageData metadata to reflect the full domain across ranks.
    static void global_min_max(const std::array<int, DIM>& min_pos, const std::array<int, DIM>& max_pos, std::array<int, DIM>& global_min, std::array<int, DIM>& global_max, int num_procs)
    //-----------------------------------------------------------------------------------
    {
      global_min = min_pos;
      global_max = max_pos;
      int mpi_running = 0;
      MPI_Initialized(&mpi_running);
      if (mpi_running && num_procs > 1) {
        std::array<int, DIM> in_min = min_pos;
        std::array<int, DIM> in_max = max_pos;
        MPI_Allreduce(in_min.data(), global_min.data(), DIM, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(in_max.data(), global_max.data(), DIM, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      }
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Convert a global min/max to ImageData cells and origin.
    // Origin stores physical coordinates while extents stay 0-based in index space.
    static void origin_and_cells(const std::array<int, DIM>& global_min, const std::array<int, DIM>& global_max, std::array<double, DIM>& origin, std::array<int, DIM>& cells)
    //-----------------------------------------------------------------------------------
    {
      for (int d=0; d<DIM; ++d) {
        origin[d] = static_cast<double>(global_min[d]);
        cells[d] = global_max[d] - global_min[d] + 1;
      }
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Shift a local min/max into the index-space used by piece extents.
    // This aligns per-rank extents with a global 0-based index space.
    static void shift_min_max(const std::array<int, DIM>& min_pos, const std::array<int, DIM>& max_pos, const std::array<int, DIM>& origin, std::array<int, DIM>& min_shift, std::array<int, DIM>& max_shift)
    //-----------------------------------------------------------------------------------
    {
      for (int d=0; d<DIM; ++d) {
        min_shift[d] = min_pos[d] - origin[d];
        max_shift[d] = max_pos[d] - origin[d];
      }
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Derive extents from node positions and validate that the nodes form a full block.
    // grid.pos is treated as physical coordinates, so origin is set to the global min.
    // Piece extents are built in a 0-based index space by subtracting that origin.
    template <typename GridT>
    static void derive_extent_from_nodes(const GridT& grid, const std::vector<int>& nodes, std::array<int, DIM>& cells, std::array<int, 6>& piece_extent, std::array<double, DIM>& origin, int num_procs)
    //-----------------------------------------------------------------------------------
    {
      if (nodes.empty()) {
        std::cerr << "  ERROR in VTK::OutputImage: Node list is empty; cannot build ImageData grid" << std::endl;
        util::safe_exit(EXIT_FAILURE);
      }
      const std::vector<int> pos = grid.pos(nodes);
      const size_t node_count = nodes.size();
      if (pos.size() != node_count*DIM) {
        std::cerr << "  ERROR in VTK::OutputImage: grid.pos(nodes) size mismatch, expected "
                  << node_count*DIM << ", got " << pos.size() << std::endl;
        util::safe_exit(EXIT_FAILURE);
      }

      std::array<int, DIM> min_pos;
      std::array<int, DIM> max_pos;
      util::find_min_max<DIM>(pos, min_pos, max_pos);
      
      std::array<int, DIM> size;
      const long long expected = util::block_size<DIM>(min_pos, max_pos, size);
      // Ensure the local node list is a full rectangular block.
      validate_full_block(pos, min_pos, size, expected);
      std::array<int, DIM> global_min;
      std::array<int, DIM> global_max;
      // Promote local bounds to global bounds for metadata.
      global_min_max(min_pos, max_pos, global_min, global_max, num_procs);
      origin_and_cells(global_min, global_max, origin, cells);
      std::array<int, DIM> min_shift;
      std::array<int, DIM> max_shift;
      // Shift the local bounds into 0-based index space for piece extents.
      shift_min_max(min_pos, max_pos, global_min, min_shift, max_shift);
      util::build_piece_extent<DIM>(min_shift, max_shift, piece_extent);
    }

    public:
    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Construct an image output from grid positions and validate rectangular coverage.
    template <typename GridT>
    OutputImage(int format, const GridT& grid, const std::vector<int>& nodes, const std::string path="out", int rank=0, int num_procs=1)
    //-----------------------------------------------------------------------------------
      : format_(format),
        grid_(std::array<int, DIM>{0}, std::array<double, DIM>{0.0}, std::array<double, DIM>{1.0}, std::array<int, 6>{0, 0, 0, 0, 0, 0}),
        path_(path), outfiles_(), get_index_(), rank_(rank), max_rank_(num_procs-1), wrappers_()
    {
      std::array<int, DIM> cells;
      std::array<int, 6> piece_extent;
      std::array<double, DIM> origin;
      derive_extent_from_nodes(grid, nodes, cells, piece_extent, origin, num_procs);
      std::array<double, DIM> spacing;
      spacing.fill(1.0);
      grid_ = ImageGrid<DIM>(cells, origin, spacing, piece_extent);
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Construct an image output with optional piece extent (defaults to whole extent).
    OutputImage(int format, const std::array<int, DIM>& cells, const std::array<double, DIM>& origin, const std::array<double, DIM>& spacing, const std::array<int, 6>& piece_extent = {0, 0, 0, 0, 0, 0}, const std::string path="out", int rank=0, int num_procs=1) 
    //-----------------------------------------------------------------------------------
      : format_(format), grid_(cells, origin, spacing, piece_extent), path_(path), outfiles_(), get_index_(), rank_(rank), max_rank_(num_procs-1), wrappers_() { }


    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Add a new output file.
    OutfileImage<T>& add_file(const std::string& name)
    //-----------------------------------------------------------------------------------
    {
      outfiles_.emplace_back(path_, name, 0);
      get_index_[name] = outfiles_.size()-1;
      return outfiles_.back();
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Add a variable from a std::vector.
    void add_variable(const std::string& name, const std::vector<T>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      // wrappers_.emplace_back(new vec_wrapper<T>(data)); // c++11 version
      wrappers_.emplace_back(std::make_unique< vec_wrapper<T> >(data));
      add_variable_(name, index, length, offset);
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Enabled only when InT != T so mismatched input types go through the cast wrapper.
    // Add a variable from a std::vector with casting.
    template <typename InT, typename std::enable_if<!std::is_same<InT, T>::value, int>::type = 0>
    void add_variable(const std::string& name, const std::vector<InT>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      wrappers_.emplace_back(std::make_unique< vec_cast_wrapper<T, InT> >(data));
      add_variable_(name, index, length, offset);
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    const std::array<int, 6>& piece_extent() const { return grid_.piece_extent(); }
    //-----------------------------------------------------------------------------------

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    const std::array<double, 3>& origin() const { return grid_.origin(); }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Add a variable from a std::valarray.
    void add_variable(const std::string& name, const std::valarray<T>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      // wrappers_.emplace_back(new arr_wrapper<T>(data));
      wrappers_.emplace_back(std::make_unique< arr_wrapper<T> >(data));
      add_variable_(name, index, length, offset);
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Enabled only when InT != T so mismatched input types go through the cast wrapper.
    // Add a variable from a std::valarray with casting.
    template <typename InT, typename std::enable_if<!std::is_same<InT, T>::value, int>::type = 0>
    void add_variable(const std::string& name, const std::valarray<InT>& data, const std::vector<int>& index=std::vector<int>(), int length=0, int offset=0)
    //-----------------------------------------------------------------------------------
    {
      wrappers_.emplace_back(std::make_unique< arr_cast_wrapper<T, InT> >(data));
      add_variable_(name, index, length, offset);
    }

    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Write all files for the current time step.
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
      std::cout << "VTK::OutputImage::write() duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << std::endl;
#endif
    }
  
  private:
    //                                     OutputImage
    //-----------------------------------------------------------------------------------
    // Register a variable and update offsets.
    void add_variable_(const std::string& name, const std::vector<int>& index, int length, int offset)
    //-----------------------------------------------------------------------------------
    {
      int size = 1;
      if (index.size() > 0) {
        size = index.size();
      } else {
        size = (*wrappers_.back()).size();
      }
      const int piece_cells = grid_.num_piece_cells();
      if (piece_cells <= 0 || (size % piece_cells) != 0) {
        std::cerr << "  ERROR in VTK::OutputImage::add_variable(" << name << "):" << std::endl;
        std::cerr << "  Variable size (" << size << ") is not compatible with the number of cells in this piece (" << piece_cells << ")" << std::endl;
        std::cerr << "  This usually means the ImageData piece extents do not match the node ordering used to build indices." << std::endl;
        util::safe_exit(EXIT_FAILURE);
      }
      int dim = size/grid_.num_piece_cells();
      if (dim != 1 and dim != DIM) {
        std::cerr << "  ERROR in VTK::OutputImage::add_variable(" << name << ", " << dim << "):" << std::endl;
        std::cerr << "  Wrong dimension: Expected " << DIM << " or 1, got " << dim << std::endl;
        util::safe_exit(EXIT_FAILURE);
      }
      if (outfiles_.empty()) {
        std::cerr << "  ERROR in VTK::OutputImage::add_variable(" << name << ", " << dim << "):" << std::endl;
        std::cerr << "  Add an output file using 'add_file(name)' before adding variables" << std::endl;
        util::safe_exit(EXIT_FAILURE);
      }
      outfiles_.back().variables().add(name, *wrappers_.back(), format_, dim, index, length, offset);
      outfiles_.back().update_offset();
    }

  };

}

#endif /* SRC_VTK_IMAGE_H_ */
