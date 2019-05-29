#include "Output.h"

//////////////////////////////////
//                              //
//   Variable class functions   //
//                              //
//////////////////////////////////


//--------------------------------------------
//
//--------------------------------------------
void Variable::set_data_array() {
  std::ostringstream oss;
  oss << "Name=\"" << name << "\" type=\"" << datatype << "\"";
  if (dim>1) {
    oss << " NumberOfComponents=\"3\"";
  }
  data_array = oss.str();
}


//////////////////////////////////
//                              //
//   File class functions   //
//                              //
//////////////////////////////////

//--------------------------------------------
//
//--------------------------------------------
void File::make_dir(std::string &dir) {
  // check that dir has a '/' at the end
  if (dir.back() != '/') {
    dir += '/';
  }
  struct stat st = {0};
  if (stat(dir.c_str(), &st) == -1) {
#ifdef _WIN32
      _mkdir(folder.c_str());
#else
      mkdir(dir.c_str(), 0700);
#endif
  }
}



//--------------------------------------------
//
//--------------------------------------------
void VTU_file::write_header(int num_points, int num_cells) {
  file_ << "<?xml version=\"1.0\"?>"                                                           << std::endl;
  file_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"   << std::endl;
  file_ << "  <UnstructuredGrid>"                                                              << std::endl;
  file_ << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells << "\">" << std::endl;
}

//--------------------------------------------
//
//--------------------------------------------
void VTU_file::write_header(VTKGrid& grid) {
  file_ << "<?xml version=\"1.0\"?>"                                                           << std::endl;
  file_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"   << std::endl;
  file_ << "  <UnstructuredGrid>"                                                              << std::endl;
  file_ << "    <Piece NumberOfPoints=\"" << grid.get_num_points() << "\" NumberOfCells=\"" << grid.get_num_cells() << "\">" << std::endl;
}

//--------------------------------------------
//
//--------------------------------------------
void VTU_file::write_grid(std::vector<double> &points, std::vector<int> &connectivity, int num_cells, int num_cell_pts) {
  file_ << "      <Points>" << std::endl;
  file_ << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << vector_as_string_list(points) << "</DataArray>" << std::endl;
  file_ << "      </Points>" << std::endl;
  file_ << "      <Cells>" << std::endl;
  file_ << "        <DataArray type=\"UInt32\"  Name=\"connectivity\" format=\"ascii\">" << vector_as_string_list(connectivity) << "</DataArray>" << std::endl;
  file_ << "        <DataArray type=\"UInt32\"  Name=\"offsets\" format=\"ascii\">" << incremental_list<int>(num_cell_pts, connectivity.size(), num_cell_pts) << "</DataArray>" << std::endl;
  file_ << "        <DataArray type=\"UInt32\" Name=\"types\" format=\"ascii\">" << repeated_list<int>(cell_type_, num_cells) << "</DataArray>" << std::endl;
  file_ << "      </Cells>" << std::endl;
}

void VTU_file::write_grid(VTKGrid& grid) {
  file_ << "      <Points>" << std::endl;
  file_ << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << grid.get_point_list() << "</DataArray>" << std::endl;
  file_ << "      </Points>" << std::endl;
  file_ << "      <Cells>" << std::endl;
  file_ << "        <DataArray type=\"UInt32\"  Name=\"connectivity\" format=\"ascii\">" << grid.get_connectivity_list() << "</DataArray>" << std::endl;
  file_ << "        <DataArray type=\"UInt32\"  Name=\"offsets\" format=\"ascii\">" << grid.get_offset_list() << "</DataArray>" << std::endl;
  file_ << "        <DataArray type=\"UInt32\" Name=\"types\" format=\"ascii\">" << grid.get_types_list() << "</DataArray>" << std::endl;
  file_ << "      </Cells>" << std::endl;
}

//--------------------------------------------
//
//--------------------------------------------
void VTU_file::write_data() {
  file_ << "      <CellData Scalars=\"" << scalar_string_ << "\" Vectors=\"" << vector_string_ << "\">" << std::endl;
  for (const auto& var : variables_) {
    file_ << "        <DataArray " << var.data_array << " format=\"ascii\">";
    for (const auto val:var.data_iter_)
      file_ << " " << *val;
    file_ << "</DataArray>" << std::endl;
  }
  file_ << "      </CellData>" << std::endl;
}


//--------------------------------------------
//
//--------------------------------------------
void VTU_file::write_footer_and_close() {
  file_ << "    </Piece>" << std::endl;
  file_ << "  </UnstructuredGrid>" << std::endl;
  file_ << "</VTKFile>" << std::endl;
  file_.close();
}


//--------------------------------------------
//
//--------------------------------------------
void VTU_file::set_filename_and_open(const int rank) {
  std::ostringstream ss;
  ss << std::setfill('0') << std::setw(4) << rank << "_" << name_ << "_" << std::setw(7) << nwrite_ << extension_;
  filename_ = ss.str();
  file_.open(path_+filename_, std::ios::out | std::ios::binary);
}

//--------------------------------------------
//
//--------------------------------------------
const std::string VTU_file::get_piece_string() const {
  std::ostringstream oss;
  oss << "    <Piece Source=\"" << folders_.back() + filename_ << "\" />";
  return oss.str();
}

//--------------------------------------------
//
//--------------------------------------------
void VTU_file::set_cell_data_string() {
  // make lists of scalars and vectors
  std::string sep_s = "", sep_v = "";
  std::ostringstream scalars, vectors, oss;
  for (const auto& v : variables_) {
    if (v.dim>1) {
      vectors << sep_v << v.name;
      sep_v = ", ";
    } else {
      scalars << sep_s << v.name;
      sep_s = ", ";
    }
  }
  oss << "Scalars=\"" << scalars.str() << "\" Vectors=\"" << vectors.str() << "\">";
  cell_data_string_ = oss.str();
}

//--------------------------------------------
//
//--------------------------------------------
void VTU_file::update_cell_data_string(Variable &var) {
  // make lists of scalars and vectors
  std::string sep = "";
  if (var.dim>1) {
    // vector
    if (num_vector_>0) {
      sep = " ,";
    }
    vector_string_ += sep + var.name;
    ++num_vector_;
  } else {
    // scalar
    if (num_scalar_>0) {
      sep = " ,";
    }
    scalar_string_ += sep + var.name;
    ++num_scalar_;
  }
}



///////////////////////////////////
//                               //
//   PVTU_file class functions   //
//                               //
///////////////////////////////////

//--------------------------------------------
//
//--------------------------------------------
std::string PVTU_file::get_timestring() {
  time_t t = time(NULL);
  tm tm = *localtime(&t);
  std::ostringstream ss;
  ss << std::setfill('0') << std::setw(2) << tm.tm_mday << "." << std::setw(2) << tm.tm_mon+1 << "."
      << std::setw(2) << tm.tm_year+1900 << " " << std::setw(2) << tm.tm_hour <<  ":" << std::setw(2)
  << tm.tm_min << ":" << std::setw(2) << tm.tm_sec;
  return ss.str();
};


//--------------------------------------------
//
//--------------------------------------------
void PVTU_file::MPI_write_piece(const std::string& piece_string, const int rank) {

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
//
//--------------------------------------------
void PVTU_file::write_header(const double time, const VTU_file &vtu) {
  file_.precision(5);
  file_ << "<?xml version=\"1.0\"?>"                                                   << std::endl;
  file_ << "<!-- Created " << get_timestring() << " -->"                               << std::endl;
  file_ << "<!-- time = " << time << " s -->"                                          << std::endl;
  file_ << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  file_ << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
  file_ << "    <PPoints>" << std::endl;
  file_ << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\" />" << std::endl;
  file_ << "    </PPoints>" << std::endl;
  file_ << "    <PCellData>" << std::endl;
  file_ << "      <PDataArray type=\"Int32\" Name=\"connectivity\" />" << std::endl;
  file_ << "      <PDataArray type=\"Int32\" Name=\"offsets\" />" << std::endl;
  file_ << "      <PDataArray type=\"Int32\" Name=\"types\" />" << std::endl;
  file_ << "    </PCellData>" << std::endl;
  file_ << "    <PCellData>"                                                          << std::endl;
  for (const auto& v : vtu.get_variables())
    file_ << "      <PDataArray " << v.data_array << " />"                             << std::endl;
  file_ << "    </PCellData>"                                                          << std::endl;
}


//--------------------------------------------
//
//--------------------------------------------
void PVTU_file::write_footer_and_close() {
  file_.open(path_+filename_, std::ios::app);
  file_ << "  </PUnstructuredGrid>" << std::endl;
  file_ << "</VTKFile>" << std::endl;
  file_.close();
}


//--------------------------------------------
//
//--------------------------------------------
void PVTU_file::set_filename() {
  //filename = get_filename(nwrite);
  std::ostringstream ss;
  ss << std::setfill('0') << name_ << "_" << std::setw(7) << nwrite_ << extension_;
  filename_ = ss.str();
}


//--------------------------------------------
//
//--------------------------------------------
void PVTU_file::open() {
  file_.open(path_+filename_, std::ios::out | std::ios::binary);
}


//--------------------------------------------
//
//--------------------------------------------
void PVTU_file::set_position_and_close() {
  file_position = file_.tellp();
  file_.close();
}


////////////////////////////////
//                            //
//   Output class functions   //
//                            //
////////////////////////////////

//--------------------------------------------
//
//--------------------------------------------
Outfile& Output::add_file(const std::string &_name) {
  //file.emplace_back(path, _name, geo_, mpi_);
  //file.emplace_back(path, _name, geo_, num_ghosts_);
  outfiles_.emplace_back(path, _name);
  get_index[_name] = outfiles_.size()-1;
  return outfiles_.back();
};

void Output::set_points_and_connectivity(std::vector<std::vector<int>>& node_pos) {
  // set 2D or 3D cell
  std::vector<std::vector<int>> cell_points = (dim_.size()>2)? cell_points_3D_ : cell_points_2D_;
  num_cell_points_ = cell_points.size();
  connectivity_.reserve(cell_points.size()*node_pos.size());

  //std::vector<std::vector<int>> nodes(node_pos.size(), std::vector<int>(dim_.size()));
  //for (auto i=0; i<node_pos.size(); ++i) {
  //  nodes[i] = std::vector<int>(node_pos[i], node_pos[i]+dim_.size());
  //}
  std::vector<int> point_index(prod(dim_+1), -1);
  //for (int i=0; i<num_nodes; ++i) {
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
  points_ = points_ - 0.5*cell_size_;
}




