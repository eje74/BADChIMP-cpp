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

//--------------------------------------------
//
//--------------------------------------------
void Variable::set_write_nodes(const std::vector<int> &_node_mode_to_write) {
  write_node.reset();  // set all to zero
  // loop through possible list of node-modes to include
  for (auto & mode : _node_mode_to_write) {
    write_node.set(mode, 1);
  }
}

//--------------------------------------------
//
//--------------------------------------------
void Variable::set_nbytes(const std::vector<int> &n) {
  //int ndata = (n[0]-2)*(n[1]-2)*std::max(n[2]-2,1);  // max needed for 2d-runs
  //nbytes = (unsigned int) sizeof(OUTPUT_DTYPE)*ndata*dim;
  nbytes = (unsigned int) sizeof(OUTPUT_DTYPE)*prod(n)*dim;
}

//--------------------------------------------
//
//--------------------------------------------
void Variable::set_offset(const Variable &prev_var) {
  offset = prev_var.offset + prev_var.nbytes + sizeof(unsigned int);
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


//////////////////////////////////
//                              //
//   VTI_file class functions   //
//                              //
//////////////////////////////////


//--------------------------------------------
//
//--------------------------------------------
void VTI_file::set_extent(const Geo &geo, const int num_ghosts) {
  std::ostringstream ss;
  std::vector<int> lb = geo.get_lower_bounds() - num_ghosts;
  std::vector<int> ub = geo.get_upper_bounds();
  if (lb.size()<3) {
    lb.push_back(0);
    ub.push_back(1);
  }
  //ss  << lb[0]-1 << " " << ub[0]  << " " << lb[1]-1  << " " << ub[1]  << " " << std::max(geo.lb_[2]-1,0) << " " << geo.ub_[2];
  ss  << lb[0] << " " << ub[0]  << " " << lb[1]  << " " << ub[1]  << " " << lb[2] << " " << ub[2];
  extent = ss.str();
}

//--------------------------------------------
//
//--------------------------------------------
void VTI_file::write_data(int*** labels, const int num_ghosts) {
  //if (buffer)
  //  update_buffer(sys);

  // opening tag
  file_ << "  <AppendedData encoding=\"raw\">" << std::endl;
  file_ << "_";

  //std::cout << n[0] << ", " << n[1] << ", " << n[2] << std::endl;

  int ng = num_ghosts;
  OUTPUT_DTYPE data;
  //std::vector<int> &n = mpi_.n_;


  //std::cout << std::endl << "WRITE_DATA:" << std::endl << std::endl;

  for (const auto& var : variables_) {
    // precede data with total number of bytes
    file_.write((char*)&(var.nbytes), sizeof(unsigned int));
    //std::vector<int> &n = var.n;
    // node loop
    //int z = (n.size()>2) ? ng : 0;  // 2D:z=0, 3D:z=num_ghosts to skip periodic/ghost rim
    int z=0, nz=0;
    if (n.size()>2) {
      z = ng;
      nz = n[2];
    }
    do {
      //std::cout << z << std::endl;
      for(int y=ng; y<n[1]-ng; ++y) {
        for(int x=ng; x<n[0]-ng; ++x) {
          //std::cout << std::vector<int>({x,y,z,nz,ng}) << std::flush;
          int nn = labels[z][y][x];
          //int nn = x + y*n[0] + z*n[0]*n[1];   // could be part of Geo-class?
          bool write_node = true; //var.write_node.test(node->mode[nn]);
          // stride loop
          for(int dim=0; dim<var.dim; ++dim) {
            //if ( (dim<n.size()) && write_node ) {
            if ( write_node ) {
              //int element = var[i].data_step*nn + dim;
              data = var.get_data<OUTPUT_DTYPE>(nn, dim);
              //std::cout << " " << nn << " ";
              //std::cout << std::setw(6) << std::setprecision(3) << data;
              //std::cout << data;
            } else {
              data = 0.0;
              //std::cout << std::setprecision(3) << 0;
            }
            file_.write((char*)&data, sizeof(OUTPUT_DTYPE));
          }

        }
      }
      ++z;
      //std::cout << z << "<" << nz << "-" << ng << std::endl;
    } while (z<nz-ng);
    //std::cout << "WHILE: " << std::endl;
  }
  // end tag
  file_ << "  </AppendedData>" << std::endl;
  //std::cout << std::endl;
}


//--------------------------------------------
//
//--------------------------------------------
void VTI_file::write_data() {

  //if (buffer)
  //  update_buffer(sys);

  // opening tag
  file_ << "  <AppendedData encoding=\"raw\">" << std::endl;
  file_ << "_";


  OUTPUT_DTYPE data;
  for (const auto& var : variables_) {
    // precede data with total number of bytes
    file_.write((char*)&(var.nbytes), sizeof(unsigned int));
    for(int nn=0; nn<prod(n); ++nn) {
      bool write_node = true; //var.write_node.test(node->mode[nn]);
      // dimension loop
      for(int dim=0; dim<var.dim; ++dim) {
        //if ( (dim<n.size()) && write_node ) {
        if (write_node) {
          data = var.get_data<OUTPUT_DTYPE>(nn, dim);
        } else {
          data = 0.0;
        }
        file_.write((char*)&data, sizeof(OUTPUT_DTYPE));
      }

    }
  }
  // end tag
  file_ << "  </AppendedData>" << std::endl;
}


//--------------------------------------------
//
//--------------------------------------------
void VTI_file::write_header() {
  file_ << "<?xml version=\"1.0\"?>"                                                           << std::endl;
  file_ << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"          << std::endl;
  file_ << "  <ImageData WholeExtent=\"" << extent << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
  file_ << "    <Piece Extent=\"" << extent << "\">"                                           << std::endl;
  file_ << "      <CellData " << cell_data_string_  << ">"                                     << std::endl;
  for (const auto& v : variables_)
    file_ << "        <DataArray " << v.data_array << " format=\"appended\" offset=\"" << v.offset << "\"/>" << std::endl;
  file_ << "      </CellData>"                                                                 << std::endl;
  file_ << "    </Piece>"                                                                      << std::endl;
  file_ << "  </ImageData>"                                                                    << std::endl;
}


//--------------------------------------------
//
//--------------------------------------------
void VTI_file::write_footer_and_close() {
  file_ << "</VTKFile>"        << std::endl;
  file_.close();
}

//--------------------------------------------
//
//--------------------------------------------
void VTI_file::set_filename_and_open(const int rank) {
  std::ostringstream ss;
  ss << std::setfill('0') << std::setw(4) << rank << "_" << name_ << "_" << std::setw(7) << nwrite_ << extension_;
  filename_ = ss.str();
  file_.open(path_+filename_, std::ios::out | std::ios::binary);
}

//--------------------------------------------
//
//--------------------------------------------
const std::string VTI_file::get_piece_extent_string() const {
  std::ostringstream oss;
  oss << "    <Piece Extent=\"" << extent << "\" Source=\"" << folders_.back() + filename_ << "\" />";
  return oss.str();
}

//--------------------------------------------
//
//--------------------------------------------
void VTI_file::set_cell_data_string() {
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
  oss << "Scalars=\"" << scalars.str() << "\" Vectors=\"" << vectors.str() << "\"";
  cell_data_string_ = oss.str();
}

//--------------------------------------------
//
//--------------------------------------------
//void VTI_file::add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
//    const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides,
//    const int num_ghosts)
//{
//  for (std::size_t i=0; i<dims.size(); ++i) {
//    variables_.emplace_back(names[i], data_ptrs[i], datasizes[i], dims[i], data_strides[i], n-2*num_ghosts);
//    std::size_t end = variables_.size()-1;
//    if (end>0)
//      // set the offset in bytes between variables in the same .vti-file
//      variables_[end].set_offset(variables_[end-1]);
//  }
//  set_cell_data_string();
//}


///////////////////////////////////
//                               //
//   PVTI_file class functions   //
//                               //
///////////////////////////////////

//--------------------------------------------
//
//--------------------------------------------
std::string PVTI_file::get_timestring() {
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
void PVTI_file::set_whole_extent(const Geo& geo, const int num_ghosts) {
  std::ostringstream ss;
  //std::vector<int> Nng = geo.global_.get_size() - 2*num_ghosts_;
  std::vector<int> Nng = geo.get_N() - 2*num_ghosts;
  if (Nng.size()<3) {
    Nng.push_back(1);
  }
  ss << "0 " << Nng[0] << " " << "0 " << Nng[1] << " " << "0 " << Nng[2];
  whole_extent = ss.str();
}


//--------------------------------------------
//
//--------------------------------------------
//void PVTI_file::MPI_write_piece(const VTI_file &vti) {
void PVTI_file::MPI_write_piece(const std::string& piece_extent_string, const int rank) {
  //std::ostringstream piece_tag;
  //piece_tag << "    <Piece Extent=\"" << vti.extent << "\" Source=\"" << vti.folders_.back() + vti.filename_ << "\" />";
  //piece_tag << "    <Piece Extent=\"" << vti.extent << "\" Source=\"" << vti.get_filename_with_path() << "\" />";
  //std::string piece_tag = vti.get_piece_tag();

  char piece[101] = {0}, piece_format[20] = {0};
  sprintf(piece_format, "%%-%ds\n", (int)(sizeof(piece))-2);
  sprintf(piece, piece_format, piece_extent_string.c_str());

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
void PVTI_file::write_header(const double time, const VTI_file &vti) {
  file_.precision(5);
  file_ << "<?xml version=\"1.0\"?>"                                                   << std::endl;
  file_ << "<!-- Created " << get_timestring() << " -->"                               << std::endl;
  file_ << "<!-- time = " << time << " s -->"                                          << std::endl;
  file_ << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  file_ << "  <PImageData WholeExtent=\"" << whole_extent << "\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
  file_ << "    <PCellData " << vti.get_cell_data_string()                             << std::endl;
  for (const auto& v : vti.get_variables())
    file_ << "      <PDataArray " << v.data_array << " />"                             << std::endl;
  file_ << "    </PCellData>"                                                          << std::endl;
}


//--------------------------------------------
//
//--------------------------------------------
void PVTI_file::write_footer_and_close() {
  file_.open(path_+filename_, std::ios::app);
  file_ << "  </PImageData>" << std::endl;
  file_ << "</VTKFile>" << std::endl;
  file_.close();
}


//--------------------------------------------
//
//--------------------------------------------
void PVTI_file::set_filename() {
  //filename = get_filename(nwrite);
  std::ostringstream ss;
  ss << std::setfill('0') << name_ << "_" << std::setw(7) << nwrite_ << extension_;
  filename_ = ss.str();
}


//--------------------------------------------
//
//--------------------------------------------
void PVTI_file::open() {
  file_.open(path_+filename_, std::ios::out | std::ios::binary);
}


//--------------------------------------------
//
//--------------------------------------------
void PVTI_file::set_position_and_close() {
  file_position = file_.tellp();
  file_.close();
}


//////////////////////////////////
//                              //
//   VTU_file class functions   //
//                              //
//////////////////////////////////


//--------------------------------------------
//
//--------------------------------------------
//void VTU_file::set_extent(const Geo &geo, const int num_ghosts) {
//  std::ostringstream ss;
//  std::vector<int> lb = geo.get_lower_bounds() - num_ghosts;
//  std::vector<int> ub = geo.get_upper_bounds();
//  if (lb.size()<3) {
//    lb.push_back(0);
//    ub.push_back(1);
//  }
//  //ss  << lb[0]-1 << " " << ub[0]  << " " << lb[1]-1  << " " << ub[1]  << " " << std::max(geo.lb_[2]-1,0) << " " << geo.ub_[2];
//  ss  << lb[0] << " " << ub[0]  << " " << lb[1]  << " " << ub[1]  << " " << lb[2] << " " << ub[2];
//  extent = ss.str();
//}

//--------------------------------------------
//
//--------------------------------------------
//void VTU_file::write_data(int*** labels, const int num_ghosts) {
//  //if (buffer)
//  //  update_buffer(sys);
//
//  // opening tag
//  file_ << "  <AppendedData encoding=\"raw\">" << std::endl;
//  file_ << "_";
//
//  //std::cout << n[0] << ", " << n[1] << ", " << n[2] << std::endl;
//
//  int ng = num_ghosts;
//  OUTPUT_DTYPE data;
//  //std::vector<int> &n = mpi_.n_;
//
//
//  //std::cout << std::endl << "WRITE_DATA:" << std::endl << std::endl;
//
//  for (const auto& var : variables_) {
//    // precede data with total number of bytes
//    file_.write((char*)&(var.nbytes), sizeof(unsigned int));
//    //std::vector<int> &n = var.n;
//    // node loop
//    //int z = (n.size()>2) ? ng : 0;  // 2D:z=0, 3D:z=num_ghosts to skip periodic/ghost rim
//    int z=0, nz=0;
//    if (n.size()>2) {
//      z = ng;
//      nz = n[2];
//    }
//    do {
//      //std::cout << z << std::endl;
//      for(int y=ng; y<n[1]-ng; ++y) {
//        for(int x=ng; x<n[0]-ng; ++x) {
//          //std::cout << std::vector<int>({x,y,z,nz,ng}) << std::flush;
//          int nn = labels[z][y][x];
//          //int nn = x + y*n[0] + z*n[0]*n[1];   // could be part of Geo-class?
//          bool write_node = true; //var.write_node.test(node->mode[nn]);
//          // stride loop
//          for(int dim=0; dim<var.dim; ++dim) {
//            //if ( (dim<n.size()) && write_node ) {
//            if ( write_node ) {
//              //int element = var[i].data_step*nn + dim;
//              data = var.get_data<OUTPUT_DTYPE>(nn, dim);
//              //std::cout << " " << nn << " ";
//              //std::cout << std::setw(6) << std::setprecision(3) << data;
//              //std::cout << data;
//            } else {
//              data = 0.0;
//              //std::cout << std::setprecision(3) << 0;
//            }
//            file_.write((char*)&data, sizeof(OUTPUT_DTYPE));
//          }
//
//        }
//      }
//      ++z;
//      //std::cout << z << "<" << nz << "-" << ng << std::endl;
//    } while (z<nz-ng);
//    //std::cout << "WHILE: " << std::endl;
//  }
//  // end tag
//  file_ << "  </AppendedData>" << std::endl;
//  //std::cout << std::endl;
//}


//--------------------------------------------
//
//--------------------------------------------
//void VTU_file::write_data() {
//
//  //if (buffer)
//  //  update_buffer(sys);
//
//  // opening tag
//  file_ << "  <AppendedData encoding=\"raw\">" << std::endl;
//  file_ << "_";
//
//
//  OUTPUT_DTYPE data;
//  for (const auto& var : variables_) {
//    // precede data with total number of bytes
//    file_.write((char*)&(var.nbytes), sizeof(unsigned int));
//    for(int nn=0; nn<prod(n); ++nn) {
//      bool write_node = true; //var.write_node.test(node->mode[nn]);
//      // dimension loop
//      for(int dim=0; dim<var.dim; ++dim) {
//        //if ( (dim<n.size()) && write_node ) {
//        if (write_node) {
//          data = var.get_data<OUTPUT_DTYPE>(nn, dim);
//        } else {
//          data = 0.0;
//        }
//        file_.write((char*)&data, sizeof(OUTPUT_DTYPE));
//      }
//
//    }
//  }
//  // end tag
//  file_ << "  </AppendedData>" << std::endl;
//}

//--------------------------------------------
//
//--------------------------------------------
void VTU_file::write_header(int num_points, int num_cells) {
  file_ << "<?xml version=\"1.0\"?>"                                                           << std::endl;
  file_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"   << std::endl;
  file_ << "  <UnstructuredGrid>"                                                              << std::endl;
  file_ << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells << "\">" << std::endl;
//  file_ << "      <CellData " << cell_data_string_  << ">"                                     << std::endl;
//  for (const auto& v : variables_)
//    file_ << "        <DataArray " << v.data_array << " format=\"appended\" offset=\"" << v.offset << "\"/>" << std::endl;
//  file_ << "      </CellData>"                                                                 << std::endl;
//  file_ << "      <Points>"                                                                    << std::endl;
//  file_ << "        <DataArray NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"/>"    << std::endl;
//  file_ << "      </Points>"                                                                   << std::endl;
//  file_ << "      <Cells>"                                                                    << std::endl;
//  file_ << "        <DataArray type=\"Int32\"  Name=\"connectivity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"/>" << std::endl;
//  file_ << "        <DataArray type=\"Int32\"  Name=\"offsets\"      NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"/>" << std::endl;
//  file_ << "        <DataArray type=\"UInt32\" Name=\"type\"         NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"/>" << std::endl;
//  file_ << "      </Cells>"                                                                   << std::endl;
//  file_ << "    </Piece>"                                                                      << std::endl;
//  file_ << "  </UnstructuredGrid>"                                                             << std::endl;
}


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

void VTU_file::write_data() {
  file_ << "      <CellData Scalars=\"" << scalar_string_ << "\" Vectors=\"" << vector_string_ << "\">" << std::endl;
  for (const auto& var : variables_) {
    //file_ << "        <DataArray " << v.data_array << " format=\"appended\" offset=\"" << v.offset << "\"/>" << std::endl;
    //file_ << "        <DataArray " << v.data_array << " format=\"ascii\">" << vector_as_string_list(*(v.data_ptr_)) << std::endl;
    file_ << "        <DataArray " << var.data_array << " format=\"ascii\">";
    for (const auto val:var.data_iter_)
      file_ << " " << *val;
    file_ << "</DataArray>" << std::endl;
  }
  //file_ << "        <DataArray " << v.data_array << " format=\"ascii\">" << 0.1234 << "</DataArray>" << std::endl; //(*(v.data_ptr_))[0] << std::endl;
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

//--------------------------------------------
//
//--------------------------------------------
//void VTU_file::add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
//    const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides)
//{
//  for (std::size_t i=0; i<dims.size(); ++i) {
//    //variables_.emplace_back(names[i], data_ptrs[i], datasizes[i], dims[i], data_strides[i], n-2*num_ghosts);
//    variables_.emplace_back(names[i], data_ptrs[i], datasizes[i], dims[i], data_strides[i]);
//    std::size_t end = variables_.size()-1;
//    if (end>0)
//      // set the offset in bytes between variables in the same .vti-file
//      variables_[end].set_offset(variables_[end-1]);
//  }
//  set_cell_data_string();
//}

//--------------------------------------------
//
//--------------------------------------------
//void VTU_file::add_variable(const std::string &name, const std::vector<double> *data_ptr,
//    const int datasize, const int dim, const int data_stride)
//{
//  variables_.emplace_back(name, data_ptr, datasize, dim, data_stride);
//  std::size_t end = variables_.size()-1;
//    if (end>0)
//      // set the offset in bytes between variables in the same .vti-file
//      variables_[end].set_offset(variables_[end-1]);
//    set_cell_data_string();
//}

//--------------------------------------------
//
//--------------------------------------------
//void VTU_file::add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
//    const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides,
//    const int num_ghosts)
//{
//  for (std::size_t i=0; i<dims.size(); ++i) {
//    variables_.emplace_back(names[i], data_ptrs[i], datasizes[i], dims[i], data_strides[i], n-2*num_ghosts);
//    std::size_t end = variables_.size()-1;
//    if (end>0)
//      // set the offset in bytes between variables in the same .vti-file
//      variables_[end].set_offset(variables_[end-1]);
//  }
//  set_cell_data_string();
//}


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
//void PVTU_file::set_whole_extent(const Geo& geo, const int num_ghosts) {
//  std::ostringstream ss;
//  //std::vector<int> Nng = geo.global_.get_size() - 2*num_ghosts_;
//  std::vector<int> Nng = geo.get_N() - 2*num_ghosts;
//  if (Nng.size()<3) {
//    Nng.push_back(1);
//  }
//  ss << "0 " << Nng[0] << " " << "0 " << Nng[1] << " " << "0 " << Nng[2];
//  whole_extent = ss.str();
//}


//--------------------------------------------
//
//--------------------------------------------
//void PVTI_file::MPI_write_piece(const VTI_file &vti) {
void PVTU_file::MPI_write_piece(const std::string& piece_string, const int rank) {
  //std::ostringstream piece_tag;
  //piece_tag << "    <Piece Extent=\"" << vti.extent << "\" Source=\"" << vti.folders_.back() + vti.filename_ << "\" />";
  //piece_tag << "    <Piece Extent=\"" << vti.extent << "\" Source=\"" << vti.get_filename_with_path() << "\" />";
  //std::string piece_tag = vti.get_piece_tag();

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


/////////////////////////////////
//                             //
//   Outfile class functions   //
//                             //
/////////////////////////////////


//--------------------------------------------
//
//--------------------------------------------
//void Outfile::add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
//    const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides)
//{
//  if ( (names.size()!=dims.size()) || (dims.size()!=data_ptrs.size()) ) {
//    std::cerr << "ERROR in Outfile::add_variables: the sizes of the input arrays must equal!" << std::endl;
//    exit(-1);
//  }
//
//  //vti_file_.add_variables(names, data_ptrs, datasizes, dims, data_strides, num_ghosts_);
//  vtu_file_.add_variables(names, data_ptrs, datasizes, dims, data_strides);
//}



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


//--------------------------------------------
//
//--------------------------------------------
//void Output::write_geo_file(const std::string &name,  Node *node, System *sys)
//{
//  static int ncall = 0;
//  ncall++;
//  if (ncall==1) {
//    geo.init(name, 1, node);
//    int stride = 1, data_step = 1;
//    // print all nodes
//    // NB the periodic nodes are not written anyway since the size is N-2
//    std::vector<int> include_modes;
//    for (int i=0; i<DISCARD; ++i)
//      include_modes.push_back(i);
//    geo.add_variable("mode", sizeof(char), stride, node->mode, data_step, sys, include_modes);
//
//  }
//  geo.name = name;
//  geo.write(sys);
//}


//--------------------------------------------
//
//--------------------------------------------
//void Outfile::update_buffer(System *sys)
//{
//  Links *links = node->links;
//  Links::Link *link;
//
//  // zero out
//  for (int n = 0; n < sys->max_n_tot; ++n)
//    buffer[n] = 0.0;
//  for (int nl = 0; nl < links->nlinks; ++nl) {
//    link = &links->list[nl];
//    if (!link->inert) { // skip inert nodes
//      buffer[link->fluid] = (OUTPUT_DTYPE) -link->pH;
//    }
//  }
//}


