#include "Output.h"


//--------------------------------------------
//
//--------------------------------------------
OUTPUT_DTYPE Variable::get_data(int nn, int dim) const {
  int element = data_stride*nn + dim;
  if (datasize==1) {
    char *char_pointer = (char *) data_pointer;
    int ival = char_pointer[element];
    return (OUTPUT_DTYPE) ival;
  } else if (datasize==8) {
    double *double_pointer = (double *) data_pointer;
    return (OUTPUT_DTYPE) double_pointer[element];
  } else {
    printf("\nERROR in Variable::get_data: unrecognized datasize: %zu\n", datasize);
    MPI_Finalize();
    return -1;
  }
}

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
  int ndata = (n[0]-2)*(n[1]-2)*std::max(n[2]-2,1);  // max needed for 2d-runs
  nbytes = (unsigned int) sizeof(OUTPUT_DTYPE)*ndata*dim;
}

//--------------------------------------------
//
//--------------------------------------------
void Variable::set_offset(const Variable &prev_var) {
  offset = prev_var.offset + prev_var.nbytes + sizeof(unsigned int);
}

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
void VTI_file::set_extent(const MPI &mpi) {
  std::ostringstream ss;
  ss  << mpi.lb_[0]-1 << " " << mpi.ub_[0]  << " " << mpi.lb_[1]-1  << " " << mpi.ub_[1]  << " " << std::max(mpi.lb_[2]-1,0) << " " << mpi.ub_[2];
  extent = ss.str();
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

  //std::cout << n[0] << ", " << n[1] << ", " << n[2] << std::endl;

  OUTPUT_DTYPE data;
  //std::vector<int> &n = mpi_.n_;
  for (const auto& var : variables_) {
    // precede data with total number of bytes
    file_.write((char*)&(var.nbytes), sizeof(unsigned int));
    //std::vector<int> &n = var.n;
    // node loop
    int z = n.size() - 2; // 2D:z=0, 3D:z=1 to skip periodic/ghost rim
    do {
      for(int y=1; y<n[1]-1; ++y) {
        for(int x=1; x<n[0]-1; ++x) {
          int nn = x + y*n[0] + z*n[0]*n[1];
          bool write_node = true; //var.write_node.test(node->mode[nn]);
          // stride loop
          for(int dim=0; dim<var.dim; ++dim) {
            if ( (dim<n.size()) && write_node ) {
              //int element = var[i].data_step*nn + dim;
              data = var.get_data(nn, dim);
              //std::cout << data;
            } else {
              data = 0.0;
            }
            file_.write((char*)&data, sizeof(OUTPUT_DTYPE));
          }

        }
      }
      ++z;
    } while (z<n[2]-1);
  }
  // end tag
  file_ << "  </AppendedData>" << std::endl;
  //std::cout << std::endl;
}


//--------------------------------------------
//
//--------------------------------------------
void VTI_file::write_header() {
  file_ << "<?xml version=\"1.0\"?>"                                                           << std::endl;
  file_ << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"          << std::endl;
  file_ << "  <ImageData WholeExtent=\"" << extent << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
  file_ << "    <Piece Extent=\"" << extent << "\">"                                           << std::endl;
  file_ << "      <CellData " << cell_data_string_                                             << std::endl;
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
void VTI_file::set_filename_and_open() {
  filename_ = get_filename();
  file_.open(path_+filename_, std::ios::out | std::ios::binary);
}

//--------------------------------------------
//
//--------------------------------------------
std::string VTI_file::get_filename() {
  std::ostringstream ss;
  ss << std::setfill('0') << std::setw(4) << rank_ << "_" << name_ << "_" << std::setw(7) << nwrite_ << extension_;
  //set_filename_with_path(ss.str());
  return ss.str();
}

//--------------------------------------------
//
//--------------------------------------------
const std::string VTI_file::get_piece_tag() const {
  std::ostringstream oss;
  oss << "    <Piece Extent=\"" << extent << "\" Source=\"" << folders_.back() + filename_ << "\" />";
  return oss.str();
}


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
void PVTI_file::set_extent(const MPI &mpi) {
  std::ostringstream ss;
  ss << "0 " << mpi.N_[0]-2 << " " << "0 " << mpi.N_[1]-2 << " " << "0 " << std::max(mpi.N_[2]-2,0);
  extent = ss.str();
}


//--------------------------------------------
//
//--------------------------------------------
void PVTI_file::MPI_write_piece(const VTI_file &vti) {
  //std::ostringstream piece_tag;
  //piece_tag << "    <Piece Extent=\"" << vti.extent << "\" Source=\"" << vti.folders_.back() + vti.filename_ << "\" />";
  //piece_tag << "    <Piece Extent=\"" << vti.extent << "\" Source=\"" << vti.get_filename_with_path() << "\" />";
  //std::string piece_tag = vti.get_piece_tag();

  char piece[101] = {0}, piece_format[20] = {0};
  sprintf(piece_format, "%%-%ds\n", (int)(sizeof(piece))-2);
  sprintf(piece, piece_format, vti.get_piece_tag().c_str());

  MPI_File mpi_file;
  MPI_Status status;
  MPI_Bcast(&file_position, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Offset seek_position = (long long) (file_position + rank_*(sizeof(piece)-1));
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
void PVTI_file::write_header(int time, const VTI_file &vti) {
  file_.precision(5);
  file_ << "<?xml version=\"1.0\"?>"                                                   << std::endl;
  file_ << "<!-- Created " << get_timestring() << " -->"                               << std::endl;
  file_ << "<!-- time = " << time << " s -->"                                          << std::endl;
  file_ << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  file_ << "  <PImageData WholeExtent=\"" << extent << "\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
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
  oss << "Scalars=\"" << scalars.str() << "\" Vectors=\"" << vectors.str() << "\">";
  cell_data_string_ = oss.str();
}

////--------------------------------------------
////
////--------------------------------------------
//void Outfile::write(double time) {
//  // each process writes its own .vti file
//  vti_file.create(nwrite);
//  vti_file.write_header(variables, cell_data);
//  vti_file.write_data(variables, node);
//  vti_file.write_footer_and_close();
//  // root process opens and writes one .pvti-file, and
//  // each process appends its own vti-file information
//  pvti_file.set_filename(nwrite);
//  if (pvti_file.rank_ == 0) {
//    pvti_file.create();
//    pvti_file.write_header(time, variables, cell_data);
//    pvti_file.set_position_and_close();
//  }
//  pvti_file.write_piece(vti_file);
//  if (pvti_file.rank_ == pvti_file.max_rank_) {
//    pvti_file.write_footer_and_close();
//  }
//}

//--------------------------------------------
//
//--------------------------------------------
void VTI_file::write() {
  // each process writes its own .vti file
  set_filename_and_open();
  write_header();
  write_data();
  write_footer_and_close();
  ++nwrite_;
}

//--------------------------------------------
//
//--------------------------------------------
void PVTI_file::write(const double time, const VTI_file &vti) {
  // root process opens and writes one .pvti-file, and
  // each process appends its own vti-file information
  set_filename();
  if (rank_ == 0) {
    open();
    write_header(time, vti);
    set_position_and_close();
  }
  MPI_write_piece(vti);
  if (rank_ == max_rank_) {
    write_footer_and_close();
  }
  ++nwrite_;
}


//--------------------------------------------
//
//--------------------------------------------
void Outfile::add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
    const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides)
{
  if ( (names.size()!=dims.size()) || (dims.size()!=data_ptrs.size()) ) {
    std::cerr << "ERROR in Outfile::add_variables: sizes of the input arrays must equal!" << std::endl;
    exit(-1);
  }

  vti_file_.add_variables(names, data_ptrs, datasizes, dims, data_strides);
}

//--------------------------------------------
//
//--------------------------------------------
void VTI_file::add_variables(const std::vector<std::string> &names, const std::vector<void*> &data_ptrs,
    const std::vector<size_t> &datasizes, const std::vector<int> &dims, const std::vector<int> &data_strides)
{
  for (std::size_t i=0; i<dims.size(); ++i) {
    variables_.emplace_back(names[i], data_ptrs[i], datasizes[i], dims[i], data_strides[i], n);
    std::size_t end = variables_.size()-1;
    if (end>0)
      // set the offset in bytes between variables in the same .vti-file
      variables_[end].set_offset(variables_[end-1]);
  }
  set_cell_data_string();
}

//--------------------------------------------
//
//--------------------------------------------
Outfile& Output::add_file(const std::string &_name) {
  file.emplace_back(path, _name, mpi_, node);
  get_index[_name] = file.size()-1;
  return file.back();
};


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


