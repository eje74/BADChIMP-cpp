/*
 * Geometry.h
 *
 *  Created on: Feb 7, 2019
 *      Author: janlv
 */

#ifndef SRC_GEO_H_
#define SRC_GEO_H_
#include <vector>
#include <cassert>
#include "Input.h"
#include "vector_func.h"
//#include "Mpi_class.h"

class Geo {
public:
  struct Limits {
    std::vector<int> n_, lb_, ub_;
    //Limits(int size) : N_(size), n_(size), lb_(size), ub_(size){};
    void assign(std::vector<int>& dim) {
      //size_t size = dim.size();
      n_.assign(dim.begin(), dim.end());
      lb_.assign(dim.size(), 0);
      ub_.assign(dim.begin(), dim.end());
    }
    //inline const std::vector<int>& get_lower_bounds() const {return lb_;}
    //inline const std::vector<int>& get_upper_bounds() const {return ub_;}
    //inline const std::vector<int>& get_lower_bounds_without_ghosts() const {return(lb_-1*num_ghost_);}
    //inline const std::vector<int>& get_upper_bounds_without_ghosts() const {return(ub_-1*num_ghost_);}
    //inline const int get_num_elements() const {return prod(n_);}
    //inline const std::vector<int>& get_size() const {return n_;};
    //inline const size_t get_dim() const {return n_.size();}
    void print_limits() {
      std::cout << "n = " << n_ << ", lb = " << lb_ << ", ub = " << ub_;
    }
  };
  Limits local_;
  Limits global_;
  int*** labels_ = nullptr;

private:
  int num_ghost_ = 0;
  int rank_ = 0;
  std::vector<int> nodes_;
  std::vector<int> node_labels_;
  double dx_ = 0.0;
  int num_fluid_ = 0;
  int num_solid_ = 0;
  int num_bulk_ = 0;
  int num_nobulk_ = 0;
  //size_t dim_ = 0;

public:
  //Geo(const std::string& filename, Mpi& mpi) : rank_(mpi.get_rank())
  Geo(const std::string& filename, int rank) : rank_(rank)
  {
    Input geo(filename); //geo.print();
    std::vector<int> dim = geo["dim"];
    local_.assign(dim);
    global_.assign(dim);
    dx_ = geo["res"];
    nodes_ = geo["geo"];
    //std::cout << nodes_ << std::endl;
    if (prod(dim)!=int(nodes_.size())) {
      //mpi.print_error_abort("Mismatch between dimension and size of geo-file");
      std::cerr << "   !!!!ERROR!!!!   " << std::endl << "   Mismatch between dimension and size of geo-file" << std::endl;
    }
    node_labels_.resize(nodes_.size());
    //set_limits(mpi);
    //add_ghost_nodes(1);
    //set_nodes<LT>();
  }
  inline int label(const int pos) const {return node_labels_[pos];}
  inline const std::vector<int>& get_labels() const {return node_labels_;}
  inline int size() const {return int(nodes_.size());};
  inline int get_num_bulk_nodes() const {return num_bulk_;};
  inline int get_num_non_bulk_nodes() const {return num_nobulk_;};
  inline int get_num_nodes() const {return nodes_.size();};
  void print_limits() {
    std::cout << "("<< rank_ << ") ghosts: " << num_ghost_ << ", global: ";
    global_.print_limits();
    std::cout << ", local: ";
    local_.print_limits();
    std::cout << std::endl;
  };
  void print_nodes() {
    std::cout << "NODES: "<< nodes_ << std::endl;
  }
  void print_labels() {
    std::cout << "LABELS: " << node_labels_ << std::endl;
  }
  //void set_labels(int*** labels) {labels_ = labels;};
  inline int get_pos(const std::vector<int>& ind, const std::vector<int>& stride) const {
    int ret = ind[0] + ind[1]*stride[0];
    if (ind.size()>2) {
      ret += ind[2]*stride[0]*stride[1];
    }
    return(ret);
    //    if (ind.size()>2)
    //      return ind[0] + ind[1]*stride[0] + ind[2]*stride[0]*stride[1];
    //    else
    //      return ind[0] + ind[1]*stride[0];
  }
  inline const std::vector<int> get_index(const int pos) const {
    std::vector<int> ind(get_dim());
    ind[0] = pos%local_.n_[0];
    ind[1] = ((pos-ind[0])/local_.n_[0]) % local_.n_[1];
    if (ind.size()>2) {
      ind[2] = int(pos/(local_.n_[1]*local_.n_[0]));
    }
    return(ind);
  }

  //inline const int global_to_local_pos(std::vector<int>& ind) { return(get_pos(ind-lb_+1, n_)); }
  inline int local_to_global_pos(const std::vector<int>& ind) { return(get_pos(ind+local_.lb_-1, global_.n_)); }
  inline int get_local_pos(const std::vector<int>& ind) { return(get_pos(ind, local_.n_)); }
  inline int get_geofile_pos(const std::vector<int>& ind) { return(get_pos(ind+local_.lb_-2, global_.n_-2)); }
  //inline const std::vector<int>& get_global_size_without_ghosts() const { return(global_.n_-2*num_ghost_); };
  inline int get_num_ghosts(){ return(num_ghost_); };
  const std::vector<int>& get_N() const { return(global_.n_);}
  const std::vector<int>& get_n() const { return(local_.n_);}
  int get_N(int axis) const { return(global_.n_[axis]);}
  int get_n(int axis) const { return(local_.n_[axis]);}
  int get_dim() const { return(global_.n_.size());}
  inline const std::vector<int>& get_lower_bounds() const {return local_.lb_;}
  inline const std::vector<int>& get_upper_bounds() const {return local_.ub_;}
  //inline char get_data(int x, int y, int z) const {return(data_[get_pos({x,y,z},local_.n_)]);};
  void export_geo_to_3D(int*** geo_3D) {
    //for (const auto& node:nodes_) {
    for (size_t pos=0; pos<nodes_.size(); ++pos) {
      std::vector<int> i = get_index(pos);
      geo_3D[i[2]][i[1]][i[0]] = nodes_[pos];
    }
  }
  void export_labels_to_3D(int*** labels_3D) {
    //for (const auto& node:nodes_) {
    for (size_t pos=0; pos<nodes_.size(); ++pos) {
      std::vector<int> i = get_index(pos);
      labels_3D[i[2]][i[1]][i[0]] = node_labels_[pos];
    }
  }


  //------------------------------------
  //   Defines nodes as fluid bulk, fluid boundary, solid bulk, solid boundary for a 3d geometry.
  //   Here it is assumend that the geometry is periodic.
  //
  //     fluid bulk     : a fluid node with only fluid nodes in the lattice neighborhood
  //     fluid boundary : a fluid node with at least one solid node in the lattice neighborhood
  //     fluid unknown  : not yet checked fluid node
  //     solid bulk     : a solid node with only solid nodes in the lattice neighborhood
  //     solid boundary : a solid node with at least one fluid node in the lattice neighborhood
  //     solid unknown  : not yet checked solid node
  //
  //     integer labels :
  //         FLUID < 3
  //         bulk fluid     : 0
  //         boundary fluid : 1
  //         fluid unknown  : 2
  //
  //         SOLID > 2
  //         boundary solid : 3
  //         bulk solid     : 4
  //         solid unknown  : 5
  //
  //
  //    usage : analyseGeometry<DXQY>(nX, nY, nZ, geo)
  //
  //    The geo-matrix should now only contain {0, 1} for fluid or {3, 4} for solid
  //
  //
  //------------------------------------
  template<typename DXQY>
  void set_node_values_v2() {
    // In the beginning all values are unkonwn
    // fluid nodes are then set to 2 and solid nodes are set to 5
    const std::vector<int>& n = local_.n_;
    int trans[] = {2, 5};  // trans[0(fluid)] = 2 and trans[1(solid)] = 5
    for (auto& node:nodes_) {
      node = trans[node];
    }
    int pos = 0;
    for (auto& node:nodes_) {
      /* Calculate the number of fluid neighbors */
      std::vector<int> xx = get_index(pos++);
      int nFluidNeig = 0;
      for (int q=0; q<DXQY::nQNonZero_; ++q) {
        std::vector<int> nn = (xx + DXQY::c(q) + n) % n;
        nFluidNeig += nodes_[get_pos(nn,n)] < 3;
        //nFluidNeig += nodes_[get_pos((xx + DXQY::c(q) + n) % n,n)] < 3;
      }
      if (node < 3) { // node is fluid
        ++num_fluid_;
        if (nFluidNeig < DXQY::nQNonZero_)
          node = 1; // (boundary fluid) neighborhood contains solid nodes
        else
          node = 0;  // (bulk fluid)
      } else {  // node is solid
        ++num_solid_;
        if (nFluidNeig > 0)
          node = 3; // (boundary solid) neighborhood contains fluid nodes
        else
          node = 4; // (bulk solid)
      }
    }
  }


  void set_labels()
/* Tags the bulk nodes in the nodeLabel matrix. Here bulk means
 *  fluid nodes means nodes that uses the standard LB update
 *  algorithm
 *
 * Input arguments:
 *  nX, nY: system size
 *  geo: geometry matrix analyzed by function 'analyseGeometry'
 *  nodeLabel: matrix crated by function newNodeLabel
 *
 * Description:
 *  - The bulk labels are number from 1 to the total number of
 *    bulk nodes.
 *  - 0 is used as a default node label, assumed to be a dummy
 *    variable
 *  - returns an integer that is the highest bulk tag.
 *  - Assumes that geo is initialized with zeros so that the
 *    bulk nodes are taged with unique labels ranging from 1 to the
 *    total number of bulk nodes, ie. the highest bulk tag.
 */
  {
    int fluid_lbl = 0;
    int non_bulk_lbl = num_fluid_;
    for (size_t pos=0; pos<nodes_.size(); ++pos) {
      assert(node_labels_[pos]==0);  // Label should not have been previously set
      if (nodes_[pos] < 3) {         // Here we have set all fluid nodes to bulk nodes
        node_labels_[pos] = ++fluid_lbl;
        ++num_bulk_;
      }
      if (nodes_[pos] == 3) {
        node_labels_[pos] = ++non_bulk_lbl;
        ++num_nobulk_;
      }
    }
  }


private:
  //void set_limits(Mpi& mpi);
  void set_nodes(const std::vector<char>& geo);
  void add_ghost_nodes(int num);
};





//------------------------------------
//
//------------------------------------
//void Mpi::set_geometry(std::vector<char>& geo_out, const std::vector<char>& geo_in) {
//  std::vector<int> geo_n = N_ - 2;
//  std::vector<int> i(get_dim());
//  i[2] = n_.size() - 2; // 2D:z=0, 3D:z=1 to skip periodic/ghost rim
//  do {
//    for(i[1]=1; i[1]<n_[1]-1; ++i[1]) {
//      for(i[0]=1; i[0]<n_[0]-1; ++i[0]) {
//        int a = get_pos(i, n_);
//        int b = get_pos(i+lb_-2, geo_n);
//        //std::cout << a << "," << b << ":" << geo_in[b] << std::endl;
//        geo_out[a] = geo_in[b];
//      }
//    }
//    ++i[2];
//  } while (i[2]<n_[2]-1);
//}




#endif /* SRC_GEO_H_ */
