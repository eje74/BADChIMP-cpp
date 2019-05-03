#ifndef FIELD_H
#define FIELD_H

#include "LBglobal.h"
#include <iostream>
#include <vector>

class Field {
public:
  Field(const int nFields, const int nNodes, const int nDim) :
    nFields_(nFields), nNodes_(nNodes), nDim_(nDim), stride_(nFields*nDim), data_(nFields*nNodes*nDim) { };
  int getNumNodes() const {return nNodes_;} // Getter for nNodes_
  int num_fields() const {return nFields_;}

  inline void set(const std::vector<lbBase_t> val) {
    if (int(val.size()) != nFields_*nDim_) {
      std::cerr << std::endl << "WARNING in Field::set(): Wrong number of constants given, "
          << val.size() << " != " << nFields_*nDim_ << std::endl << std::endl;
    }
    for (int f=0; f<nFields_; ++f) {
      for (size_t i=f+1; i<data_.size(); i+=stride_) {
        for (int d=0; d<nDim_; ++d) {
          data_[i] = val[f*nDim_+d];
        }
      }
    }
  };

  inline void set(const int fieldNo, const std::vector<lbBase_t> val) {
    // add 1 due to dummy node
    if (fieldNo > nFields_-1) {
      std::cerr << std::endl << "WARNING in Field::set(): Number of fields exceeded, "
          << fieldNo << " > " << nFields_-1 << std::endl << std::endl;
    }
    for (size_t i=fieldNo+1; i<data_.size(); i+=stride_) {
      for (int d=0; d<nDim_; ++d) {
        data_[i] = val[d];
      }
    }
  };

  inline void set(const std::vector<int> nodes, const std::vector<lbBase_t> val) {
    for (const auto& n:nodes) {
      int nn = stride_*n;
      for (int f=0; f<nFields_; ++f) {
        data_[nn + f] = val[f];
      }
    }
  };

  void print() { std::cout << data_ << std::endl; };

protected:
  int nFields_ = 0;
  int nNodes_ = 0;
  int nDim_ = 0;
  int stride_ = 0;
  std::vector<lbBase_t> data_;
};

class Scalar_field : public Field {
public:
  Scalar_field(const int nFields, const int nNodes) : Field(nFields, nNodes, 1) {};

  inline const lbBase_t& operator()(const int fieldNo, const int nodeNo) const {
    return data_[stride_*nodeNo + fieldNo];
  }
  inline lbBase_t& operator()(const int fieldNo, const int nodeNo) {
    return data_[stride_*nodeNo + fieldNo];
  }


};

class Vector_field : public Field {
public:
  Vector_field(const int nFields, const int nNodes, const int nDim) : Field(nFields, nNodes, nDim) {};

  inline const lbBase_t& operator()(const int fieldNo, const int nodeNo, const int dim) const {
    return data_[stride_*nodeNo + nDim_*fieldNo + dim];
  }
  inline std::vector<lbBase_t> operator()(const int fieldNo, const int nodeNo) {
    std::vector<lbBase_t>::iterator it = data_.begin() + (stride_*nodeNo+nDim_*fieldNo);
    std::vector<lbBase_t> vec(it, it+nDim_);
    return vec;
  }
};

class LB_field : public Field {
public:
  LB_field(const int nFields, const int nNodes, const int nDim) : Field(nFields, nNodes, nDim) {};

  inline const lbBase_t& operator()(const int fieldNo, const int nodeNo, const int dir) const {
    return data_[stride_*nodeNo + nDim_*fieldNo + dir];
  }
  inline std::vector<lbBase_t> operator()(const int fieldNo, const int nodeNo) {
    std::vector<lbBase_t>::iterator it = data_.begin() + (stride_*nodeNo+nDim_*fieldNo);
    return std::vector<lbBase_t>(it, it+nDim_);
  }
  inline void swapData(LB_field& field) {
    data_.swap(field.data_);
  }
};


#endif // FIELD_H
