/*
 * vector_func.h
 *
 *  Created on: Dec 4, 2018
 *      Author: janlv
 */

#ifndef VECTOR_FUNC_H_
#define VECTOR_FUNC_H_
#include <vector>
//#include <iomanip>

//------------------------------------
// vec + vec operator for std::vector
//------------------------------------
template <typename T, typename U>
std::vector<T> operator+(const std::vector<T> &vec1, const std::vector<U> &vec2){
  std::vector<T> ret(vec1.size());
  for (size_t i=0; i<ret.size(); ++i) {
    ret[i] = vec1[i] + vec2[i];
  }
  return ret;
};

//------------------------------------
// vec + a operator for std::vector
//------------------------------------
template <typename T, typename U>
std::vector<T> operator+(const std::vector<T> &vec, const U a){
  std::vector<T> ret = vec;
  for (auto& v:ret)
    v += a;
  return ret;
};

//------------------------------------
// vec - a operator for std::vector
//------------------------------------
template <typename T, typename U>
std::vector<T> operator-(const std::vector<T> &vec, const U a){
  std::vector<T> ret = vec;
  for (auto& v:ret)
    v -= a;
  return ret;
};

//------------------------------------
// << operator for std::vector
//------------------------------------
template <typename T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &v) {
  //os << "(";
  std::string sep = "";
  //os << std::setfill(' ') << std::setw(2);
  for (const auto& i : v) {
    os << sep << i;
    sep = ",";
  }
  //os << ")";
  return os;
}





#endif /* VECTOR_FUNC_H_ */
