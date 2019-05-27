/*
 * vector_func.h
 *
 *  Created on: Dec 4, 2018
 *      Author: janlv
 */

#ifndef VECTOR_FUNC_H_
#define VECTOR_FUNC_H_
#include <vector>
#include <sstream>
//#include <iomanip>

//------------------------------------
// vec1 + vec2 operator for std::vector
//------------------------------------
template <typename T, typename U>
inline std::vector<T> operator+(const std::vector<T> &vec1, const std::vector<U> &vec2){
  std::vector<T> ret(vec1.size());
  for (size_t i=0; i<ret.size(); ++i) {
    ret[i] = vec1[i] + vec2[i];
  }
  return ret;
};

//------------------------------------
// vec1 * vec2 operator for std::vector
//------------------------------------
template <typename T, typename U>
inline std::vector<T> operator*(const std::vector<T> &vec1, const std::vector<U> &vec2){
  std::vector<T> ret(vec1.size());
  for (size_t i=0; i<ret.size(); ++i) {
    ret[i] = vec1[i] * vec2[i];
  }
  return ret;
};

//------------------------------------
// vec1 % vec2 operator for std::vector
//------------------------------------
template <typename T, typename U>
inline std::vector<T> operator%(const std::vector<T> &vec1, const std::vector<U> &vec2){
  std::vector<T> ret(vec1.size());
  for (size_t i=0; i<ret.size(); ++i) {
    ret[i] = vec1[i] % vec2[i];
  }
  return ret;
};

//------------------------------------
// vec + a operator for std::vector
//------------------------------------
template <typename T, typename U>
inline std::vector<T> operator+(const std::vector<T> &vec, const U a){
  std::vector<T> ret = vec;
  for (auto& v:ret)
    v += a;
  return ret;
};

//------------------------------------
// vec - a operator for std::vector
//------------------------------------
template <typename T, typename U>
inline std::vector<T> operator-(const std::vector<T> &vec, const U a){
  std::vector<T> ret = vec;
  for (auto& v:ret)
    v -= a;
  return ret;
};

//------------------------------------
// vec = a operator for std::vector
//------------------------------------
template <typename T>
inline std::vector<T> set_const(std::vector<T> &vec, T a){
  std::vector<T> ret = vec;
  for (auto& v:ret)
    v = a;
  return ret;
};

//------------------------------------
// product of vector elementsvec1 * vec2 operator for std::vector
//------------------------------------
template <typename T>
inline T prod(const std::vector<T> &vec){
  T ret = 1;
  for (const auto& v : vec) {
    ret *= v;
  }
  return ret;
};

//------------------------------------
// << operator for std::vector
//------------------------------------
template <typename T>
inline std::ostream& operator<<(std::ostream &os, const std::vector<T> &v) {
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

//------------------------------------
// << operator for std::vector
//------------------------------------
//template <typename T>
inline std::ostream& operator<<(std::ostream &os, const std::vector<std::vector<int>::iterator> &itr) {
  //os << "(";
  std::string sep = "";
  //os << std::setfill(' ') << std::setw(2);
  for (const auto& i : itr) {
    os << sep << *i;
    sep = ",";
  }
  //os << ")";
  return os;
}


//------------------------------------
// make a string list of a vector separated by space (default)
//------------------------------------
template <typename T>
inline const std::string vector_as_string_list(const std::vector<T> &vec, const std::string sep=" ") {
  std::stringstream ss;
  for (const auto& v : vec) {
    ss << sep << v;
  }
  //os << ")";
  return ss.str();
}

//------------------------------------
// repeat an int value n times and return as a new string
//------------------------------------
template <typename T>
inline const std::string repeated_list(T value, int n) {
  std::stringstream ss;
  for (int i=0; i<n; ++i) {
    ss << value << " ";
  }
  return ss.str();
}

//------------------------------------
// return a string of an incremental list repeat an int value n times and return as a new string
//------------------------------------
template <typename T>
inline const std::string incremental_list(T start, T stop, T increment) {
  std::stringstream ss;
  T value = start;
  while (value <= stop) {
    ss << value << " ";
    value += increment;
  }
  return ss.str();
}









#endif /* VECTOR_FUNC_H_ */
