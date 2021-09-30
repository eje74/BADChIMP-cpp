/*
 * vector_func.h
 *
 *  Created on: Dec 4, 2018
 *      Author: janlv
 */

#ifndef VECTOR_FUNC_H_
#define VECTOR_FUNC_H_
#include <vector>
#include <array>
#include <sstream>
#include <numeric>
#include <math.h>
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
// vec1 - vec2 operator for std::vector
//------------------------------------
template <typename T, typename U>
inline std::vector<T> operator-(const std::vector<T> &vec1, const std::vector<U> &vec2){
  std::vector<T> ret(vec1.size());
  for (size_t i=0; i<ret.size(); ++i) {
    ret[i] = vec1[i] - vec2[i];
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
  std::vector<T> ret(vec);
  for (auto& v:ret)
    v += a;
  return ret;
};

//------------------------------------
// vec - a operator for std::vector
//------------------------------------
template <typename T, typename U>
inline std::vector<T> operator-(const std::vector<T> &vec, const U a){
  std::vector<T> ret(vec);
  for (auto& v:ret)
    v -= a;
  return ret;
};

//------------------------------------
// vec / a operator for std::vector
//------------------------------------
template <typename T, typename U>
inline std::vector<T> operator/(const std::vector<T> &vec, const U a){
  std::vector<T> ret = vec;
  for (auto& v:ret)
    v /= a;
  return ret;
};

//------------------------------------
// vec * a operator for std::array
//------------------------------------
template <typename T, typename U>
inline std::array<T,3> operator*(const std::array<T,3> &vec, const U a){
  std::array<T,3> ret = vec;
  for (auto& v:ret)
    v /= a;
  return ret;
};

//------------------------------------
// norm for std::vector
//------------------------------------
template <typename T>
inline T norm(const std::vector<T> &vec){
  T ret = 0.0; 
  for (auto& v:vec)
    ret += v*v;
  return sqrt(ret);
};

//------------------------------------
// norm for std::vector
//------------------------------------
template <typename T>
inline T norm(T x, T y, T z){
  T ret = 0; 
  std::vector<T> vec {x,y,z};
  for (auto& v:vec)
    ret += v*v;
  return sqrt(ret);
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
inline T prod(const std::vector<T> &v1, const std::vector<T> &v2){
  T ret = 0;
  for (int i=0; i<v1.size(); i++) {
    ret += v1[i]*v2[i];
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
// make a string list of a vector separated by space (default)
//------------------------------------
template <typename T>
inline const std::string vector_as_string_list(const std::vector<T> &vec, int start, int stop, double min=-1, const std::string sep=" ") {
  std::stringstream ss;
  for (int ind=start; ind<stop; ++ind) {
    if (min>0 && std::abs(vec[ind])<min)
      ss << sep << 0.0;
    else
      ss << sep << vec[ind];
  } 
  return ss.str();
}

//------------------------------------
// make a string list of a vector separated by space (default)
//------------------------------------
template <typename T>
inline const std::string vector_as_string_list(const std::vector<T> &vec, const std::vector<int> &index, double min=-1, const std::string sep=" ") {
  std::stringstream ss;
  for (const auto& ind : index) {
    if (min>0 && std::abs(vec[ind])<min)
      ss << sep << 0.0;
    else
      ss << sep << vec[ind];
  } 
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

//------------------------------------
// C++ version of cumsum in numpy
//------------------------------------
template <typename T>
inline void cumsum(std::vector<T>& vec) {
  std::partial_sum(vec.begin(), vec.end(), vec.begin(), std::plus<T>());
}

//------------------------------------
// 
//------------------------------------
template <typename T>
inline std::vector<T> linspace(T start, T stop, T inc) {
  std::vector<T> vec((int)((stop-start)/inc)+1, inc);
  vec[0] = start;
  std::partial_sum(vec.begin(), vec.end(), vec.begin(), std::plus<T>());
  return vec;
}






#endif /* VECTOR_FUNC_H_ */
