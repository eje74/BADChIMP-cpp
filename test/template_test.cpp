#include <iostream>


template<int N>
inline double cu(const double *vec)
{
  return 0.0;
}

template<>
inline double cu<1>(const double *vec)
{
  return vec[0];
}

template<>
inline double cu<2>(const double *vec)
{
  return vec[1];
}

template<>
inline double cu<3>(const double *vec)
{
  return -vec[0];
}


template<>
inline double cu<4>(const double *vec)
{
  return -vec[1];
}


template<int N>
inline void forLoop(const double *vec) {

  std::cout << N << std::endl;
  std::cout << cu<N>(vec) << std::endl << std::endl; 

  // Recursive 
  forLoop<N-1>(vec);
}


template<>
void forLoop<0>(const double *vec) {
  std::cout << "0" << std::endl;
}


int main() 
{
  double vec[2] = {1.0, 2.0};

  forLoop<4>(vec);
  return 0;
}
