#ifndef _FTK_PRINT_H
#define _FTK_PRINT_H

#include <iostream>

namespace ftk {

template <typename T>
inline void print_polynomial(const std::string &name, const T P[], int n)
{
  std::cerr << name << "(x) = ";
  for (int i = 0; i <= n; i ++)
    std::cerr << P[i] << "*x^" << i << " ";
  std::cerr << std::endl;
}

template <typename T, int m, int n>
inline void print_matrix(std::ostream &os, const T A[m][n])
{
  for (int i = 0; i < m; i ++) {
    if (i == 0) os << "[";
    for (int j = 0; j < n; j ++) {
      if (j == 0) os << "[";
      os << A[i][j];
      if (j == n-1) os << "]";
      else os << ", ";
    }
    if (i == m-1) os << "]";
    else os << ", ";
  }
}

template <typename T, int m, int n>
inline void print_matrix(const std::string& name, const T A[m][n])
{
  std::cerr << name << "=";
  print_matrix<T, m, n>(std::cerr, A);
  std::cerr << std::endl;
}

template <typename T>
inline void print2x2(const std::string &name, const T m[2][2])
{
  print_matrix<T, 2, 2>(name, m);
}

template <typename T>
inline void print2x2(const std::string &name, const T m[])
{
  fprintf(stderr, "%s=[[%.10f, %.10f], [%.10f, %.10f]]\n",
      name.c_str(), m[0], m[1], m[2], m[3]);
}

template <typename T>
inline void print3x3(const std::string &name, const T m[])
{
  fprintf(stderr, "%s=[[%.10f, %.10f, %.10f], [%.10f, %.10f, %.10f], [%.10f, %.10f, %.10f]]\n",
      name.c_str(), m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
}

template <typename T>
inline void print3(const std::string &name, const T m[3])
{
  fprintf(stderr, "%s=[%.10f, %.10f, %.10f]\n", 
      name.c_str(), m[0], m[1], m[2]);
}

template <typename T>
inline void print2x3(const std::string &name, const T m[2][3])
{
  print_matrix<T, 2, 3>(name, m);
}

template <typename T>
inline void print3x2(const std::string &name, const T m[3][2])
{
  print_matrix<T, 3, 2>(name, m);
}

template <typename T>
inline void print3x3(const std::string &name, const T m[3][3])
{
  print_matrix<T, 3, 3>(name, m);
}

template <typename T>
inline void print4x3(const std::string &name, const T A[4][3])
{
  print_matrix<T, 4, 3>(name, A);
}

template <typename T>
inline void print4x4(const std::string &name, const T m[])
{
  fprintf(stderr, "%s=[[%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f]]\n",
      name.c_str(), m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);
}

template <typename T>
inline void print4x4(const std::string &name, const T m[4][4])
{
  print_matrix<T, 4, 4>(name, m);
}

}

#endif
