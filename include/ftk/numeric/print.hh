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

template <typename T>
inline void print2x2(const std::string &name, const T m[2][2])
{
  fprintf(stderr, "%s=[[%.10f, %.10f], [%.10f, %.10f]]\n",
      name.c_str(), m[0][0], m[0][1], m[1][0], m[1][1]);
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
inline void print3x3(const std::string &name, const T m[3][3])
{
  fprintf(stderr, "%s=[[%.10f, %.10f, %.10f], [%.10f, %.10f, %.10f], [%.10f, %.10f, %.10f]]\n",
      name.c_str(), m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]);
}

template <typename T>
inline void print3x4(const std::string &name, const T m[])
{
  fprintf(stderr, "%s=[[%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f]\n",
      name.c_str(), m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
}

template <typename T>
inline void print4x4(const std::string &name, const T m[])
{
  fprintf(stderr, "%s=[[%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f]]\n",
      name.c_str(), m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);
}

}

#endif
