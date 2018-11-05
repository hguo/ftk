#ifndef _FTK_PRINT_H
#define _FTK_PRINT_H

#include <iostream>

namespace ftk {

template <typename ValueType>
inline void print2x2(const std::string &name, const ValueType m[])
{
  fprintf(stderr, "%s=[[%.10f, %.10f], [%.10f, %.10f]\n",
      name.c_str(), m[0], m[1], m[2], m[3]);
}

template <typename ValueType>
inline void print3x3(const std::string &name, const ValueType m[])
{
  fprintf(stderr, "%s=[[%.10f, %.10f, %.10f], [%.10f, %.10f, %.10f], [%.10f, %.10f, %.10f]]\n",
      name.c_str(), m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
}

template <typename ValueType>
inline void print4x4(const std::string &name, const ValueType m[])
{
  fprintf(stderr, "%s=[[%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f], [%.10f, %.10f, %.10f, %.10f]]\n",
      name.c_str(), m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);
}

}

#endif
