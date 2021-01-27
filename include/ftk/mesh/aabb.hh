#ifndef _FTK_AABB_H
#define _FTK_AABB_H

#include <ftk/config.hh>

namespace ftk {

template <typename T>
inline T min3(T x, T y, T z) {
  return std::min(std::min(x, y), z);
}

template <typename T>
inline T max3(T x, T y, T z) {
  return std::max(std::max(x, y), z);
}

template <typename T>
struct AABB {
  int id = 0;
  T A[2] = {std::numeric_limits<T>::max(), std::numeric_limits<T>::max()}, 
    B[2] = {-std::numeric_limits<T>::max(), -std::numeric_limits<T>::max()};
  T C[2] = {0, 0}; // centroid

  bool contains(const T X[]) const {
    return X[0] >= A[0] && X[0] < B[0] && X[1] >= A[1] && X[1] < B[1];
  }

  void update_centroid() {
    C[0] = (A[0] + B[0]) / 2;
    C[1] = (A[1] + B[1]) / 2;
  }

  void print() const {
    fprintf(stderr, "A={%f, %f}, B={%f, %f}, centroid={%f, %f}\n", 
        A[0], A[1], B[0], B[1], C[0], C[1]);
  }
};

}

#endif
