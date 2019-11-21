#ifndef _FTK_CRITICAL_POINT_T_HH
#define _FTK_CRITICAL_POINT_T_HH

namespace ftk {

template <int N, typename T>
struct critical_point_t {
  T operator[](size_t i) const {if (i >= N) return 0; else return x[i];}
  T x[N];
  T scalar = T(0);
  unsigned int type = 0;
  unsigned long long tag = 0;
};

}

#endif
