#ifndef _FTK_CRITICAL_POINT_T_HH
#define _FTK_CRITICAL_POINT_T_HH

#include <ftk/external/diy/serialization.hpp>

namespace ftk {

template <int N/*dimensionality*/, typename ValueType=double, typename IntegerType=unsigned long long>
struct critical_point_t {
  ValueType operator[](size_t i) const {if (i >= N) return 0; else return x[i];}
  ValueType x[N];
  ValueType scalar = ValueType(0);
  unsigned int type = 0;
  IntegerType tag = 0;
};

}

// serialization
namespace diy {
  // struct Serialization<ftk::critical_point_t<N, V, I>> {
    template <int N, typename V, typename I> 
    static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t<N, V, I> &cp) {
      for (int i = 0; i < N; i ++)
        diy::save(bb, cp.x[i]);
      diy::save(bb, cp.scalar);
      diy::save(bb, cp.type);
      diy::save(bb, cp.tag);
    }

    template <int N, typename V, typename I> 
    static void load(diy::BinaryBuffer& bb, ftk::critical_point_t<N, V, I> &cp) {
      for (int i = 0; i < N; i ++)
        diy::load(bb, cp.x[i]);
      diy::load(bb, cp.scalar);
      diy::load(bb, cp.type);
    }
  // };
}

#endif
