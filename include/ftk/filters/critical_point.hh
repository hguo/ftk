#ifndef _FTK_CRITICAL_POINT_T_HH
#define _FTK_CRITICAL_POINT_T_HH

#include <ftk/external/diy/serialization.hpp>

namespace ftk {

// template <int N/*dimensionality*/, typename ValueType=double, typename IntegerType=unsigned long long>
struct critical_point_t {
  double operator[](size_t i) const {return x[i];}
  double x[4] = {0}; // coordinates in (spacetime) cartisian grid
  double rx[4] = {0}; // coordinates in transformed (e.g. curvilinear) grid, if eligible
  double scalar[FTK_CP_MAX_NUM_VARS] = {0};
  unsigned int type = 0;
  unsigned long long tag = 0;

  // constexpr size_t size() const noexcept { return sizeof(critical_point_t<N, ValueType, IntegerType>); }
  constexpr size_t size() const noexcept { return sizeof(critical_point_t); }
};

}

// serialization
namespace diy {
  // template <int N, typename V, typename I> 
  // static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t<N, V, I> &cp) {
  static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t &cp) {
    for (int i = 0; i < 4; i ++) {
      diy::save(bb, cp.x[i]);
      diy::save(bb, cp.rx[i]);
    }
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++)
      diy::save(bb, cp.scalar[i]);
    diy::save(bb, cp.type);
    diy::save(bb, cp.tag);
  }

  // template <int N, typename V, typename I> 
  // static void load(diy::BinaryBuffer& bb, ftk::critical_point_t<N, V, I> &cp) {
  static void load(diy::BinaryBuffer& bb, ftk::critical_point_t &cp) {
    for (int i = 0; i < 4; i ++) {
      diy::load(bb, cp.x[i]);
      diy::load(bb, cp.rx[i]);
    }
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++)
      diy::load(bb, cp.scalar[i]);
    diy::load(bb, cp.type);
    diy::save(bb, cp.tag);
  }
}

#endif
