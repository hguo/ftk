#ifndef _FTK_CRITICAL_POINT_T_HH
#define _FTK_CRITICAL_POINT_T_HH

#include <ftk/external/diy/serialization.hpp>

namespace ftk {

// template <int N/*dimensionality*/, typename ValueType=double, typename IntegerType=unsigned long long>
struct critical_point_t {
  double operator[](size_t i) const {return x[i];}
  double x[4] = {0}; // coordinates in (spacetime) cartisian grid
  // double rx[4] = {0}; // coordinates in transformed (e.g. curvilinear) grid, if eligible
  double scalar[FTK_CP_MAX_NUM_VARS] = {0};
  unsigned int type = 0;
  unsigned long long tag = 0;
  bool ordinal = false;

  // constexpr size_t size() const noexcept { return sizeof(critical_point_t<N, ValueType, IntegerType>); }
  constexpr size_t size() const noexcept { return sizeof(critical_point_t); }
};

struct critical_point_traj_t : public std::vector<critical_point_t>
{
  bool complete = false;
  std::array<double, FTK_CP_MAX_NUM_VARS> max, min, persistence;
  std::array<double, 4> bbmin, bbmax; // bounding box
  unsigned int consistent_type = 0; // 0 if no consistent type
  
  void update_statistics() {
    max.fill( std::numeric_limits<double>::lowest() );
    min.fill( std::numeric_limits<double>::max() );
    bbmax.fill( std::numeric_limits<double>::lowest() );
    bbmin.fill( std::numeric_limits<double>::max() );

    for (auto i = 0; i < size(); i ++) {
      for (int k = 0; k < FTK_CP_MAX_NUM_VARS; k ++) {
        max[k] = std::max(max[k], at(i).scalar[k]);
        min[k] = std::min(min[k], at(i).scalar[k]);
      }
      for (int k = 0; k < 4; k ++) {
        bbmax[k] = std::max(bbmax[k], at(i).x[k]);
        bbmin[k] = std::min(bbmin[k], at(i).x[k]);
      }
    }

    for (int k = 0; k < FTK_CP_MAX_NUM_VARS; k ++)
      persistence[k] = max[k] - min[k];

    consistent_type = at(0).type;
    for (auto i = 0; i < size(); i ++)
      if (consistent_type != at(i).type) {
        consistent_type = 0;
        break;
      }
  }
};

}

// serialization
namespace diy {
  // template <int N, typename V, typename I> 
  // static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t<N, V, I> &cp) {
  static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t &cp) {
    for (int i = 0; i < 4; i ++) {
      diy::save(bb, cp.x[i]);
      // diy::save(bb, cp.rx[i]);
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
      // diy::load(bb, cp.rx[i]);
    }
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++)
      diy::load(bb, cp.scalar[i]);
    diy::load(bb, cp.type);
    diy::save(bb, cp.tag);
  }
  
  static void save(diy::BinaryBuffer& bb, const ftk::critical_point_traj_t &t) {
    diy::save(bb, t.complete);
    diy::save(bb, t.max);
    diy::save(bb, t.min);
    diy::save(bb, t.persistence);
    diy::save(bb, t.bbmin);
    diy::save(bb, t.bbmax);
    diy::save(bb, t.consistent_type);
    diy::save<std::vector<ftk::critical_point_t>>(bb, t);
  }
  
  static void load(diy::BinaryBuffer& bb, ftk::critical_point_traj_t &t) {
    diy::load(bb, t.complete);
    diy::load(bb, t.max);
    diy::load(bb, t.min);
    diy::load(bb, t.persistence);
    diy::load(bb, t.bbmin);
    diy::load(bb, t.bbmax);
    diy::load(bb, t.consistent_type);
    diy::load<std::vector<ftk::critical_point_t>>(bb, t);
  }
}

#endif
