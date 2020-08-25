#ifndef _FTK_CRITICAL_POINT_T_HH
#define _FTK_CRITICAL_POINT_T_HH

#include <ftk/external/diy/serialization.hpp>

namespace ftk {

// template <int N/*dimensionality*/, typename ValueType=double, typename IntegerType=unsigned long long>
struct critical_point_t {
  double operator[](size_t i) const {return x[i];}
  double x[3] = {0}; // coordinates 
  double t = 0.0; // time
  int timestep = 0; 
  // double rx[4] = {0}; // coordinates in transformed (e.g. curvilinear) grid, if eligible
  double scalar[FTK_CP_MAX_NUM_VARS] = {0};
  unsigned int type = 0;
  unsigned long long tag = 0;
  bool ordinal = false;

  // constexpr size_t size() const noexcept { return sizeof(critical_point_t<N, ValueType, IntegerType>); }
  constexpr size_t size() const noexcept { return sizeof(critical_point_t); }

  std::ostream& print(std::ostream& os, const int cpdims, const std::vector<std::string>& scalar_components) const {
    if (cpdims == 2) os << "x=(" << x[0] << ", " << x[1] << "), ";
    else os << "x=(" << x[0] << ", " << x[1] << ", " << x[2] << "), ";
    os << "t=" << t << ", ";

    for (int k = 0; k < scalar_components.size(); k ++)
      os << scalar_components[k] << "=" << scalar[k] << ", ";
    
    os << "type=" << critical_point_type_to_string(cpdims, type, scalar_components.size()) << ", "; 
    os << "timestep=" << timestep << ", ";
    os << "ordinal=" << ordinal << ", ";
    os << "tag=" << tag; // << std::endl;
    return os;
  }
};

struct critical_point_traj_t : public std::vector<critical_point_t>
{
  unsigned long long identifier;
  bool complete = false, loop = false;
  std::array<double, FTK_CP_MAX_NUM_VARS> max, min, persistence;
  std::array<double, 3> bbmin, bbmax; // bounding box
  double tmin, tmax; // time bounding box
  unsigned int consistent_type = 0; // 0 if no consistent type

  void discard_interval_points() {
    critical_point_traj_t traj;
    traj.identifier = identifier;
    traj.loop = loop;
    for (auto i = 0; i < size(); i ++) {
      if (at(i).ordinal) 
        traj.push_back(at(i));
    }
    traj.update_statistics();
    *this = traj;
  }
  
  void discard_degenerate_points() {
    critical_point_traj_t traj;
    traj.identifier = identifier;
    traj.loop = loop;
    for (auto i = 0; i < size(); i ++) {
      if (at(i).type != 1 && at(i).type != 0) // degenerate or unknown
        traj.push_back(at(i));
    }
    traj.update_statistics();
    *this = traj;
  }

  void update_statistics() {
    if (empty()) return; // nothing to do

    max.fill( std::numeric_limits<double>::lowest() );
    min.fill( std::numeric_limits<double>::max() );
    bbmax.fill( std::numeric_limits<double>::lowest() );
    bbmin.fill( std::numeric_limits<double>::max() );
    tmax = std::numeric_limits<double>::lowest();
    tmin = std::numeric_limits<double>::max();

    for (auto i = 0; i < size(); i ++) {
      for (int k = 0; k < FTK_CP_MAX_NUM_VARS; k ++) {
        max[k] = std::max(max[k], at(i).scalar[k]);
        min[k] = std::min(min[k], at(i).scalar[k]);
      }
      for (int k = 0; k < 3; k ++) {
        bbmax[k] = std::max(bbmax[k], at(i).x[k]);
        bbmin[k] = std::min(bbmin[k], at(i).x[k]);
      }
      tmax = std::max(tmax, at(i).t);
      tmin = std::min(tmin, at(i).t);
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
 
  std::vector<critical_point_traj_t> to_consistent_sub_traj() const {
    std::vector<critical_point_traj_t> results;
    critical_point_traj_t subtraj;
    unsigned int current_type;

    for (auto i = 0; i < size(); i ++) {
      if (subtraj.empty())
        current_type = at(i).type;

      if (at(i).type == current_type)
        subtraj.push_back(at(i));
      else {
        subtraj.update_statistics();
        results.push_back(subtraj);
        subtraj.clear();
      }
    }

    return results;
  }
};

}

// serialization
namespace diy {
  // template <int N, typename V, typename I> 
  // static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t<N, V, I> &cp) {
  static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t &cp) {
    for (int i = 0; i < 3; i ++)
      diy::save(bb, cp.x[i]);
    diy::save(bb, cp.t);
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++)
      diy::save(bb, cp.scalar[i]);
    diy::save(bb, cp.type);
    diy::save(bb, cp.tag);
    diy::save(bb, cp.ordinal);
  }

  // template <int N, typename V, typename I> 
  // static void load(diy::BinaryBuffer& bb, ftk::critical_point_t<N, V, I> &cp) {
  static void load(diy::BinaryBuffer& bb, ftk::critical_point_t &cp) {
    for (int i = 0; i < 4; i ++) 
      diy::load(bb, cp.x[i]);
    diy::save(bb, cp.t);
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++)
      diy::load(bb, cp.scalar[i]);
    diy::load(bb, cp.type);
    diy::load(bb, cp.tag);
    diy::load(bb, cp.ordinal);
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
