#ifndef _FTK_CRITICAL_POINT_T_HH
#define _FTK_CRITICAL_POINT_T_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/critical_point_lite.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/json.hh>

namespace ftk {

using nlohmann::json;

// template <int N/*dimensionality*/, typename ValueType=double, typename IntegerType=unsigned long long>
struct critical_point_t {
  critical_point_t() {}
  critical_point_t(const critical_point_lite_t& cp) {
    for (int i = 0; i < 3; i ++)
      x[i] = cp.x[i];
    t = cp.t;
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++)
      scalar[i] = cp.scalar[i];
    type = cp.type;
    tag = cp.tag;
  }

  double operator[](size_t i) const {return x[i];}
  std::array<double, 3> x; // double x[3] = {0}; // coordinates 
  double t = 0.0; // time
  int timestep = 0; 
  // double rx[4] = {0}; // coordinates in transformed (e.g. curvilinear) grid, if eligible
  std::array<double, FTK_CP_MAX_NUM_VARS> scalar; // double scalar[FTK_CP_MAX_NUM_VARS] = {0};
  unsigned int type = 0;
  bool ordinal = false;
  unsigned long long tag = 0, id = 0;

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
    os << "tag=" << tag << ", "; 
    os << "id=" << id;  // << std::endl;
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

  void relabel(unsigned long long i) {
    identifier = i;
    for (i = 0; i < size(); i ++)
      at(i).id = i;
  }

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

      if (at(i).type != current_type || i == size() - 1) {
        if (subtraj.size() > 0) {
          subtraj.update_statistics();
          results.push_back(subtraj);
          subtraj.clear();
        }
      } 
    }

    return results;
  }

  std::vector<int/*idx in orignal traj*/> to_ordinals() const {
    std::vector<int> result;
    for (auto i = 0; i < size(); i ++)
      if (at(i).ordinal)
        result.push_back(i);
    return result;
  }

  void rotate() { // if the traj is a loop and the type is inconsistent, rotate the traj before splitting
    if (!loop) return;
    if (front().type != back().type) return;
    int i;
    for (i = 0; i < size(); i ++)
      if (front().type != at(i).type) break;
    if (i < size()) {
      // fprintf(stderr, "rotating, i=%d, size=%zu\n", i, size());
      std::rotate(begin(), begin()+i, end());
    }
  }

  void reorder() { // reorder by timestep
    if (loop) return;
    bool reverse = false;
    if (front().timestep == back().timestep) {
      if (front().t > back().t) reverse = true;
    } else {
      if (front().timestep > back().timestep)
        reverse = true;
    }
    if (reverse)
      std::reverse(std::begin(*this), std::end(*this));
  }

  void adjust_time() { // assuming traj has only one single branch and is organized in ascending order
    // ascending loop
    for (auto i = 0; i < size(); i ++) {
      if (i == 0 || at(i).ordinal) continue;
      at(i).t = std::max(at(i-1).t, at(i).t);
    }
    // descending loop
    for (auto i = size(); i -- > 0; ) {
      if (i == size() - 1 || at(i).ordinal) continue;
      at(i).t = std::min(at(i+1).t, at(i).t);
    }
  }

  void smooth_ordinal_types(const int half_window_size=2) {
    auto ordinals = to_ordinals();
    if (ordinals.size() < half_window_size*2+1) return;

    // 1. smooth ordinal types
    std::map<int/*idx in orginal traj*/, unsigned int/*target type*/> pending_changes;
    
    for (int i = half_window_size; i < ordinals.size() - half_window_size; i ++) {
      unsigned int local_consistent_type = at(ordinals[i - half_window_size]).type;
      for (int j = i - half_window_size; j <= i + half_window_size; j ++)  {
        if (j == i) continue;
        else if (local_consistent_type != at(ordinals[j]).type) {
          local_consistent_type = 0; // type is inconsistent within the window
          break;
        }
      }
        
      if (local_consistent_type != 0 && at(ordinals[i]).type != local_consistent_type)
        pending_changes[ordinals[i]] = local_consistent_type;
    }

    for (const auto &kv : pending_changes)
      at(kv.first).type = kv.second;
  }

  void smooth_interval_types() {
    const auto ordinals = to_ordinals();
    if (ordinals.empty()) return;

    // front terminal
    const unsigned int front_terminal_type = at(ordinals[0]).type;
    for (int i = 0; i < ordinals[0]; i ++) 
      at(i).type = front_terminal_type;

    // back terminal
    const unsigned int back_terminal_type = at(ordinals[ordinals.size()-1]).type;
    for (int i = ordinals[ordinals.size()-1]; i < size(); i ++)
      at(i).type = back_terminal_type;

    // intervals; preventing type oscilation within intervals
    for (int i = 0; i < ordinals.size()-1; i ++) {
      if (at(ordinals[i]).type == at(ordinals[i+1]).type) {
        const unsigned int interval_type = at(ordinals[i]).type;
        for (int j  = ordinals[i]; j < ordinals[i+1]; j ++)
          at(j).type = interval_type;
      } else { // find the first point that changes type, and then smooth the rest in the interval
        const unsigned int ltype = at(ordinals[i]).type, rtype = at(ordinals[i+1]).type;
        int j;
        for (j = ordinals[i] ; j < ordinals[i+1]; j ++)
          if (at(j).type != ltype) break;
        for (; j < ordinals[i+1]; j ++)
          at(j).type = rtype;
      }
    }
  }
  
  int locate(double t, bool cap=false/*lower (false) or higher (true) cap*/) const { // locate interval id of given t; assuming the traj is already reordered
    return 0; // TODO; WIP
  }

  critical_point_traj_t intercept(int t0, int t1) // assuming the traj is already reordered
  {
    critical_point_traj_t result;
    if (t0 > back().t || t1 < front().t) return result;
    return result;

    int i0, i1;
    // for (i0 = 0; i0 < size(); i ++)
  }
};

}

// serialization w/ json
namespace nlohmann
{
  using namespace ftk;

  template <>
  struct adl_serializer<critical_point_t> {
    static void to_json(json &j, const critical_point_t& cp) {
      j["x"] = cp.x; // {cp.x[0], cp.x[1], cp.x[2]};
      j["t"] = cp.t;
      j["timestep"] = cp.timestep;
      j["scalar"] = cp.scalar; // std::vector<double>(cp.scalar, cp.scalar+FTK_CP_MAX_NUM_VARS);
      j["type"] = cp.type;
      j["ordinal"] = cp.ordinal;
      j["tag"] = cp.tag;
      j["id"] = cp.id;
    }

    static void from_json(const json& j,critical_point_t& cp) {
      cp.x = j["x"];  // for (int i = 0; i < 3; i ++) cp.x[i] = j["x"][i];
      cp.t = j["t"];
      cp.timestep = j["timestep"];
      cp.scalar = j["scalar"];  // for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++) cp.scalar[i] = j["scalar"][i];
      cp.type = j["type"];
      cp.ordinal = j["ordinal"];
      cp.tag = j["tag"];
      cp.id = j["id"];
    }
  };
  
  template <>
  struct adl_serializer<critical_point_traj_t> {
    static void to_json(json &j, const critical_point_traj_t& t) {
      j = {
        {"id", t.identifier},
        {"max", t.max},
        {"min", t.min},
        {"persistence", t.persistence},
        {"bbmin", t.bbmin},
        {"bbmax", t.bbmax},
        {"tmin", t.tmin},
        {"tmax", t.tmax},
        {"consistent_type", t.consistent_type},
        {"traj", static_cast<std::vector<critical_point_t>>(t)}
      };
    }

    static void from_json(const json&j, critical_point_traj_t& t) {
      t.identifier = j["id"];
      t.max = j["max"];
      t.min = j["min"];
      t.persistence = j["persistence"];
      t.bbmin = j["bbmin"];
      t.bbmax = j["bbmax"];
      t.tmin = j["tmin"];
      t.tmax = j["tmax"];
      t.consistent_type = j["consistent_type"];
      std::vector<critical_point_t> traj = j["traj"];
      t.clear();
      t.insert(t.begin(), traj.begin(), traj.end());
    }
  };
}



// serialization
namespace diy {
  // template <int N, typename V, typename I> 
  // static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t<N, V, I> &cp) {
  static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t &cp) {
    diy::save(bb, cp.x); // for (int i = 0; i < 3; i ++) diy::save(bb, cp.x[i]);
    diy::save(bb, cp.t);
    diy::save(bb, cp.scalar); // for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++) diy::save(bb, cp.scalar[i]);
    diy::save(bb, cp.type);
    diy::save(bb, cp.ordinal);
    diy::save(bb, cp.tag);
    diy::save(bb, cp.id);
  }

  // template <int N, typename V, typename I> 
  // static void load(diy::BinaryBuffer& bb, ftk::critical_point_t<N, V, I> &cp) {
  static void load(diy::BinaryBuffer& bb, ftk::critical_point_t &cp) {
    diy::load(bb, cp.x); // for (int i = 0; i < 4; i ++) diy::load(bb, cp.x[i]);
    diy::load(bb, cp.t);  
    diy::load(bb, cp.scalar); // for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++) diy::load(bb, cp.scalar[i]);
    diy::load(bb, cp.type);
    diy::load(bb, cp.ordinal);
    diy::load(bb, cp.tag);
    diy::load(bb, cp.id);
  }
  
  static void save(diy::BinaryBuffer& bb, const ftk::critical_point_traj_t &t) {
    diy::save(bb, t.complete);
    diy::save(bb, t.max);
    diy::save(bb, t.min);
    diy::save(bb, t.persistence);
    diy::save(bb, t.bbmin);
    diy::save(bb, t.bbmax);
    diy::save(bb, t.tmin);
    diy::save(bb, t.tmax);
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
    diy::load(bb, t.tmin);
    diy::load(bb, t.tmax);
    diy::load(bb, t.consistent_type);
    diy::load<std::vector<ftk::critical_point_t>>(bb, t);
  }
}

#endif
