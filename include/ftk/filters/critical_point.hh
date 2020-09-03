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
  std::array<double, 3> x = {0}; // double x[3] = {0}; // coordinates 
  double t = 0.0; // time
  int timestep = 0; 
  // double rx[4] = {0}; // coordinates in transformed (e.g. curvilinear) grid, if eligible
  std::array<double, FTK_CP_MAX_NUM_VARS> scalar = {0}; // double scalar[FTK_CP_MAX_NUM_VARS] = {0};
  std::array<double, 3> v = {0};
  unsigned int type = 0;
  bool ordinal = false;
  unsigned long long tag = 0, id = 0;

  // constexpr size_t size() const noexcept { return sizeof(critical_point_t<N, ValueType, IntegerType>); }
  constexpr size_t size() const noexcept { return sizeof(critical_point_t); }

  double vmag() const { // velocity magnitude
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }

  std::ostream& print(std::ostream& os, const int cpdims, const std::vector<std::string>& scalar_components) const {
    if (cpdims == 2) os << "x=(" << x[0] << ", " << x[1] << "), ";
    else os << "x=(" << x[0] << ", " << x[1] << ", " << x[2] << "), ";
    os << "t=" << t << ", ";

    for (int k = 0; k < scalar_components.size(); k ++)
      os << scalar_components[k] << "=" << scalar[k] << ", ";

    os << "v=";
    if (cpdims == 2) os << "(" << v[0] << ", " << v[1] << "), ";
    else os << "(" << v[0] << ", " << v[1] << ", " << v[2] << "), ";
    
    os << "type=" << critical_point_type_to_string(cpdims, type, scalar_components.size()) << ", "; 
    os << "timestep=" << timestep << ", ";
    os << "ordinal=" << ordinal << ", ";
    os << "tag=" << tag << ", "; 
    os << "id=" << id;  // << std::endl;
    return os;
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
      j["v"] = cp.v;
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
      cp.v = j["v"];
      cp.type = j["type"];
      cp.ordinal = j["ordinal"];
      cp.tag = j["tag"];
      cp.id = j["id"];
    }
  };
}

// serialization
namespace diy {
  static void save(diy::BinaryBuffer& bb, const ftk::critical_point_t &cp) {
    diy::save(bb, cp.x); // for (int i = 0; i < 3; i ++) diy::save(bb, cp.x[i]);
    diy::save(bb, cp.t);
    diy::save(bb, cp.timestep);
    diy::save(bb, cp.scalar); // for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++) diy::save(bb, cp.scalar[i]);
    diy::save(bb, cp.v);
    diy::save(bb, cp.type);
    diy::save(bb, cp.ordinal);
    diy::save(bb, cp.tag);
    diy::save(bb, cp.id);
  }

  static void load(diy::BinaryBuffer& bb, ftk::critical_point_t &cp) {
    diy::load(bb, cp.x); // for (int i = 0; i < 4; i ++) diy::load(bb, cp.x[i]);
    diy::load(bb, cp.t);  
    diy::load(bb, cp.timestep);
    diy::load(bb, cp.scalar); // for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++) diy::load(bb, cp.scalar[i]);
    diy::load(bb, cp.v);
    diy::load(bb, cp.type);
    diy::load(bb, cp.ordinal);
    diy::load(bb, cp.tag);
    diy::load(bb, cp.id);
  }
}

#endif
