#ifndef _FTK_CRITICAL_POINT_T_HH
#define _FTK_CRITICAL_POINT_T_HH

#include <ftk/config.hh>
#include <ftk/features/feature_point_lite.hh>
#include <ftk/features/mpas_particle.hh>
// #include <ftk/numeric/critical_point_type.hh>
#include <ftk/external/diy/serialization.hpp>
#include <yaml-cpp/yaml.h>

namespace ftk {

using json = YAML::Node;

/**
 * @brief Represents a detected feature point (critical point, particle, etc.)
 *
 * This structure stores comprehensive information about a single feature point
 * detected in a scientific dataset. Feature points can represent:
 * - Critical points in scalar/vector fields (minima, maxima, saddles)
 * - Particles in flow fields
 * - Intersection points in contour tracking
 * - Blob centroids in threshold-based tracking
 *
 * The structure includes:
 * - Spatial coordinates (x) and time (t)
 * - Scalar values at the feature location
 * - Velocity/motion vectors (v)
 * - Type classification (minimum, maximum, saddle, etc.)
 * - Tracking metadata (tag, id, timestep)
 *
 * Feature points can be serialized to JSON, binary, or text formats and
 * are the building blocks of feature_curve_t trajectories.
 */
struct feature_point_t {
  feature_point_t() {}
  feature_point_t(const feature_point_t& p) {
    x = p.x;
    t = p.t;
    // cond = p.cond;
    timestep = p.timestep;
    scalar = p.scalar;
    v = p.v;
    type = p.type;
    ordinal = p.ordinal;
    tag = p.tag;
    id = p.id;
  }
  feature_point_t& operator=(const feature_point_t& p) {
    x = p.x;
    t = p.t;
    // cond = p.cond;
    timestep = p.timestep;
    scalar = p.scalar;
    v = p.v;
    type = p.type;
    ordinal = p.ordinal;
    tag = p.tag;
    id = p.id;
    return *this;
  }
  feature_point_t(const feature_point_lite_t& cp) {
    for (int i = 0; i < 3; i ++)
      x[i] = cp.x[i];
    t = cp.t;
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++)
      scalar[i] = cp.scalar[i];
    type = cp.type;
    tag = cp.tag;
  }
  feature_point_t(const mpas_particle_t<>& p) {
    double x[3];
    p.get_x(x);
    for (int i = 0; i < 3; i ++) {
      this->x[i] = x[i];
      this->v[i] = p.v[i];
    }
    this->t = p.t;
    this->scalar[0] = p.vv;
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS-1; i ++)
      this->scalar[i+1] = p.scalar[i];
  }

  feature_point_lite_t to_lite() const {
    feature_point_lite_t l;
    for (int i = 0; i < 3; i ++)
      l.x[i] = x[i];
    l.t = t;
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++)
      l.scalar[i] = scalar[i];
    l.type = type;
    l.tag = tag;
    return l;
  }

  mpas_particle_t<> to_mpas() const {
    mpas_particle_t<> p;
    p.set_x( this->x.data() );
    p.t = this->t;
    p.hint_c = this->tag;
    p.hint_l = this->type;
    p.vv = this->scalar[0];
    for (int i = 0; i < FTK_CP_MAX_NUM_VARS; i ++)
      p.scalar[i] = this->scalar[i];

    return p;
  }

  /**
   * @brief Access spatial coordinate by index (const)
   * @param i Coordinate index (0=x, 1=y, 2=z)
   * @return Coordinate value
   */
  double operator[](size_t i) const {return x[i];}

  /**
   * @brief Access spatial coordinate by index (mutable)
   * @param i Coordinate index (0=x, 1=y, 2=z)
   * @return Reference to coordinate value
   */
  double &operator[](size_t i) {return x[i];}

  // constexpr size_t size() const noexcept { return sizeof(feature_point_t); }

  /**
   * @brief Compute velocity magnitude
   * @return Magnitude of the velocity vector
   */
  double vmag() const { // velocity magnitude
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }

  /**
   * @brief Convert Cartesian coordinates to geographic coordinates
   * @param R0 Reference radius (default: Earth radius in meters)
   * @return Tuple of (longitude in degrees, latitude in degrees, altitude)
   */
  std::tuple<double, double, double> lonlatz(const double R0 = 6371229.0) const {
    const double R = std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    return std::make_tuple(
      std::atan2(x[1], x[0]) / M_PI * 180, // longitude
      std::asin(x[2] / R) / M_PI * 180, // latitude
      R - R0);
  }

  /**
   * @brief Print feature point information to stream
   * @param os Output stream
   * @param scalar_components Names of scalar components
   * @return Output stream reference
   */
  std::ostream& print(std::ostream& os, const std::vector<std::string>& scalar_components) const {
    os << "x=(" << x[0] << ", " << x[1] << ", " << x[2] << "), ";
    os << "t=" << t << ", ";

    // os << "cond=" << cond << ", ";

    for (int k = 0; k < scalar_components.size(); k ++)
      os << scalar_components[k] << "=" << scalar[k] << ", ";

    os << "v=";
    os << "(" << v[0] << ", " << v[1] << ", " << v[2] << "), ";
    
    os << "type=" << type << ", "; 
    os << "timestep=" << timestep << ", ";
    os << "ordinal=" << ordinal << ", ";
    os << "tag=" << tag << ", "; 
    os << "id=" << id;  // << std::endl;
    return os;
  }

  /**
   * @brief Compute spatiotemporal distance between two feature points
   * @param a First feature point
   * @param b Second feature point
   * @return Distance in 4D spacetime
   */
  friend double dist(const feature_point_t& a, const feature_point_t& b) {
    return std::exp2(a.x[0] - b.x[0])
      + std::exp2(a.x[1] - b.x[1])
      + std::exp2(a.x[2] - b.x[2])
      + std::exp2(a.t - b.t);
  }

public:
  std::array<double, 3> x = {0}; ///< Spatial coordinates (x, y, z)
  double t = 0.0; ///< Time coordinate
  // double cond = 0.0; // condition number
  int timestep = 0; ///< Discrete timestep index
  // double rx[4] = {0}; // coordinates in transformed (e.g. curvilinear) grid, if eligible
  std::array<double, FTK_CP_MAX_NUM_VARS> scalar = {0}; ///< Scalar field values at the feature point
  std::array<double, 3> v = {0}; ///< Velocity or motion vector
  unsigned int type = 0; ///< Feature type classification (e.g., min, max, saddle)
  bool ordinal = false; ///< True if point occurs exactly at a timestep (not interpolated)
  unsigned long long tag = 0; ///< Mesh element tag/identifier
  unsigned long long id = 0; ///< Trajectory ID that this point belongs to
};

}

// serialization w/ yaml-cpp
namespace YAML {
  using namespace ftk;

  template<>
  struct convert<feature_point_t> {
    static Node encode(const feature_point_t& cp) {
      Node node;
      node["x"] = cp.x;
      node["t"] = cp.t;
      // node["cond"] = cp.cond;
      node["timestep"] = cp.timestep;
      node["scalar"] = cp.scalar;
      node["v"] = cp.v;
      node["type"] = cp.type;
      node["ordinal"] = cp.ordinal;
      node["tag"] = cp.tag;
      node["id"] = cp.id;
      return node;
    }

    static bool decode(const Node& node, feature_point_t& cp) {
      if(!node.IsMap()) return false;
      auto x_vec = node["x"].as<std::vector<double>>();
      std::copy_n(x_vec.begin(), std::min(x_vec.size(), cp.x.size()), cp.x.begin());
      cp.t = node["t"].as<double>();
      // cp.cond = node["cond"].as<bool>();
      cp.timestep = node["timestep"].as<int>();
      auto scalar_vec = node["scalar"].as<std::vector<double>>();
      std::copy_n(scalar_vec.begin(), std::min(scalar_vec.size(), cp.scalar.size()), cp.scalar.begin());
      auto v_vec = node["v"].as<std::vector<double>>();
      std::copy_n(v_vec.begin(), std::min(v_vec.size(), cp.v.size()), cp.v.begin());
      cp.type = node["type"].as<int>();
      cp.ordinal = node["ordinal"].as<size_t>();
      cp.tag = node["tag"].as<unsigned long long>();
      cp.id = node["id"].as<unsigned long long>();
      return true;
    }
  };
}

// serialization
namespace diy {
  template <> struct Serialization<ftk::feature_point_t> {
    static void save(diy::BinaryBuffer& bb, const ftk::feature_point_t &cp) {
      diy::save(bb, cp.x); 
      diy::save(bb, cp.t);
      // diy::save(bb, cp.cond);
      diy::save(bb, cp.timestep);
      diy::save(bb, cp.scalar); 
      diy::save(bb, cp.v);
      diy::save(bb, cp.type);
      diy::save(bb, cp.ordinal);
      diy::save(bb, cp.tag);
      diy::save(bb, cp.id);
    }

    static void load(diy::BinaryBuffer& bb, ftk::feature_point_t &cp) {
      diy::load(bb, cp.x); 
      diy::load(bb, cp.t);  
      // diy::load(bb, cp.cond);
      diy::load(bb, cp.timestep);
      diy::load(bb, cp.scalar); 
      diy::load(bb, cp.v);
      diy::load(bb, cp.type);
      diy::load(bb, cp.ordinal);
      diy::load(bb, cp.tag);
      diy::load(bb, cp.id);
    }
  };
}

#endif
