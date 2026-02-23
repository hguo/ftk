#ifndef _FTK_FEATURE_CURVE_T_HH
#define _FTK_FEATURE_CURVE_T_HH

#include <ftk/features/feature_point.hh>

namespace ftk {

/**
 * @brief Represents a feature trajectory (curve) in spacetime
 *
 * A feature curve is a sequence of feature points that traces the temporal
 * evolution of a feature (e.g., a critical point moving through space over time).
 * It inherits from std::vector<feature_point_t> and adds trajectory-specific
 * metadata and operations.
 *
 * Feature curves can represent:
 * - Critical point trajectories showing feature motion over time
 * - Particle paths in flow fields
 * - Evolution of blob centroids
 * - Loops (periodic features) or open curves
 *
 * The class provides methods for:
 * - Statistical analysis (min/max values, persistence, bounding boxes)
 * - Simplification (discarding interval points, removing noise)
 * - Type smoothing (ensuring consistent feature classification)
 * - Temporal slicing and filtering
 * - Geographic coordinate conversion
 */
struct feature_curve_t : public std::vector<feature_point_t>
{
  feature_curve_t() {}
  feature_curve_t(const feature_curve_t&);
  feature_curve_t& operator=(const feature_curve_t&);

  /**
   * @brief Assign a new ID to the trajectory and all its points
   * @param i New trajectory ID
   */
  void relabel(int i);

  /**
   * @brief Compute trajectory statistics (min/max, persistence, bounding boxes)
   */
  void update_statistics();

  /**
   * @brief Derive velocity from position changes between ordinal points
   */
  void derive_velocity();

  /**
   * @brief Convert Cartesian coordinates to geographic coordinates
   */
  void convert_to_geo();

  /**
   * @brief Split trajectory at longitude discontinuities (for geographic data)
   * @return Vector of sub-trajectories
   */
  std::vector<feature_curve_t> split_geo() const;

  /**
   * @brief Split trajectory into segments with consistent type
   * @return Vector of sub-trajectories, each with a single consistent type
   */
  std::vector<feature_curve_t> split() const;

  /**
   * @brief Get indices of ordinal (non-interpolated) points
   * @return Vector of indices into the original trajectory
   */
  std::vector<int> to_ordinals() const;

  /**
   * @brief Select points matching a predicate
   * @param f Predicate function
   * @return Indices of matching points
   */
  std::vector<int> select(std::function<bool(const feature_point_t&)> f);

  /**
   * @brief Rotate loop trajectory to avoid type inconsistency at endpoints
   */
  void rotate();

  /**
   * @brief Reorder trajectory points in ascending time order
   */
  void reorder();

  /**
   * @brief Ensure time increases monotonically by adjusting interval points
   */
  void adjust_time();

  /**
   * @brief Smooth type classifications of ordinal points using a sliding window
   * @param half_window_size Half-width of the smoothing window (default: 2)
   */
  void smooth_ordinal_types(const int half_window_size=2);

  /**
   * @brief Make type classifications consistent within temporal intervals
   */
  void smooth_interval_types();

  /**
   * @brief Remove points matching a predicate
   * @param f Predicate function (returns true for points to discard)
   */
  void discard(std::function<bool(const feature_point_t&)> f);

  /**
   * @brief Remove all interval (interpolated) points, keeping only ordinal points
   */
  void discard_interval_points();

  /**
   * @brief Remove degenerate critical points (type 0 or 1)
   */
  void discard_degenerate_points();

  /**
   * @brief Locate the interval containing a given time
   * @param t Time value
   * @param cap If true, cap to upper bound; otherwise lower bound
   * @return Interval index
   */
  int locate(double t, bool cap=false) const;

  /**
   * @brief Copy trajectory metadata from another trajectory
   * @param t Source trajectory
   */
  void copy_info(feature_curve_t& t) {id = t.id; complete = t.complete; loop = t.loop;}

  /**
   * @brief Extract sub-trajectory within a time range
   * @param t0 Start time
   * @param t1 End time
   * @return Sub-trajectory containing points in [t0, t1]
   */
  feature_curve_t intercept(int t0, int t1) const;

  /**
   * @brief Unwrap trajectory in a periodic domain
   * @tparam k Coordinate dimension to unwrap (default: 2 for z)
   * @param period Period of the domain
   */
  template <int k=2> void unwrap(const double period);

public:
  int id; ///< Trajectory identifier
  bool complete = false; ///< True if trajectory is complete (reached end of time)
  bool loop = false; ///< True if trajectory forms a closed loop
  std::array<double, FTK_CP_MAX_NUM_VARS> max; ///< Maximum scalar values along trajectory
  std::array<double, FTK_CP_MAX_NUM_VARS> min; ///< Minimum scalar values along trajectory
  std::array<double, FTK_CP_MAX_NUM_VARS> persistence; ///< Persistence (max - min) for each scalar
  std::array<double, 3> bbmin; ///< Spatial bounding box minimum
  std::array<double, 3> bbmax; ///< Spatial bounding box maximum
  double tmin; ///< Minimum time value
  double tmax; ///< Maximum time value
  double vmmin; ///< Minimum velocity magnitude
  double vmmax; ///< Maximum velocity magnitude
  unsigned int consistent_type = 0; ///< Type if all points have same type, 0 otherwise
};

/////////

inline feature_curve_t::feature_curve_t(const feature_curve_t& c)
{
  id = c.id;
  complete = c.complete;
  loop = c.loop;
  max = c.max;
  min = c.min;
  persistence = c.persistence;
  bbmin = c.bbmin;
  bbmax = c.bbmax;
  tmin = c.tmin;
  tmax = c.tmax;
  vmmin = c.vmmin;
  vmmax = c.vmmax;
  consistent_type = c.consistent_type;
  for (const feature_point_t& p : c) 
    this->push_back(p);
}
  
inline feature_curve_t& feature_curve_t::operator=(const feature_curve_t& c)
{
  id = c.id;
  complete = c.complete;
  loop = c.loop;
  max = c.max;
  min = c.min;
  persistence = c.persistence;
  bbmin = c.bbmin;
  bbmax = c.bbmax;
  tmin = c.tmin;
  tmax = c.tmax;
  vmmin = c.vmmin;
  vmmax = c.vmmax;
  consistent_type = c.consistent_type;
  clear();
  for (const feature_point_t& p : c) 
    this->push_back(p);

  return *this;
}

inline void feature_curve_t::relabel(int i)
{
  id = i;
  for (i = 0; i < size(); i ++)
    at(i).id = id;
}

inline std::vector<int> feature_curve_t::select(std::function<bool(const feature_point_t&)> f)
{
  std::vector<int> result;
  for (int i = 0; i < size(); i ++) 
    if (f(at(i)))
      result.push_back(i);
  return result;
}

inline void feature_curve_t::discard(std::function<bool(const feature_point_t&)> f)
{
  feature_curve_t traj;
  traj.copy_info(*this);

  for (const auto &cp : *this)
    if (!f(cp))
      traj.push_back(cp);
  
  traj.update_statistics();
  *this = traj;
}

inline void feature_curve_t::discard_interval_points() {
  discard([&](const feature_point_t& cp) { return !cp.ordinal; });
}
  
inline void feature_curve_t::discard_degenerate_points() {
  discard([&](const feature_point_t& cp) { 
    return cp.type == 0 || cp.type == 1; 
  });
}

// inline void feature_curve_t::discard_high_cond(double threshold)
// {
//   discard([&](const feature_point_t& cp) {
//     if (std::isinf(cp.cond) || std::isnan(cp.cond) || cp.cond > threshold) return true;
//     else return false;
//   });
// }
  
inline void feature_curve_t::update_statistics() {
  if (empty()) return; // nothing to do

  max.fill( std::numeric_limits<double>::lowest() );
  min.fill( std::numeric_limits<double>::max() );
  bbmax.fill( std::numeric_limits<double>::lowest() );
  bbmin.fill( std::numeric_limits<double>::max() );
  tmax = std::numeric_limits<double>::lowest();
  tmin = std::numeric_limits<double>::max();
  vmmax = std::numeric_limits<double>::lowest();
  vmmin = std::numeric_limits<double>::max();

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
    vmmax = std::max(vmmax, at(i).vmag());
    vmmin = std::min(vmmin, at(i).vmag());
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

inline void feature_curve_t::convert_to_geo()
{
  for (size_t i = 0; i < size(); i ++) {
    auto &p = at(i);
    auto lonlatz = p.lonlatz();

    p.x[0] = std::get<0>(lonlatz);
    p.x[1] = std::get<1>(lonlatz);
    p.x[2] = std::get<2>(lonlatz);
  }
}

inline std::vector<feature_curve_t> feature_curve_t::split_geo() const
{
  std::vector<feature_curve_t> results;
 
  feature_curve_t subtraj;
  subtraj.id = this->id;

  for (size_t i = 0; i < size(); i ++) {
    subtraj.push_back(at(i));

    if (i < size() - 1) {
      if (std::abs(at(i+1).x[0] - at(i).x[0]) > 180.0) {
        results.push_back(subtraj);
        subtraj.clear();
      }
    } else { // the last one
      if (subtraj.size())
        results.push_back(subtraj);
    }
  }
  
  // fprintf(stderr, "split from %zu to %zu\n", this->size(), results.size());
  return results;
}

inline std::vector<feature_curve_t> feature_curve_t::split() const 
{
  std::vector<feature_curve_t> results;
  feature_curve_t subtraj;
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
  
inline std::vector<int> feature_curve_t::to_ordinals() const 
{
  std::vector<int> result;
  for (auto i = 0; i < size(); i ++)
    if (at(i).ordinal)
      result.push_back(i);
  return result;
}
  
inline void feature_curve_t::rotate() //! if the traj is a loop and the type is inconsistent, rotate the traj before splitting
{ 
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

inline void feature_curve_t::reorder() { // reorder by timestep
  if (empty()) return;
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
  
inline void feature_curve_t::adjust_time() { //! assuming traj has only one single branch and is organized in ascending order, ensure that the time increases monotonously
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
  
inline void feature_curve_t::smooth_ordinal_types(const int half_window_size) {
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
  
inline void feature_curve_t::smooth_interval_types() { //! make types in intervals consistent
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
  
inline int feature_curve_t::locate(double t, bool cap) const 
{
  return 0; // TODO: WIP
}
  
inline feature_curve_t feature_curve_t::intercept(int t0, int t1) const //! assuming the traj is already reordered
{
  feature_curve_t result;
  if (t0 > back().t || t1 < front().t) 
    return result;

  for (int i = 0; i < size(); i ++)
    if (at(i).t >= t0 && at(i).t <= t1)
      result.push_back(at(i));

  if (result.size() > 0) {
    result.relabel(id);
    result.update_statistics();
  }

  return result;
}
  
template <int k> void feature_curve_t::unwrap(const double period)
{
  for (int i = 0; i < size()-1; i ++) {
    const double delta = at(i+1).x[k] - at(i).x[k];
    if (delta >= period/2) {
      for (int j = i+1; j < size(); j ++)
        at(j).x[k] -= period;
    } else if (delta <= -period/2) {
      for (int j = i+1; j < size(); j ++)
        at(j).x[k] += period;
    }
  }
  
  // update_statistics();
  bbmin[k] = std::numeric_limits<double>::max();
  bbmax[k] = std::numeric_limits<double>::lowest();
  for (auto i = 0; i < size(); i ++) {
    bbmax[k] = std::max(bbmax[k], at(i).x[k]);
    bbmin[k] = std::min(bbmin[k], at(i).x[k]);
  }
}

inline void feature_curve_t::derive_velocity() // const std::vector<double> &kernel)
{
  if (size() < 2) return;
  // const int half_kernel_size = kernel.size() / 2;
  for (int k = 0; k < 3 /*cpdims*/; k ++) {
    for (int i = 0; i < size(); i ++) {
      if (i == 0) at(i).
        v[k] = at(i+1).x[k] - at(i).x[k];
      else if (i == size()-1) 
        at(i).v[k] = at(i).x[k] - at(i-1).x[k];
      else at(i).v[k] = 
        0.5 * (at(i+1).x[k] - at(i-1).x[k]);
    }
#if 0
    // preparation & padding
    std::vector<double> f(size() + half_kernel_size * 2);
    for (int i = 0; i < half_kernel_size; i ++) {
      f[i] = front().x[k];
      f[half_kernel_size + size() + i] = back().x[k];
    }
    for (int i = 0; i < size(); i ++)
      f[half_kernel_size + i] = at(i).x[k];

    // convolution
    for (int i = 0; i < size(); i ++) {
      double v = 0;
      for (int j = 0; j < kernel.size(); j ++)
        v += f[i + half_kernel_size + j] * kernel[j];
      at(i).velocity[k] = v;
    }
#endif
  }
}

} // namespace ftk


// serialization w/ yaml-cpp
namespace YAML {
  using namespace ftk;
  template <>
  struct convert<feature_curve_t> {
    static Node encode(const feature_curve_t& t) {
      Node node;
      node["id"] = t.id;
      node["max"] = std::vector<double>(t.max.begin(), t.max.end());
      node["min"] = std::vector<double>(t.min.begin(), t.min.end());
      node["persistence"] = std::vector<double>(t.persistence.begin(), t.persistence.end());
      node["bbmin"] = std::vector<double>(t.bbmin.begin(), t.bbmin.end());
      node["bbmax"] = std::vector<double>(t.bbmax.begin(), t.bbmax.end());
      node["tmin"] = t.tmin;
      node["tmax"] = t.tmax;
      node["consistent_type"] = t.consistent_type;
      for (const auto& p : static_cast<std::vector<feature_point_t>>(t)) {
        node["traj"].push_back(p);
      }
      return node;
    }

    static bool decode(const Node& node, feature_curve_t& t) {
      if(!node.IsMap()) return false;
      t.id = node["id"].as<size_t>();
      auto max_vec = node["max"].as<std::vector<double>>();
      std::copy_n(max_vec.begin(), std::min(max_vec.size(), t.max.size()), t.max.begin());
      auto min_vec = node["min"].as<std::vector<double>>();
      std::copy_n(min_vec.begin(), std::min(min_vec.size(), t.min.size()), t.min.begin());
      auto persistence_vec = node["persistence"].as<std::vector<double>>();
      std::copy_n(persistence_vec.begin(), std::min(persistence_vec.size(), t.persistence.size()), t.persistence.begin());
      auto bbmin_vec = node["bbmin"].as<std::vector<double>>();
      std::copy_n(bbmin_vec.begin(), std::min(bbmin_vec.size(), t.bbmin.size()), t.bbmin.begin());
      auto bbmax_vec = node["bbmax"].as<std::vector<double>>();
      std::copy_n(bbmax_vec.begin(), std::min(bbmax_vec.size(), t.bbmax.size()), t.bbmax.begin());
      t.tmin = node["tmin"].as<double>();
      t.tmax = node["tmax"].as<double>();
      t.consistent_type = node["consistent_type"].as<int>();
      std::vector<feature_point_t> traj = node["traj"].as<std::vector<feature_point_t>>();
      t.clear();
      t.insert(t.begin(), traj.begin(), traj.end());
      return true;
    }
  };
}

// serialization
namespace diy {
  template <> struct Serialization<ftk::feature_curve_t> {
    static void save(diy::BinaryBuffer& bb, const ftk::feature_curve_t &t) {
      diy::save(bb, t.complete);
      diy::save(bb, t.max);
      diy::save(bb, t.min);
      diy::save(bb, t.persistence);
      diy::save(bb, t.bbmin);
      diy::save(bb, t.bbmax);
      diy::save(bb, t.tmin);
      diy::save(bb, t.tmax);
      diy::save(bb, t.consistent_type);
      diy::save(bb, t.size());
      for (auto i = 0; i < t.size(); i ++)
        diy::save(bb, t[i]);
      // diy::save<std::vector<ftk::feature_point_t>>(bb, t);
    }
    
    static void load(diy::BinaryBuffer& bb, ftk::feature_curve_t &t) {
      diy::load(bb, t.complete);
      diy::load(bb, t.max);
      diy::load(bb, t.min);
      diy::load(bb, t.persistence);
      diy::load(bb, t.bbmin);
      diy::load(bb, t.bbmax);
      diy::load(bb, t.tmin);
      diy::load(bb, t.tmax);
      diy::load(bb, t.consistent_type);
      size_t s;
      diy::load(bb, s);
      t.resize(s);
      for (auto i = 0; i < t.size(); i ++)
        diy::load(bb, t[i]);
      // diy::load<std::vector<ftk::feature_point_t>>(bb, t);
    }
  };
} // namespace diy

#endif
