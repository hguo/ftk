#ifndef _FTK_RIDGE_VALLEY_TRACKER_3D_REGULAR_HH
#define _FTK_RIDGE_VALLEY_TRACKER_3D_REGULAR_HH

#include <ftk/config.hh>
#include <ftk/filters/sujudi_haimes_tracker_3d_regular.hh>

namespace ftk {

struct ridge_valley_tracker_3d_regular : public sujudi_haimes_tracker_3d_regular 
{
  ridge_valley_tracker_3d_regular(diy::mpi::communicator comm) : 
    sujudi_haimes_tracker_3d_regular(comm),
    critical_line_tracker(comm),
    regular_tracker(comm, 3), 
    tracker(comm) {}
  virtual ~ridge_valley_tracker_3d_regular() {};

protected:
  bool check_simplex(const element_t& s, feature_point_t& cp) const;
  void push_field_data_snapshot(const ndarray<float>& data);
  
  void simplex_scalars(
      const std::vector<std::vector<int>>& vertices,
      float scalars[3]) const;

  std::vector<std::string> varnames() const { return {"residue", "discJ", "scalar"}; }
};

void ridge_valley_tracker_3d_regular::simplex_scalars(
      const std::vector<std::vector<int>>& vertices,
      float scalars[3]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    const auto &f = field_data_snapshots[iv];

    const auto idx = f.scalar.index(std::vector<size_t>({
          vertices[i][0] - local_array_domain.start(0), 
          vertices[i][1] - local_array_domain.start(1), 
          vertices[i][2] - local_array_domain.start(2)}));

    scalars[i] = f.scalar[idx];
  }
}

inline bool ridge_valley_tracker_3d_regular::check_simplex(const element_t& e, feature_point_t& cp) const
{
  if (!sujudi_haimes_tracker_3d_regular::check_simplex(e, cp)) return false;
  
  float mu[3] = {(float)cp.v[0], (float)cp.v[1], (float)cp.v[2]};
  float scalars[3];
  simplex_scalars(e.vertices(m), scalars);
 
  float scalar = lerp_s2(scalars, mu);
  cp.scalar[2] = scalar;

  return true;
}

inline void ridge_valley_tracker_3d_regular::push_field_data_snapshot(const ndarray<float> &data)
{
  field_data_snapshot_t snapshot; 
  // snapshot.uv = data;
 
  snapshot.scalar = data;
  snapshot.v = gradient3D(data);
  snapshot.J = jacobian3D(snapshot.v);
  snapshot.uv = cross_product3D(snapshot.v, Jv_dot_v(snapshot.v));

  field_data_snapshots.emplace_back(snapshot);
}

}

#endif
