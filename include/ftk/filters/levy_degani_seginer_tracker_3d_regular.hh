#ifndef _FTK_LEVY_DEGANI_SEGINER_TRACKER_3D_REGULAR_HH
#define _FTK_LEVY_DEGANI_SEGINER_TRACKER_3D_REGULAR_HH

#include <ftk/config.hh>
#include <ftk/filters/critical_line_tracker_3d_regular.hh>

namespace ftk {

struct levy_degani_seginer_tracker_3d_regular : public critical_line_tracker_3d_regular 
{
  levy_degani_seginer_tracker_3d_regular(diy::mpi::communicator comm) : 
    critical_line_tracker_3d_regular(comm), 
    critical_line_tracker(comm),
    regular_tracker(comm, 3), 
    tracker(comm) {}
  virtual ~levy_degani_seginer_tracker_3d_regular() {};

protected:
  bool check_simplex(const element_t& s, feature_point_t& cp) const;

  void simplex_residue_vorts(
      const std::vector<std::vector<int>>& vertices,
      float residues[3], float vorts[3][3]) const;

  void push_field_data_snapshot(const ndarray<float>& data);

  std::vector<std::string> varnames() const { return {"residue", "vortmag"}; }
};

void levy_degani_seginer_tracker_3d_regular::simplex_residue_vorts(
      const std::vector<std::vector<int>>& vertices,
      float residues[3], float vorts[3][3]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    const auto &f = field_data_snapshots[iv];

    const auto idx = f.uv.index(std::vector<size_t>({
          0,
          vertices[i][0] - local_array_domain.start(0), 
          vertices[i][1] - local_array_domain.start(1), 
          vertices[i][2] - local_array_domain.start(2)}));

    residues[i] = f.uv[idx+2];
    for (int j = 0; j < 3; j ++)
      vorts[i][j] = f.vorticity[idx+j];
  }

}

inline bool levy_degani_seginer_tracker_3d_regular::check_simplex(const element_t& e, feature_point_t& cp) const
{
  if (!critical_line_tracker_3d_regular::check_simplex(e, cp)) return false;
  
  float mu[3] = {(float)cp.v[0], (float)cp.v[1], (float)cp.v[2]};
  float residues[3], vorts[3][3];
  simplex_residue_vorts(e.vertices(m), residues, vorts);
 
  float residue = lerp_s2(residues, mu);
  cp.scalar[0] = residue;

  float vort[3];
  lerp_s2v3(vorts, mu, vort);
  cp.scalar[1] = vector_2norm_2(vort);

  // fprintf(stderr, "residue=%f, detJ=%f\n", residue, detJ);
  // print3x3("J", J);

  return true;
}

inline void levy_degani_seginer_tracker_3d_regular::push_field_data_snapshot(const ndarray<float> &data)
{
  field_data_snapshot_t snapshot; 
  // snapshot.uv = data;
  
  snapshot.vorticity = vorticity3D(data);
  snapshot.uv = cross_product3D(data, snapshot.vorticity);

  field_data_snapshots.emplace_back(snapshot);
}

}

#endif
