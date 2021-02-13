#ifndef _FTK_SUJUDI_HAIMES_TRACKER_3D_REGULAR_HH
#define _FTK_SUJUDI_HAIMES_TRACKER_3D_REGULAR_HH

#include <ftk/config.hh>
#include <ftk/filters/critical_line_tracker_3d_regular.hh>

namespace ftk {

struct sujudi_haimes_tracker_3d_regular : public critical_line_tracker_3d_regular 
{
  sujudi_haimes_tracker_3d_regular(diy::mpi::communicator comm=MPI_COMM_WORLD) : critical_line_tracker_3d_regular(comm) {}
  virtual ~critical_line_tracker_3d_regular();

  bool check_simplex(const element_t& s, feature_point_t& cp) const;

  void simplex_resudue_J(
      const std::vector<std::vector<int>>& vertices,
      float residues[3], float Js[3][3][3]) const;

  void push_field_data_snapshot(const ndarray<float>& data);
};

void sujudi_haimes_tracker_3d_regular::simplex_resudue_J(
      const std::vector<std::vector<int>>& vertices,
      float residue[3], float Js[3][3][3]) const
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

    const auto idxJ = f.J.index(std::vector<size_t>({
          0, 0,
          vertices[i][0] - local_array_domain.start(0), 
          vertices[i][1] - local_array_domain.start(1), 
          vertices[i][2] - local_array_domain.start(2)}));

    for (int j = 0; j < 3; j ++)
      for (int k = 0; k < 3; k ++)
        Js[i][j][k] = f.J[idxJ + j*3 + k];
  }

}

inline bool sujudi_haimes_tracker_3d_regular::check_simplex(const element_t& e, feature_point_t& cp) const
{
  if (!critical_line_tracker_3d_regular::check_simplex(e, cp)) return false;

  float residues[3], Js[3][3][3];
  
  simplex_resudue_J(e.vertices(m), residues, Js);

  return true;
}

inline void sujudi_haimes_tracker_3d_regular::push_field_data_snapshot(const ndarray<float> &data)
{
  field_data_snapshot_t snapshot; 
  // snapshot.uv = data;

  snapshot.J = jacobian3D(data);
  snapshot.uv = cross_product3D(data, Jv_dot_v(data));

  field_data_snapshots.emplace_back(snapshot);
}

}

#endif
