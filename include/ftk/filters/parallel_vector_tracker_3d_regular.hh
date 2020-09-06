#ifndef _FTK_PARALLEL_VECTOR_TRACKER_HH
#define _FTK_PARALLEL_VECTOR_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/parallel_vector_curve_set.hh>
#include <deque>

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

namespace ftk {

struct parallel_vector_tracker_3d_regular : public filter 
{
  parallel_vector_tracker_3d_regular();
  
  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void push_field_data_snapshot(
      const ndarray<double> &v, 
      const ndarray<double> &w);

  void update_timestep();

protected:
  simplicial_regular_mesh m;
  typedef simplicial_regular_mesh_element element_t;
  typedef parallel_vector_point_t pv_t;
  std::map<element_t, parallel_vector_point_t> discrete_pvs; // discrete parallel vector points
  
  struct field_data_snapshot_t {
    ndarray<double> v,  w;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;

  template <int n, typename T>
  void simplex_values(const element_t &e, 
      T X[n][4], T V[n][3], T W[n][3]) const;

private:
  bool check_simplex(const element_t& s, pv_t& cp);

protected:
  lattice domain, array_domain;
  int current_timestep = 0;
};

/////
void parallel_vector_tracker_3d_regular::update_timestep()
{
  fprintf(stderr, "current_timestep = %d\n", current_timestep);

  auto func = [=](element_t e) {
    pv_t pv;
    if (check_simplex(e, pv)) {
      std::lock_guard<std::mutex> guard(mutex);
      discrete_pvs[e] = pv;
      fprintf(stderr, "x={%f, %f, %f}, t=%f, lambda=%f\n", pv.x[0], pv.x[1], pv.x[2], pv.t, pv.lambda);
    }
  };
}

inline void parallel_vector_tracker_3d_regular::push_field_data_snapshot(
    const ndarray<double>& v,
    const ndarray<double>& w)
{
  field_data_snapshot_t snapshot;
  snapshot.v = v;
  snapshot.w = w;

  field_data_snapshots.emplace_back(snapshot);
}

template <int n/*number of vertices*/, typename T>
void parallel_vector_tracker_3d_regular::simplex_values(
    const element_t &e, 
    T X[n][4], T V[n][3], T W[n][3]) const
{
  const auto &vertices = e.vertices(m);
  assert(n == vertices.size());
  
  for (int i = 0; i < n; i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    const auto &data = field_data_snapshots[iv];
   
    // v and w
    for (int j = 0; j < 3; j ++) {
      V[i][j] = field_data_snapshots[iv].v(j, vertices[i][0], vertices[i][1], vertices[i][2]);
      W[i][j] = field_data_snapshots[iv].w(j, vertices[i][0], vertices[i][1], vertices[i][2]);
    }
    
    // coordinates
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
  }
}

bool parallel_vector_tracker_3d_regular::check_simplex(
    const simplicial_regular_mesh_element& e, // 2-simplex
    pv_t& cp)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m);

  double X[3][4], V[3][3], W[3][3];
  simplex_values<3>(e, X, V, W);

  return false;
}

}

#endif
