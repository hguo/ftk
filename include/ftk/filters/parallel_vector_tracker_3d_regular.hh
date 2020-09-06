#ifndef _FTK_PARALLEL_VECTOR_TRACKER_HH
#define _FTK_PARALLEL_VECTOR_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/parallel_vector_curve_set.hh>
#include <ftk/numeric/parallel_vector_solver3.hh>
#include <ftk/ndarray/grad.hh>
#include <deque>

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

namespace ftk {

struct parallel_vector_tracker_3d_regular : public filter 
{
  parallel_vector_tracker_3d_regular() : m(4) {};
  
  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void push_field_data_snapshot(
      const ndarray<double> &v, 
      const ndarray<double> &w);

  void update_timestep();

  void initialize();
  void update() {}

protected:
  simplicial_regular_mesh m;
  typedef simplicial_regular_mesh_element element_t;
  typedef parallel_vector_point_t pv_t;
  std::multimap<element_t, parallel_vector_point_t> discrete_pvs; // discrete parallel vector points
  
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
  int start_timestep = 0, end_timestep = std::numeric_limits<int>::max();
};

/////
void parallel_vector_tracker_3d_regular::initialize()
{
  // initializing bounds
  m.set_lb_ub({
      static_cast<int>(domain.start(0)),
      static_cast<int>(domain.start(1)),
      static_cast<int>(domain.start(2)),
      start_timestep
    }, {
      static_cast<int>(domain.size(0)),
      static_cast<int>(domain.size(1)),
      static_cast<int>(domain.size(2)),
      end_timestep
    });
}

void parallel_vector_tracker_3d_regular::update_timestep()
{
  fprintf(stderr, "current_timestep = %d\n", current_timestep);

  auto func = [=](element_t e) {
    pv_t pv;
    if (check_simplex(e, pv)) {
      std::lock_guard<std::mutex> guard(mutex);
      discrete_pvs.insert(std::pair<element_t, pv_t>(e, pv));
      fprintf(stderr, "x={%f, %f, %f}, t=%f, lambda=%f\n", pv.x[0], pv.x[1], pv.x[2], pv.t, pv.lambda);
    }
  };
  
  m.element_for_ordinal(2, current_timestep, func);
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
    pv_t& pv)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m);

  double X[3][4], V[3][3], W[3][3];
  simplex_values<3>(e, X, V, W);

  double lambda[3] = {0}, mu[3][3];
  int ns = solve_pv_s2v3(V, W, lambda, mu);

  if (ns>1) fprintf(stderr, "ns=%d, lambda0=%f, lambda1=%f, lambda2=%f\n", 
      ns, lambda[0], lambda[1], lambda[2]);

#if 0
  for (int k = 0; k < std::min(1, ns); k ++) { // only handle the first eigval in this impl
    pv.lambda = lambda[k];
  }
#endif

  return false;
}

}

#endif
