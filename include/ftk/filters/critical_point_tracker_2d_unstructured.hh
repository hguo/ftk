#ifndef _FTK_CRITICAL_POINT_TRACKER_2D_UNSTRUCTURED_HH
#define _FTK_CRITICAL_POINT_TRACKER_2D_UNSTRUCTURED_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
// #include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/filters/critical_point_tracker_regular.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_2d_mesh.hh>
#include <ftk/external/diy/serialization.hpp>

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

#if FTK_HAVE_VTK
#include <vtkUnsignedIntArray.h>
#include <vtkVertex.h>
#endif

namespace ftk {

// typedef critical_point_t<3, double> critical_point_t;

struct critical_point_tracker_2d_unstructured : public critical_point_tracker_regular
{
  // critical_point_tracker_2d_unstructured(const simplicial_unstructured_extruded_2d_mesh<>& m) : m(m) {}
  // critical_point_tracker_2d_unstructured() {}
  critical_point_tracker_2d_unstructured(const simplicial_unstructured_2d_mesh<>& m) : m(simplicial_unstructured_extruded_2d_mesh<>(m)) {}
  virtual ~critical_point_tracker_2d_unstructured() {};
  
  int cpdims() const { return 2; }

  void initialize() {}
  void finalize();
  void reset() {}

  void update_timestep();

  void push_scalar_field_snapshot(const ndarray<double>&) {} // TODO
  void push_vector_field_snapshot(const ndarray<double>&) {} // TODO

  std::vector<critical_point_t> get_critical_points() const;

protected:
  bool check_simplex(int, critical_point_t& cp);

  template <int n, typename T> void simplex_values(
      const int verts[n], // vertices
      T X[n][3], // coordinates
      T f[n][FTK_CP_MAX_NUM_VARS], // scalars
      T v[n][2], // vectors
      T J[n][2][2] // jacobians
  ) const;

protected:
  const simplicial_unstructured_extruded_2d_mesh<> m;
  
  std::map<int, critical_point_t> discrete_critical_points;
  // std::vector<std::vector<critical_point_t>> traced_critical_points;
};

////////////////////////

template <int n, typename T>
inline void critical_point_tracker_2d_unstructured::simplex_values(
    const int verts[n], T X[n][3], T f[n][FTK_CP_MAX_NUM_VARS], T v[n][2], T J[n][2][2]) const
{
  for (int i = 0; i < n; i ++) {
    const int iv = m.flat_vertex_time(verts[i]) == current_timestep ? 0 : 1;
    const int k = m.flat_vertex_id(verts[i]);
    const auto &data = field_data_snapshots[iv];
    m.get_coords(verts[i], X[i]);
    if (!data.scalar.empty()) {
      for (int j = 0; j < get_num_scalar_components(); j ++)
        f[i][j] = data.scalar(j, k);
    }
    for (int j = 0; j < 2; j ++) {
      v[i][j] = data.vector(j, k);
      for (int j1 = 0; j1 < 2; j1 ++)
        J[i][j][j1] = data.jacobian(j1, j, k);
    }
  }
}

inline bool critical_point_tracker_2d_unstructured::check_simplex(int i, critical_point_t& cp)
{
#if FTK_HAVE_GMP
  typedef mpf_class fp_t;
#else
  typedef ftk::fixed_point<> fp_t;
#endif
  
  int tri[3];
  m.get_simplex(2, i, tri); 

  double X[3][3], f[3][FTK_CP_MAX_NUM_VARS], V[3][2], Js[3][2][2];
  simplex_values<3, double>(tri, X, f, V, Js);

  fp_t Vf[3][2];
  for (int k = 0; k < 3; k ++) 
    for (int j = 0; j < 2; j ++)
      Vf[k][j] = V[k][j];
   
  bool succ = ftk::robust_critical_point_in_simplex2(Vf, tri);
  if (!succ) return false;

  double mu[3], x[3];
  bool succ2 = ftk::inverse_lerp_s2v2(V, mu);
  if (!succ2) { // clamp to normal range
    fprintf(stderr,  "mu =%f, %f, %f\n", mu[0], mu[1], mu[2]);
    if (std::isnan(mu[0]) || std::isinf(mu[0])) 
      mu[0] = mu[1] = mu[2] = 1.0 / 3;
    else {
      mu[0] = std::max(std::min(1.0, mu[0]), 0.0);
      mu[1] = std::max(std::min(1.0-mu[0], mu[1]), 0.0);
      mu[2] = 1.0 - mu[0] - mu[1];
    }
    fprintf(stderr,  "mu'=%f, %f, %f\n", mu[0], mu[1], mu[2]);
  }
  ftk::lerp_s2v3(X, mu, x);
  cp.x[0] = x[0];
  cp.x[1] = x[1];
  cp.t = x[2];
 
  for (int k = 0; k < get_num_scalar_components(); k ++)
    cp.scalar[k] = f[0][k] * mu[0] + f[1][k] * mu[1] + f[2][k] * mu[2];

  double H[2][2]; // hessian or jacobian
  ftk::lerp_s2m2x2(Js, mu, H);

  cp.type = ftk::critical_point_type_2d(H, true);
  cp.tag = i;
  cp.ordinal = m.is_ordinal(2, i);
  cp.timestep = current_timestep;

  return true;
}

inline void critical_point_tracker_2d_unstructured::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  auto func = [&](int i) {
    critical_point_t cp;
    if (check_simplex(i, cp)) {
      std::lock_guard<std::mutex> guard(mutex);

      if (enable_ignoring_degenerate_points) {
        if (cp.type == 0 || cp.type == 1) return;
      }

      if (filter_critical_point_type(cp)) {
        discrete_critical_points[i] = cp;
        // cp.ordinal = e.is_ordinal(m);
      }
    }
  };

  m.element_for_ordinal(2, current_timestep, func, nthreads);
  if (field_data_snapshots.size() >= 2)
    m.element_for_interval(2, current_timestep, func, nthreads);

  if (enable_streaming_trajectories) {
    // grow trajectories
    grow_trajectories<int>(
        traced_critical_points, 
        discrete_critical_points,
        [&](int f) {
          std::set<int> neighbors;
          const auto cells = m.side_of(2, f);
          for (const auto c : cells)
            for (const auto f1 : m.sides(3, c))
              neighbors.insert(f1);
          return neighbors;
        },
        [](unsigned long long i) {return  i;}
    );
  }
}

inline void critical_point_tracker_2d_unstructured::finalize()
{
  fprintf(stderr, "finalizing...\n");

  if (enable_streaming_trajectories)  {
    // already done
  } else {
    // Convert connected components to geometries
    auto neighbors = [&](int f) {
      std::set<int> neighbors;
      const auto cells = m.side_of(2, f);
      // fprintf(stderr, "face=%d\n", f);
      // int vf[3];
      // m.get_simplex(2, f, vf);
      // fprintf(stderr, "face.simplex=%d, %d, %d\n", vf[0], vf[1], vf[2]);

      for (const auto c : cells) {
        // fprintf(stderr, "--cell=%d\n", c);
        // int vc[4];
        // m.get_simplex(3, c, vc);
        // fprintf(stderr, "--cell.simplex=%d, %d, %d, %d\n", vc[0], vc[1], vc[2], vc[3]);

        const auto elements = m.sides(3, c);
        for (const auto f1 : elements) {
          // fprintf(stderr, "----face=%d\n", f1);
          neighbors.insert(f1);
        }
      }
      // fprintf(stderr, "size_neighbors=%zu\n", neighbors.size());
      // assert(neighbors.find(f) != neighbors.end());
      return neighbors;
    };

    traced_critical_points = trace_critical_points_offline<int>(
        discrete_critical_points, neighbors);
    
    fprintf(stderr, "np=%zu, nc=%zu\n", discrete_critical_points.size(), traced_critical_points.size());
  }
  
  if (enable_discarding_interval_points)
    for (auto& traj : traced_critical_points)
      traj.discard_interval_points();

  if (enable_discarding_degenerate_points)
    for (auto& traj : traced_critical_points)
      traj.discard_degenerate_points();
    
  update_traj_statistics();
}

inline std::vector<critical_point_t> critical_point_tracker_2d_unstructured::get_critical_points() const
{
  std::vector<critical_point_t> results;
  for (const auto &kv : discrete_critical_points) 
    results.push_back(kv.second);
  return results;
}

}

#endif
