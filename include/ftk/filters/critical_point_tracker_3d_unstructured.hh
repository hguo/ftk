#ifndef _FTK_CRITICAL_POINT_TRACKER_3D_UNSTRUCTURED_HH
#define _FTK_CRITICAL_POINT_TRACKER_3D_UNSTRUCTURED_HH

#include <ftk/config.hh>
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
#include <ftk/filters/unstructured_3d_tracker.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_3d_mesh.hh>
#include <ftk/external/diy/serialization.hpp>

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

namespace ftk {

// typedef critical_point_t<3, double> critical_point_t;

struct critical_point_tracker_3d_unstructured : 
  public critical_point_tracker, public unstructured_3d_tracker
{
  // critical_point_tracker_3d_unstructured(const simplicial_unstructured_extruded_3d_mesh<>& m) : m(m) {}
  // critical_point_tracker_3d_unstructured() {}
  critical_point_tracker_3d_unstructured(diy::mpi::communicator comm, const simplicial_unstructured_3d_mesh<>& m) : 
    critical_point_tracker(comm), unstructured_3d_tracker(comm, m), tracker(comm) {}
  virtual ~critical_point_tracker_3d_unstructured() {};
  
  // int cpdims() const { return 2; }

  void initialize() {}
  void finalize();
  void reset() {}

  void update_timestep();

  void push_scalar_field_snapshot(const ndarray<double>&) {} // TODO
  void push_vector_field_snapshot(const ndarray<double>&) {} // TODO

public:
  std::vector<feature_point_t> get_critical_points() const;
  void put_critical_points(const std::vector<feature_point_t>&);

protected:
  bool check_simplex(int, feature_point_t& cp);

  template <int n, typename T> void simplex_values(
      const int verts[n], // vertices
      T X[n][4], // coordinates
      T f[n][FTK_CP_MAX_NUM_VARS], // scalars
      T v[n][3], // vectors
      T J[n][3][3] // jacobians
  ) const;

protected:
  std::map<int, feature_point_t> discrete_critical_points;
  // std::vector<std::vector<critical_point_t>> traced_critical_points;
};

////////////////////////

template <int n, typename T>
inline void critical_point_tracker_3d_unstructured::simplex_values(
    const int verts[n], T X[n][4], T f[n][FTK_CP_MAX_NUM_VARS], T v[n][3], T J[n][3][3]) const
{
  for (int i = 0; i < n; i ++) {
    const int iv = m.flat_vertex_time(verts[i]) == current_timestep ? 0 : 1;
    const int k = m.flat_vertex_id(verts[i]);
    const auto &data = field_data_snapshots[iv];
    m.get_coords(verts[i], X[i]);

    if (!data.scalar.empty())
      for (int j = 0; j < get_num_scalar_components(); j ++)
        f[i][j] = data.scalar(j, k);

    for (int j = 0; j < 3; j ++) 
      v[i][j] = data.vector(j, k);

    if (!data.jacobian.empty())
      for (int j = 0; j < 3; j ++) 
        for (int j1 = 0; j1 < 3; j1 ++)
          J[i][j][j1] = data.jacobian(j1, j, k);
  }
}

inline bool critical_point_tracker_3d_unstructured::check_simplex(int i, feature_point_t& cp)
{
  int tet[4];
  m.get_simplex(3, i, tet); 

  double X[4][4], f[4][FTK_CP_MAX_NUM_VARS], V[4][3], Js[4][3][3];
  simplex_values<4, double>(tet, X, f, V, Js);

  // fprintf(stderr, "tet=%d, %d, %d, %d\n", tet[0], tet[1], tet[2], tet[3]);
  // print4x3("V", V);

#if FTK_HAVE_GMP
  typedef mpf_class fp_t;
  fp_t Vf[4][3];
  for (int k = 0; k < 4; k ++) 
    for (int j = 0; j < 3; j ++)
      Vf[k][j] = V[k][j];
#else
  int64_t Vf[4][3];
  for (int k = 0; k < 4; k ++) 
    for (int j = 0; j < 3; j ++)
      Vf[k][j] = V[k][j] * vector_field_scaling_factor;
#endif
   
  bool succ = ftk::robust_critical_point_in_simplex3(Vf, tet);
  if (!succ) return false;

  double mu[4], x[4], cond;
  bool succ2 = ftk::inverse_lerp_s3v3(V, mu, &cond);
  if (!succ2) { // clamp to normal range
    // return false;
    if (std::isnan(mu[0]) || std::isinf(mu[0])) 
      return false;
      // mu[0] = mu[1] = mu[2] = 1.0 / 3;
    else {
      mu[0] = std::max(std::min(1.0, mu[0]), 0.0);
      mu[1] = std::max(std::min(1.0-mu[0], mu[1]), 0.0);
      mu[2] = std::max(std::min(1.0-mu[1]-mu[0], mu[2]), 0.0);
      mu[2] = 1.0 - mu[0] - mu[1] - mu[2];
      // fprintf(stderr,  "mu =%f, %f, %f\n", mu[0], mu[1], mu[2]);
    }
  }
  
  ftk::lerp_s3v4(X, mu, x);
  cp.x[0] = x[0];
  cp.x[1] = x[1];
  cp.x[2] = x[2];
  cp.t = x[3];

  fprintf(stderr, "X=%f, %f, %f, %f, tet=%d, tettype=%d\n", x[0], x[1], x[2], x[3], i, m.tet_type(i));

  if (!field_data_snapshots[0].scalar.empty())
    for (int k = 0; k < get_num_scalar_components(); k ++)
      cp.scalar[k] = f[0][k] * mu[0] + f[1][k] * mu[1] + f[2][k] * mu[2] + f[3][k] * mu[3];

  if (!field_data_snapshots[0].jacobian.empty()) {
    double H[3][3]; // hessian or jacobian
    ftk::lerp_s3m3x3(Js, mu, H);
    cp.type = ftk::critical_point_type_3d(H, true);
  }
  
  cp.tag = i;
  cp.ordinal = m.is_ordinal(3, i);
  cp.timestep = current_timestep;

  return true;
}

inline void critical_point_tracker_3d_unstructured::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

#ifndef FTK_HAVE_GMP
  update_vector_field_scaling_factor();
#endif
  
  auto func = [&](int i) {
    feature_point_t cp;
    if (check_simplex(i, cp)) {
      std::lock_guard<std::mutex> guard(mutex);

      if (enable_ignoring_degenerate_points) {
        if (cp.type == 0 || cp.type == 1) return;
      }

      if (filter_critical_point_type(cp)) {
        // if (discrete_critical_points.find(i) != discrete_critical_points.end())
        //   fprintf(stderr, "FATAL: overwritting cp!!\n");
        discrete_critical_points[i] = cp;
      }
    }
  };

  m.element_for_ordinal(3, current_timestep, func, xl, nthreads, enable_set_affinity);
  // fprintf(stderr, "#dcp=%zu\n", discrete_critical_points.size());
  if (field_data_snapshots.size() >= 2)
    m.element_for_interval(3, current_timestep, func, xl, nthreads, enable_set_affinity);
  // fprintf(stderr, "#dcp=%zu\n", discrete_critical_points.size());

  if (enable_streaming_trajectories) {
    // grow trajectories
    trace_critical_points_online<int>(
        traced_critical_points, 
        discrete_critical_points,
        [&](int f) {
          std::set<int> neighbors;
          const auto cells = m.side_of(3, f);
          for (const auto c : cells)
            for (const auto f1 : m.sides(4, c))
              neighbors.insert(f1);
          return neighbors;
        },
        [](unsigned long long i) {return  i;},
        [](unsigned long long i) {return  i;}
    );
  }
}

inline void critical_point_tracker_3d_unstructured::finalize()
{
  fprintf(stderr, "finalizing...\n");

  if (enable_streaming_trajectories)  {
    // already done
  } else {
    // Convert connected components to geometries
    auto neighbors = [&](int f) {
      std::set<int> neighbors;
      const auto cells = m.side_of(3, f);
#if 0
      fprintf(stderr, "tet=%d\n", f);
      int vf[4];
      m.get_simplex(3, f, vf);
      fprintf(stderr, "tet.simplex=%d, %d, %d, %d\n", vf[0], vf[1], vf[2], vf[3]);
      for (const auto c : cells) 
        fprintf(stderr, "tet.sideof=%d\n", c);
      // fprintf(stderr, "tet.sideof#=%zu\n", cells.size());
#endif
      for (const auto c : cells) {
#if 0
        fprintf(stderr, "--pent=%d\n", c);
        int vc[5];
        m.get_simplex(4, c, vc);
        fprintf(stderr, "--pent.simplex=%d, %d, %d, %d, %d\n", vc[0], vc[1], vc[2], vc[3], vc[4]);
#endif
        const auto elements = m.sides(4, c);
        for (const auto f1 : elements) {
          neighbors.insert(f1);
          // fprintf(stderr, "----tet=%d\n", f1);
        }
      }
      // fprintf(stderr, "size_neighbors=%zu\n", neighbors.size());
      // for (const auto neighbor : neighbors) 
      //   fprintf(stderr, "--neighbor=%d\n", neighbor);
      // assert(neighbors.find(f) != neighbors.end());
      return neighbors;
    };

    fprintf(stderr, "##dcp=%zu\n", discrete_critical_points.size());

    traced_critical_points.add(trace_critical_points_offline<int>(
        discrete_critical_points, neighbors));
    
    fprintf(stderr, "np=%zu, nc=%zu\n", discrete_critical_points.size(), traced_critical_points.size());
  }
  
  if (enable_discarding_interval_points)
    traced_critical_points.foreach([](feature_curve_t& traj) {
      traj.discard_interval_points();
    });

  if (enable_discarding_degenerate_points)
    traced_critical_points.foreach([](feature_curve_t& traj) {
      traj.discard_degenerate_points();
    });
    
  update_traj_statistics();
}

inline std::vector<feature_point_t> critical_point_tracker_3d_unstructured::get_critical_points() const
{
  std::vector<feature_point_t> results;
  for (const auto &kv : discrete_critical_points) 
    results.push_back(kv.second);
  return results;
}

inline void critical_point_tracker_3d_unstructured::put_critical_points(const std::vector<feature_point_t>& data)
{
  // fprintf(stderr, "##cps=%zu\n", data.size());
  for (const auto& cp : data)
    discrete_critical_points[cp.tag] = cp;
}

}

#endif
