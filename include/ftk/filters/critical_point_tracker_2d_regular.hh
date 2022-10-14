#ifndef _FTK_CRITICAL_POINT_TRACKER_2D_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_2D_REGULAR_HH

#include <ftk/config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/clamp.hh>
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
#include <ftk/numeric/critical_point_degree.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/filters/critical_point_tracker_regular.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/utils/gather.hh>

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

extern std::vector<ftk::feature_point_lite_t> // <3, double>> 
extract_cp2dt_cuda(
    int scope, int current_timestep, 
    const ftk::lattice& domain, // 3D
    const ftk::lattice& core, // 3D
    const ftk::lattice& ext, // 2D, array dimension
    const double *Vc, // current timestep
    const double *Vn, // next timestep
    const double *Jc, // jacobian of current timestep
    const double *Jn, // jacobian of next timestep
    const double *Sc, // scalar of current timestep
    const double *Sn, // scalar of next timestep
    bool use_explicit_coords,
    const double *coords // coords of vertices
  ); 

extern std::vector<ftk::feature_point_lite_t> // <3, double>> 
extract_cp2dt_sycl(
    int scope, int current_timestep, 
    const ftk::lattice& domain, // 3D
    const ftk::lattice& core, // 3D
    const ftk::lattice& ext, // 2D, array dimension
    const double *Vc, // current timestep
    const double *Vn, // next timestep
    const double *Jc, // jacobian of current timestep
    const double *Jn, // jacobian of next timestep
    const double *Sc, // scalar of current timestep
    const double *Sn, // scalar of next timestep
    bool use_explicit_coords,
    const double *coords // coords of vertices
  ); 

static std::vector<ftk::feature_point_lite_t>
extract_cp2dt_xl_wrapper(
    int xl,
    int scope, int current_timestep, 
    const ftk::lattice& domain, // 3D
    const ftk::lattice& core, // 3D
    const ftk::lattice& ext, // 2D, array dimension
    const double *Vc, // current timestep
    const double *Vn, // next timestep
    const double *Jc, // jacobian of current timestep
    const double *Jn, // jacobian of next timestep
    const double *Sc, // scalar of current timestep
    const double *Sn, // scalar of next timestep
    bool use_explicit_coords,
    const double *coords // coords of vertices
  )
{
  using namespace ftk;
  if (xl == FTK_XL_CUDA) {
#if FTK_HAVE_CUDA
    return extract_cp2dt_cuda(scope, current_timestep, domain, core, ext, Vc, Vn, Jc, Jn, Sc, Sn, use_explicit_coords, coords);
#else
    fatal(FTK_ERR_NOT_BUILT_WITH_CUDA);
    return std::vector<ftk::feature_point_lite_t>();
#endif
  } else if (xl == FTK_XL_SYCL) {
#if FTK_HAVE_SYCL
    return extract_cp2dt_sycl(scope, current_timestep, domain, core, ext, Vc, Vn, Jc, Jn, Sc, Sn, use_explicit_coords, coords);
#else
    fatal(FTK_ERR_NOT_BUILT_WITH_HIPSYCL);
    return std::vector<ftk::feature_point_lite_t>();
#endif
  } else {
    fatal(FTK_ERR_ACCELERATOR_UNSUPPORTED);
    return std::vector<ftk::feature_point_lite_t>();
  }
}

namespace ftk {

// typedef critical_point_t<3, double> critical_point_t;

struct critical_point_tracker_2d_regular : public critical_point_tracker_regular {
  critical_point_tracker_2d_regular(diy::mpi::communicator comm) : 
    critical_point_tracker_regular(comm, 2),
    tracker(comm)
  {}
  virtual ~critical_point_tracker_2d_regular() {}

  int cpdims() const { return 2; }

  void finalize();
  void reset();

  void update_timestep();

  void push_scalar_field_snapshot(const ndarray<double>&);
  void push_vector_field_snapshot(const ndarray<double>&);

protected:
  typedef simplicial_regular_mesh_element element_t;
  
protected:
  bool check_simplex(const element_t& s, feature_point_t& cp);
  void trace_intersections();
  void trace_connected_components();

  virtual void simplex_coordinates(const std::vector<std::vector<int>>& vertices, double X[][4]) const;
  template <typename T=double> void simplex_vectors(const std::vector<std::vector<int>>& vertices, T v[][2]) const;
  virtual void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[]) const;
  virtual void simplex_jacobians(const std::vector<std::vector<int>>& vertices, 
      double Js[][2][2]) const;
  
  void put_critical_points(const std::vector<feature_point_t>&);
};


////////////////////
inline void critical_point_tracker_2d_regular::finalize()
{
  double max_accumulated_kernel_time;
  diy::mpi::reduce(comm, accumulated_kernel_time, max_accumulated_kernel_time, get_root_proc(), diy::mpi::maximum<double>());
  if (comm.rank() == get_root_proc())
    fprintf(stderr, "max_accumulated_kernel_time=%f\n", accumulated_kernel_time);

  if (enable_streaming_trajectories) {
    // done
  } else {
    // fprintf(stderr, "rank=%d, root=%d, #cp=%zu\n", comm.rank(), get_root_proc(), discrete_critical_points.size());
    // diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, get_root_proc());

    if (1) { // if (comm.rank() == get_root_proc()) {
      // fprintf(stderr, "finalizing...\n");
#if 0
      duf<uint64_t> uf(comm);

      for (const auto &kv : discrete_critical_points) {
        std::set<element_t> neighbors;
        const auto cells = kv.first.side_of(m);
        for (const auto c : cells) {
          const auto elements = c.sides(m);
          for (const auto f1 : elements)
            neighbors.insert(f1);
        }

        for (const auto &n : neighbors)
          // if (true) // (discrete_critical_points.find(n) != discrete_critical_points.end())
          if (kv.first != n) {
            uint64_t i0 = kv.first.to_integer<uint64_t>(m), 
                     i1 = n.to_integer<uint64_t>(m);
            // fprintf(stderr, "uniting %lld, %lld, rank=%d\n", i0, i1, comm.rank());
            uf.unite(kv.first.to_integer<uint64_t>(m), n.to_integer<uint64_t>(m));
          }
      }
      // fprintf(stderr, "dUF sync...\n");
      uf.sync();
      // fprintf(stderr, "dUF done.\n");
      fprintf(stderr, "dUF done., #pts=%zu, #roots=%zu\n", discrete_critical_points.size(), uf.get_roots().size());

      comm.barrier();
      exit(1);
#endif

      traced_critical_points.add( trace_critical_points_offline<element_t>(discrete_critical_points, 
          [&](element_t f) {
            std::set<element_t> neighbors;
            const auto cells = f.side_of(m);
            for (const auto c : cells) {
              const auto elements = c.sides(m);
              for (const auto f1 : elements)
                neighbors.insert(f1);
            }
            return neighbors;
#if 0
            std::set<element_t> strict_neighbors;
            if (discrete_critical_points.find(f) != discrete_critical_points.end()) {
              for (const auto neighbor : neighbors) {
                if (discrete_critical_points.find(neighbor) != discrete_critical_points.end()) {
                  if (discrete_critical_points.at(neighbor).type == discrete_critical_points.at(f).type)
                    strict_neighbors.insert(neighbor);
                }
              }
            }

            fprintf(stderr, "size_strict_neighbors=%ld\n", strict_neighbors.size());
            return strict_neighbors;
#endif
      }));

      // trace_intersections();
      // trace_connected_components();
    }
  }
  
  if (enable_discarding_interval_points)
    traced_critical_points.foreach([](feature_curve_t& traj) {
      traj.discard_interval_points();
    });

  update_traj_statistics();
}

inline void critical_point_tracker_2d_regular::reset()
{
  current_timestep = 0;

  field_data_snapshots.clear();
  discrete_critical_points.clear();
  traced_critical_points.clear();

  critical_point_tracker::reset();
}

inline void critical_point_tracker_2d_regular::push_scalar_field_snapshot(const ndarray<double>& s)
{
  field_data_snapshot_t snapshot;

  snapshot.scalar = s;
  if (vector_field_source == SOURCE_DERIVED) {
    snapshot.vector = gradient2D(s);
    if (jacobian_field_source == SOURCE_DERIVED)
      snapshot.jacobian = jacobian2D<double, true>(snapshot.vector);
  }

  field_data_snapshots.emplace_back( snapshot );
}

inline void critical_point_tracker_2d_regular::push_vector_field_snapshot(const ndarray<double>& v)
{
  field_data_snapshot_t snapshot;
 
  snapshot.vector = v;
  if (jacobian_field_source == SOURCE_DERIVED)
    snapshot.jacobian = jacobian2D(snapshot.vector);

  field_data_snapshots.emplace_back( snapshot );
}

inline void critical_point_tracker_2d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

#ifndef FTK_HAVE_GMP
  update_vector_field_scaling_factor();
#endif

  typedef std::chrono::high_resolution_clock clock_type;
  auto t0 = clock_type::now();

  // scan 2-simplices
  // fprintf(stderr, "tracking 2D critical points...\n");
  auto func2 = [=](element_t e) {
      feature_point_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        if (filter_critical_point_type(cp)) {
          discrete_critical_points[e] = cp;
          // std::cerr << "tag=" << cp.tag << ", " << e << "\t" << element_t(m, 2, cp.tag) << std::endl;
          // assert(element_t(m, 2, cp.tag) == e);
        }
      }
  };

  auto grow = [&]() {
    trace_critical_points_online<element_t>(
        traced_critical_points, 
        discrete_critical_points, 
        [&](element_t f) {
          std::set<element_t> neighbors;
          const auto cells = f.side_of(m);
          for (const auto c : cells) {
            const auto elements = c.sides(m);
            for (const auto f1 : elements)
              neighbors.insert(f1);
          }
          return neighbors;
#if 0
          std::set<element_t> strict_neighbors;
          if (discrete_critical_points.find(f) != discrete_critical_points.end()) {
            for (const auto neighbor : neighbors) {
              if (discrete_critical_points.find(neighbor) != discrete_critical_points.end()) {
                if (discrete_critical_points.at(neighbor).type == discrete_critical_points.at(f).type)
                  strict_neighbors.insert(neighbor);
              }
            }
          }

          fprintf(stderr, "size_strict_neighbors=%ld\n", strict_neighbors.size());
          return strict_neighbors;
#endif
        }, 
        [&](unsigned long long tag) {
          return element_t(m, 2, tag);
        }, 
        [&](element_t e) {
          return e.to_integer(m);
        });
  };

  if (xl == FTK_XL_NONE) {
    element_for_ordinal(2, func2);
    if (field_data_snapshots.size() >= 2) { // interval
      element_for_interval(2, func2);

      if (enable_streaming_trajectories)
        grow();
    }
  } else { //  if (xl == FTK_XL_CUDA) {
    ftk::lattice domain3({
          domain.start(0), 
          domain.start(1), 
          0
        }, {
          domain.size(0)-1,
          domain.size(1)-1,
          std::numeric_limits<int>::max()
        });

    ftk::lattice ordinal_core({
          local_domain.start(0), 
          local_domain.start(1), 
          static_cast<size_t>(current_timestep), 
        }, {
          local_domain.size(0), 
          local_domain.size(1), 
          1
        });

    ftk::lattice interval_core({
          local_domain.start(0), 
          local_domain.start(1), 
          // static_cast<size_t>(current_timestep-1), 
          static_cast<size_t>(current_timestep), 
        }, {
          local_domain.size(0), 
          local_domain.size(1), 
          1
        });

    ftk::lattice ext({0, 0}, 
        {field_data_snapshots[0].vector.dim(1), 
        field_data_snapshots[0].vector.dim(2)});
    
    // ordinal
    auto results = extract_cp2dt_xl_wrapper(
        xl,
        ELEMENT_SCOPE_ORDINAL, 
        current_timestep, 
        domain3,
        ordinal_core,
        ext,
        field_data_snapshots[0].vector.data(),
        NULL, // V[0].data(),
        field_data_snapshots[0].jacobian.data(),
        NULL, // gradV[0].data(),
        field_data_snapshots[0].scalar.data(),
        NULL, // scalar[0].data(),
        use_explicit_coords, 
        coords.data()
      );
    
    fprintf(stderr, "ordinal_results#=%zu\n", results.size());
    for (auto lcp : results) {
      feature_point_t cp(lcp);
      element_t e(3, 2);
      e.from_work_index(m, cp.tag, ordinal_core, ELEMENT_SCOPE_ORDINAL);
      cp.tag = e.to_integer(m);
      cp.ordinal = true;
      cp.timestep = current_timestep;
      discrete_critical_points[e] = cp;
    }

    if (field_data_snapshots.size() >= 2) { // interval
      fprintf(stderr, "processing interval %d, %d\n", current_timestep, current_timestep+1);
      auto results = extract_cp2dt_xl_wrapper(
          xl,
          ELEMENT_SCOPE_INTERVAL, 
          current_timestep,
          domain3, 
          interval_core,
          ext,
          field_data_snapshots[0].vector.data(), // current
          field_data_snapshots[1].vector.data(), // next
          field_data_snapshots[0].jacobian.data(), 
          field_data_snapshots[1].jacobian.data(),
          field_data_snapshots[0].scalar.data(),
          field_data_snapshots[1].scalar.data(),
          use_explicit_coords, 
          coords.data()
        );
      fprintf(stderr, "interal_results#=%zu\n", results.size());
      for (auto lcp : results) {
        feature_point_t cp(lcp);
        element_t e(3, 2);
        e.from_work_index(m, cp.tag, interval_core, ELEMENT_SCOPE_INTERVAL);
        cp.tag = e.to_integer(m);
        cp.ordinal = false;
        cp.timestep = current_timestep;
        discrete_critical_points[e] = cp;
      }
      
      if (enable_streaming_trajectories)
        grow();
    }
  }

  auto t1 = clock_type::now();
  accumulated_kernel_time += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9;
}

inline void critical_point_tracker_2d_regular::trace_intersections()
{
  // scan 3-simplices to get connected components
  union_find<element_t> uf;
  for (const auto &kv : discrete_critical_points) 
    uf.add(kv.first);

  m.element_for(3, [&](const simplicial_regular_mesh_element& f) {
    const auto sides = f.sides(m);
    std::set<element_t> intersected_sides;

    for (const auto& side : sides)
      if (uf.has(side)) 
        intersected_sides.insert(side);

    if (intersected_sides.size() > 1) {// the size of intersected_size should be 0 or 2
      for (auto it = std::next(intersected_sides.begin(), 1); it != intersected_sides.end(); it ++) {
        std::lock_guard<std::mutex> guard(mutex);
        uf.unite(*intersected_sides.begin(), *it);
      }
    }
  });
  uf.get_sets(connected_components);
}

inline void critical_point_tracker_2d_regular::trace_connected_components()
{
  // Convert connected components to geometries
  auto neighbors = [&](element_t f) {
    std::set<element_t> neighbors;
    const auto cells = f.side_of(m);
    for (const auto c : cells) {
      const auto elements = c.sides(m);
      for (const auto f1 : elements)
        neighbors.insert(f1);
    }
    return neighbors;
  };

  std::set<element_t> elements;
  for (const auto &kv : discrete_critical_points)
    elements.insert(kv.first);
  connected_components = extract_connected_components<element_t, std::set<element_t>>(
      neighbors, elements);

  for (const auto &component : connected_components) {
    // std::vector<std::vector<double>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(component, neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      feature_curve_t traj; 
      for (int k = 0; k < linear_graphs[j].size(); k ++)
        traj.push_back(discrete_critical_points[linear_graphs[j][k]]);
      traced_critical_points.add(traj);
      // const auto subtrajs = traj.to_consistent_sub_traj();
      // traced_critical_points.insert(traced_critical_points.end(), subtrajs.begin(), subtrajs.end());
    }
  }
}

inline void critical_point_tracker_2d_regular::simplex_coordinates(
    const std::vector<std::vector<int>>& vertices, double X[][4]) const
{
  if (mode_phys_coords == REGULAR_COORDS_SIMPLE) {
    for (int i = 0; i < vertices.size(); i ++) {
      X[i][0] = vertices[i][0]; // x
      X[i][1] = vertices[i][1]; // y
      X[i][2] = 0.0; // z
      X[i][3] = vertices[i][2]; // t
    }
  } else if (mode_phys_coords == REGULAR_COORDS_BOUNDS) {
    for (int i = 0; i < vertices.size(); i ++) {
      X[i][0] = ((vertices[i][0] - array_domain.lower_bound(0)) / double(array_domain.size(0)-1)) * (bounds_coords[1] - bounds_coords[0]) + bounds_coords[0] ; // x
      X[i][1] = ((vertices[i][1] - array_domain.lower_bound(1)) / double(array_domain.size(1)-1)) * (bounds_coords[3] - bounds_coords[2]) + bounds_coords[2] ; // y
      X[i][2] = 0.0; // z
      X[i][3] = vertices[i][2]; // t
    }
  } else if (mode_phys_coords == REGULAR_COORDS_RECTILINEAR) {
    for (int i = 0; i < vertices.size(); i ++) {
      X[i][0] = rectilinear_coords[0][ vertices[i][0] ]; // x
      X[i][1] = rectilinear_coords[1][ vertices[i][1] ]; // y
      X[i][2] = 0.0; // z
      X[i][3] = vertices[i][2]; // t
    }
  } else if (mode_phys_coords == REGULAR_COORDS_EXPLICIT) {
    for (int i = 0; i < vertices.size(); i ++) {
      X[i][0] = explicit_coords(0, vertices[i][0], vertices[i][1]); // x
      X[i][1] = explicit_coords(1, vertices[i][0], vertices[i][1]); // y
      if (explicit_coords.dim(0) > 2) // z
        X[i][2] = explicit_coords(2, vertices[i][0], vertices[i][1]);
      else 
        X[i][2] = 0.0;
      X[i][3] = vertices[i][2]; // t
    }
  }
#if 0
  if (use_explicit_coords) {
    for (int i = 0; i < vertices.size(); i ++) {
      for (int j = 0; j < 2; j ++) 
        X[i][j] = coords(j, vertices[i][0], vertices[i][1]);
      X[i][2] = vertices[i][2];
    }
  } else {
    for (int i = 0; i < vertices.size(); i ++)
      for (int j = 0; j < 3; j ++)
        X[i][j] = vertices[i][j];
  }
#endif
}

template <typename T>
inline void critical_point_tracker_2d_regular::simplex_vectors(
    const std::vector<std::vector<int>>& vertices, T v[][2]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][2] == current_timestep ? 0 : 1;
    for (int j = 0; j < 2; j ++)
      v[i][j] = field_data_snapshots[iv].vector(j, 
          vertices[i][0] - local_array_domain.start(0), 
          vertices[i][1] - local_array_domain.start(1));
  }
}

inline void critical_point_tracker_2d_regular::simplex_scalars(
    const std::vector<std::vector<int>>& vertices, double values[]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][2] == current_timestep ? 0 : 1;
    values[i] = field_data_snapshots[iv].scalar(
        vertices[i][0] - local_array_domain.start(0), 
        vertices[i][1] - local_array_domain.start(1));
  }
}

inline void critical_point_tracker_2d_regular::simplex_jacobians(
    const std::vector<std::vector<int>>& vertices, 
    double Js[][2][2]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][2] == current_timestep ? 0 : 1;
    for (int j = 0; j < 2; j ++) {
      for (int k = 0; k < 2; k ++) {
        Js[i][j][k] = field_data_snapshots[iv].jacobian(k, j, 
            vertices[i][0] - local_array_domain.start(0), 
            vertices[i][1] - local_array_domain.start(1));
      }
    }
  }
}

inline bool critical_point_tracker_2d_regular::check_simplex(
    const simplicial_regular_mesh_element& e,
    feature_point_t& cp)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex
 
  double v[3][2]; // obtain vector values
  simplex_vectors(vertices, v);

#if FTK_HAVE_GMP
  typedef mpf_class fp_t;
  // typedef double fp_t;

  fp_t vf[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++) {
      const double x = v[i][j];
      if (std::isnan(x) || std::isinf(x)) return false;
      else vf[i][j] = v[i][j];
    }
#else
  // typedef fixed_point<> fp_t;
  const double factor = 1.0 / vector_field_resolution;
  // fprintf(stderr, "factor=%f\n", factor);
  int64_t vf[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++) {
      const double x = v[i][j];
      if (std::isnan(x) || std::isinf(x)) return false;
      else vf[i][j] = v[i][j] * vector_field_scaling_factor;
    }
#endif

  // robust critical point test
  int indices[3];
  simplex_indices(vertices, indices);
  bool succ = robust_critical_point_in_simplex2(vf, indices);
  if (!succ) return false;

  double mu[3]; // check intersection
  double cond;
  bool succ2 = inverse_lerp_s2v2(v, mu, &cond);
  // if (!succ2) return false;
  // if (std::isnan(mu[0]) || std::isnan(mu[1]) || std::isnan(mu[2])) return false;
  // fprintf(stderr, "mu=%f, %f, %f\n", mu[0], mu[1], mu[2]);

  if (!succ2) clamp_barycentric<3>(mu);

  double X[3][4], x[4]; // position
  simplex_coordinates(vertices, X);
  lerp_s2v4(X, mu, x);
  cp.x[0] = x[0];
  cp.x[1] = x[1];
  cp.x[2] = x[2];
  cp.t = x[3];
  // cp.cond = cond;
  // fprintf(stderr, "x=%f, %f, %f\n", cp.x[0], cp.x[1], cp.x[2]);

  if (scalar_field_source != SOURCE_NONE) {
    double values[3];
    simplex_scalars(vertices, values);
    cp.scalar[0] = lerp_s2(values, mu);
  }

  cp.tag = e.to_integer(m);
  cp.ordinal = e.is_ordinal(m);
  cp.timestep = current_timestep;

  if (enable_computing_degrees) {
    if (cp.ordinal) {
      auto deg = positive2(vf, indices);
      int chi = e.type == 4 ? 1 : -1;
      deg *= chi;
      // fprintf(stderr, "deg=%d, chi=%d\n", deg, chi);
      if (deg == 1) cp.type = 1; 
      else cp.type = 2;
    } else 
      cp.type = 0;
  } else {
    double J[2][2] = {0}; // jacobian
    if (jacobian_field_source != SOURCE_NONE) { // lerp jacobian
      double Js[3][2][2];
      simplex_jacobians(vertices, Js);
      lerp_s2m2x2(Js, mu, J);
      ftk::make_symmetric2x2(J); // TODO
    } else {
      // TODO: jacobian is not given
    }
#if 0
    double X2[3][2];
    for (int i = 0; i < 3; i ++)
      for (int j = 0; j < 2; j ++)
        X2[i][j] = X[i][j];
    jacobian_3dsimplex2(X2, v, J);
    ftk::make_symmetric2x2(J); // TODO
#endif
    cp.type = critical_point_type_2d(J, is_jacobian_field_symmetric);
  }

  return true;
} 

inline void critical_point_tracker_2d_regular::put_critical_points(const std::vector<feature_point_t>& data) 
{
  for (const auto& cp : data) {
    element_t e(m, 2, cp.tag);
    discrete_critical_points[e] = cp;
  }
}


}

#endif
