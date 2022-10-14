#ifndef _FTK_CRITICAL_POINT_TRACKER_3D_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_3D_REGULAR_HH

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
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/features/feature_point.hh>
#include <ftk/filters/critical_point_tracker_regular.hh>
#include <ftk/external/diy/serialization.hpp>

#if FTK_HAVE_VTK
#include <vtkSmartPointer.h>
#include <vtkUnsignedIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>
#endif

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

#if FTK_HAVE_CUDA
extern std::vector<ftk::feature_point_lite_t> // <4, double>> 
extract_cp3dt_cuda(
    int scope, int current_timestep, 
    const ftk::lattice& domain4,
    const ftk::lattice& core4, 
    const ftk::lattice& ext3,
    const double *Vc, // current timestep
    const double *Vl,  // last timestep
    const double *Jc, // jacobian of current timestep
    const double *Jl, // jacobian of last timestep
    const double *Sc, // scalar of current timestep
    const double *Sl  // scalar of last timestep
  );
#endif

namespace ftk {

struct critical_point_tracker_3d_regular : public critical_point_tracker_regular {
  critical_point_tracker_3d_regular(diy::mpi::communicator comm) : critical_point_tracker_regular(comm, 3), tracker(comm) {}
  virtual ~critical_point_tracker_3d_regular() {}
  
  // int cpdims() const { return 3; }

  void finalize();

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
  virtual void simplex_vectors(const std::vector<std::vector<int>>& vertices, double v[4][3]) const;
  virtual void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[4]) const;
  virtual void simplex_jacobians(const std::vector<std::vector<int>>& vertices, 
      double Js[4][3][3]) const;
  
  void put_critical_points(const std::vector<feature_point_t>&);
};


////////////////////
inline void critical_point_tracker_3d_regular::finalize()
{
  double max_accumulated_kernel_time;
  diy::mpi::reduce(comm, accumulated_kernel_time, max_accumulated_kernel_time, get_root_proc(), diy::mpi::maximum<double>());
  if (comm.rank() == get_root_proc())
    fprintf(stderr, "max_accumulated_kernel_time=%f\n", accumulated_kernel_time);
 
  if (enable_streaming_trajectories) {
    // already done
  } else {
    // diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, get_root_proc());

    if (1) { // comm.rank() == 0) {
      // fprintf(stderr, "finalizing...\n");
      // trace_intersections();
      // trace_connected_components();
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
      }));
    }
  }
  
  update_traj_statistics();
}

inline void critical_point_tracker_3d_regular::push_scalar_field_snapshot(const ndarray<double>& s)
{
  field_data_snapshot_t snapshot;
  
  snapshot.scalar = s;
  if (vector_field_source == SOURCE_DERIVED) {
    snapshot.vector = gradient3D(s);
    if (jacobian_field_source == SOURCE_DERIVED)
      snapshot.jacobian = jacobian3D(snapshot.vector);
  }

  field_data_snapshots.emplace_back( snapshot );
}

inline void critical_point_tracker_3d_regular::push_vector_field_snapshot(const ndarray<double>& v)
{
  field_data_snapshot_t snapshot;
 
  snapshot.vector = v;
  if (jacobian_field_source == SOURCE_DERIVED)
    snapshot.jacobian = jacobian3D(snapshot.vector);

  field_data_snapshots.emplace_back( snapshot );
}

inline void critical_point_tracker_3d_regular::update_timestep()
{
  if (comm.rank() == 0) 
    fprintf(stderr, "current_timestep = %d\n", current_timestep);
  
// #ifndef FTK_HAVE_GMP
  update_vector_field_scaling_factor();
// #endif
  
  typedef std::chrono::high_resolution_clock clock_type;
  auto t0 = clock_type::now();

  // scan 3-simplices
  // fprintf(stderr, "tracking 3D critical points...\n");
  auto func3 = [=](element_t e) {
      feature_point_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        discrete_critical_points[e] = cp;
        // fprintf(stderr, "x={%f, %f, %f}, t=%f, cond=%f, type=%d\n", cp[0], cp[1], cp[2], cp.t, cp.cond, cp.type);
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
        }, 
        [&](unsigned long long tag) {
          return element_t(m, 3, tag);
        },
        [&](element_t e) {
          return e.to_integer(m);
        });
  };

  if (xl == FTK_XL_NONE) {
    element_for_ordinal(3, func3);
    if (field_data_snapshots.size() >= 2) { // interval
      element_for_interval(3, func3);
      
      if (enable_streaming_trajectories)
        grow();
    }
  } else if (xl == FTK_XL_CUDA) {
#if FTK_HAVE_CUDA
    ftk::lattice domain4({
          domain.start(0), 
          domain.start(1), 
          domain.start(2), 
          0
        }, {
          domain.size(0)-1,
          domain.size(1)-1,
          domain.size(2)-1,
          std::numeric_limits<int>::max()
        });

    ftk::lattice ordinal_core({
          local_domain.start(0), 
          local_domain.start(1), 
          local_domain.start(2), 
          static_cast<size_t>(current_timestep), 
        }, {
          local_domain.size(0), 
          local_domain.size(1), 
          local_domain.size(2), 
          1
        });

    ftk::lattice interval_core({
          local_domain.start(0), 
          local_domain.start(1), 
          local_domain.start(2), 
          // static_cast<size_t>(current_timestep-1), 
          static_cast<size_t>(current_timestep), 
        }, {
          local_domain.size(0), 
          local_domain.size(1), 
          local_domain.size(2), 
          1
        });

    ftk::lattice ext({0, 0, 0}, 
        {field_data_snapshots[0].vector.dim(1), 
         field_data_snapshots[0].vector.dim(2),
         field_data_snapshots[0].vector.dim(3)});

    // ordinal
    auto results = extract_cp3dt_cuda(
        ELEMENT_SCOPE_ORDINAL, 
        current_timestep, 
        domain4,
        ordinal_core,
        ext,
        field_data_snapshots[0].vector.data(),
        NULL, // V[0].data(),
        field_data_snapshots[0].jacobian.data(),
        NULL, // gradV[0].data(),
        field_data_snapshots[0].scalar.data(),
        NULL // scalar[0].data(),
      );
   
    for (auto lcp : results) {
      feature_point_t cp(lcp);
      element_t e(4, 3);
      e.from_work_index(m, cp.tag, ordinal_core, ELEMENT_SCOPE_ORDINAL);
      cp.tag = e.to_integer(m);
      cp.ordinal = true;
      cp.timestep = current_timestep;
      discrete_critical_points[e] = cp;
    }

    if (field_data_snapshots.size() >= 2) { // interval
      fprintf(stderr, "processing interval %d, %d\n", current_timestep - 1, current_timestep);
      auto results = extract_cp3dt_cuda(
          ELEMENT_SCOPE_INTERVAL, 
          current_timestep,
          domain4,
          interval_core,
          ext,
          field_data_snapshots[0].vector.data(), // current
          field_data_snapshots[1].vector.data(), // next
          field_data_snapshots[0].jacobian.data(), 
          field_data_snapshots[1].jacobian.data(),
          field_data_snapshots[0].scalar.data(),
          field_data_snapshots[0].scalar.data()
        );
      fprintf(stderr, "interval_results#=%zu\n", results.size());
      for (auto lcp : results) {
        feature_point_t cp(lcp);
        element_t e(4, 3);
        e.from_work_index(m, cp.tag, interval_core, ELEMENT_SCOPE_INTERVAL);
        cp.tag = e.to_integer(m);
        cp.ordinal = false;
        cp.timestep = current_timestep;
        discrete_critical_points[e] = cp;
      }
      
      if (enable_streaming_trajectories)
        grow();
    }
#else
    assert(false);
#endif
  }
  
  auto t1 = clock_type::now();
  accumulated_kernel_time += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9;
}

inline void critical_point_tracker_3d_regular::trace_connected_components()
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
    std::vector<std::vector<double>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(component, neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      feature_curve_t traj; 
      for (int k = 0; k < linear_graphs[j].size(); k ++)
        traj.push_back(discrete_critical_points[linear_graphs[j][k]]);
      traced_critical_points.add(traj);
    }
  }
}

inline void critical_point_tracker_3d_regular::simplex_coordinates(
    const std::vector<std::vector<int>>& vertices, double X[][4]) const
{
#if 0 // legacy
  for (int i = 0; i < vertices.size(); i ++)
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
#endif

  if (mode_phys_coords == REGULAR_COORDS_SIMPLE) {
    for (int i = 0; i < vertices.size(); i ++) {
      X[i][0] = vertices[i][0]; // x
      X[i][1] = vertices[i][1]; // y
      X[i][2] = vertices[i][2]; // z
      X[i][3] = vertices[i][3]; // t
    }
  } else if (mode_phys_coords == REGULAR_COORDS_BOUNDS) {
    for (int i = 0; i < vertices.size(); i ++) {
      X[i][0] = ((vertices[i][0] - array_domain.lower_bound(0)) / double(array_domain.size(0)-1)) * (bounds_coords[1] - bounds_coords[0]) + bounds_coords[0] ; // x
      X[i][1] = ((vertices[i][1] - array_domain.lower_bound(1)) / double(array_domain.size(1)-1)) * (bounds_coords[3] - bounds_coords[2]) + bounds_coords[2] ; // y
      X[i][2] = ((vertices[i][2] - array_domain.lower_bound(2)) / double(array_domain.size(2)-1)) * (bounds_coords[5] - bounds_coords[4]) + bounds_coords[4]; // z
      X[i][3] = vertices[i][3]; // t
    }
  } else if (mode_phys_coords == REGULAR_COORDS_RECTILINEAR) {
    for (int i = 0; i < vertices.size(); i ++) {
      X[i][0] = rectilinear_coords[0][ vertices[i][0] ]; // x
      X[i][1] = rectilinear_coords[1][ vertices[i][1] ]; // y
      X[i][2] = rectilinear_coords[2][ vertices[i][2] ]; // z
      X[i][3] = vertices[i][3]; // t
    }
  } else if (mode_phys_coords == REGULAR_COORDS_EXPLICIT) {
    for (int i = 0; i < vertices.size(); i ++) {
      X[i][0] = explicit_coords(0, vertices[i][0], vertices[i][1]); // x
      X[i][1] = explicit_coords(1, vertices[i][0], vertices[i][1]); // y
      X[i][2] = explicit_coords(2, vertices[i][0], vertices[i][1]);
      X[i][3] = vertices[i][2]; // t
    }
  }
}

inline void critical_point_tracker_3d_regular::simplex_vectors(
    const std::vector<std::vector<int>>& vertices, double v[4][3]) const
{
  for (int i = 0; i < 4; i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    for (int j = 0; j < 3; j ++)
      v[i][j] = field_data_snapshots[iv].vector(j, 
          vertices[i][0] - local_array_domain.start(0), 
          vertices[i][1] - local_array_domain.start(1),
          vertices[i][2] - local_array_domain.start(2));
  }
}

inline void critical_point_tracker_3d_regular::simplex_scalars(
    const std::vector<std::vector<int>>& vertices, double values[4]) const
{
  for (int i = 0; i < 4; i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    values[i] = field_data_snapshots[iv].scalar(
        vertices[i][0] - local_array_domain.start(0), 
        vertices[i][1] - local_array_domain.start(1), 
        vertices[i][2] - local_array_domain.start(2));
  }
}

inline void critical_point_tracker_3d_regular::simplex_jacobians(
    const std::vector<std::vector<int>>& vertices, 
    double Js[4][3][3]) const
{
  for (int i = 0; i < 4; i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    for (int j = 0; j < 3; j ++) {
      for (int k = 0; k < 3; k ++) {
        Js[i][j][k] = field_data_snapshots[iv].jacobian(k, j, 
            vertices[i][0] - local_array_domain.start(0), 
            vertices[i][1] - local_array_domain.start(1), 
            vertices[i][2] - local_array_domain.start(2));
      }
    }
  }
}


inline bool critical_point_tracker_3d_regular::check_simplex(
    const simplicial_regular_mesh_element& e,
    feature_point_t& cp)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m);

  double v[4][3]; // vector values on vertices
  simplex_vectors(vertices, v);
  
  double mu[4]; // check intersection
  double cond; // condition number
  bool succ2 = ftk::inverse_lerp_s3v3(v, mu, &cond);
  // if (!succ2) return false; // TODO
  // for (int i = 0; i < 4; i ++)
  //   if (std::isnan(mu[i]) || std::isinf(mu[i])) return false;

  if (enable_robust_detection) {
#if 0 // FTK_HAVE_GMP
    typedef mpf_class fp_t;
    fp_t vf[4][3];
    for (int i = 0; i < 4; i ++)
      for (int j = 0; j < 3; j ++) {
        const double x = v[i][j];
        if (std::isnan(x) || std::isinf(x)) return false;
        else vf[i][j] = v[i][j];
      }
#else
    int64_t vf[4][3];
    for (int i = 0; i < 4; i ++)
      for (int j = 0; j < 3; j ++) {
        const double x = v[i][j];
        if (std::isnan(x) || std::isinf(x)) return false;
        else vf[i][j] = v[i][j] * vector_field_scaling_factor;
      }
#endif

    int indices[4];
    simplex_indices(vertices, indices);
    bool succ = robust_critical_point_in_simplex3(vf, indices);
    if (!succ) return false;
  } else {
    if (!succ2) return false;
  }
      
  clamp_barycentric<4>(mu);

  double X[4][4], x[4]; // position
  simplex_coordinates(vertices, X);
  lerp_s3v4(X, mu, x);
  cp.x[0] = x[0];
  cp.x[1] = x[1];
  cp.x[2] = x[2];
  cp.t = x[3];
  // cp.cond = cond;

  if (scalar_field_source != SOURCE_NONE) {
    double values[4];
    simplex_scalars(vertices, values);
    cp.scalar[0] = lerp_s3(values, mu);
  }

  double Js[4][3][3], J[3][3];
  simplex_jacobians(vertices, Js);
  ftk::lerp_s3m3x3(Js, mu, J);

  cp.type = ftk::critical_point_type_3d(J, is_jacobian_field_symmetric);
  cp.tag = e.to_integer(m);
  cp.ordinal = e.is_ordinal(m);
  cp.timestep = current_timestep;

  // fprintf(stderr, "mu=%f, %f, %f, x=%f, %f, %f, %f\n", 
  //     mu[0], mu[1], mu[2], cp.x[0], cp.x[1], cp.x[2], cp.t);

  return true; // TODO
#if 0 // legacy
  double J[3][3] = {0}; // jacobian or hessian
  if (jacobian_field_source != SOURCE_NONE) {
    double Js[4][3][3];
    simplex_jacobians(vertices, Js);
    ftk::lerp_s3m3x3(Js, mu, J);
  } else {
    // TODO: jacobian not given
  }

  cp.type = critical_point_type_3d(J, is_jacobian_field_symmetric);
  if (filter_critical_point_type(cp)) return true; 
  else return false;
#endif
} 

inline void critical_point_tracker_3d_regular::put_critical_points(const std::vector<feature_point_t>& data) 
{
  for (const auto& cp : data) {
    element_t e(m, 3, cp.tag);
    discrete_critical_points[e] = cp;
  }
}

}

#endif
