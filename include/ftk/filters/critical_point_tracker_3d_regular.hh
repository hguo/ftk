#ifndef _FTK_CRITICAL_POINT_TRACKER_3D_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_3D_REGULAR_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
// #include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/filters/critical_point.hh>
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

#if FTK_HAVE_CUDA
extern std::vector<ftk::critical_point_t> // <4, double>> 
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
  critical_point_tracker_3d_regular() : m(4) {}
  virtual ~critical_point_tracker_3d_regular() {}
  
  int cpdims() const { return 3; }

  void initialize();
  void finalize();

  void update_timestep();
  
  void push_scalar_field_snapshot(const ndarray<double>&);
  void push_vector_field_snapshot(const ndarray<double>&);

  std::vector<critical_point_t> get_critical_points() const;
  
protected:
  simplicial_regular_mesh m;
  
  typedef simplicial_regular_mesh_element element_t;
  
  std::map<element_t, critical_point_t> discrete_critical_points;
  std::vector<std::set<element_t>> connected_components;
  // std::vector<std::vector<critical_point_t>> traced_critical_points;

protected:
  bool check_simplex(const element_t& s, critical_point_t& cp);
  void trace_intersections();
  void trace_connected_components();

  virtual void simplex_positions(const std::vector<std::vector<int>>& vertices, double X[4][4]) const;
  virtual void simplex_vectors(const std::vector<std::vector<int>>& vertices, double v[4][3]) const;
  virtual void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[4]) const;
  virtual void simplex_jacobians(const std::vector<std::vector<int>>& vertices, 
      double Js[4][3][3]) const;
};


////////////////////
void critical_point_tracker_3d_regular::initialize()
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

  if (use_default_domain_partition) {
    lattice_partitioner partitioner(domain);
    
    // a ghost size of 2 is necessary for jacobian derivaition; 
    // even if jacobian is not necessary, a ghost size of 1 is 
    // necessary for accessing values on boundaries
    partitioner.partition(comm.size(), {}, {2, 2, 2});

    local_domain = partitioner.get_core(comm.rank());
    local_array_domain = partitioner.get_ext(comm.rank());
  }

  if (!is_input_array_partial)
    local_array_domain = array_domain;
}

void critical_point_tracker_3d_regular::finalize()
{
  diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, get_root_proc());

  if (comm.rank() == 0) {
    fprintf(stderr, "finalizing...\n");
    // trace_intersections();
    trace_connected_components();
  }
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
  fprintf(stderr, "current_timestep = %d\n", current_timestep);

  // scan 3-simplices
  // fprintf(stderr, "tracking 3D critical points...\n");
  auto func3 = [=](element_t e) {
      critical_point_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        discrete_critical_points[e] = cp;
        fprintf(stderr, "%f, %f, %f, %f, type=%d\n", cp[0], cp[1], cp[2], cp[3], cp.type);
      }
    };

  if (xl == FTK_XL_NONE) {
    m.element_for(3, lattice({ // ordinal
          local_domain.start(0), 
          local_domain.start(1), 
          local_domain.start(2), 
          static_cast<size_t>(current_timestep), 
        }, {
          local_domain.size(0), 
          local_domain.size(1), 
          local_domain.size(2), 
          1
        }), 
        ftk::ELEMENT_SCOPE_ORDINAL, 
        func3, nthreads);

    if (field_data_snapshots.size() >= 2) { // interval
      m.element_for(3, lattice({
            local_domain.start(0), 
            local_domain.start(1), 
            local_domain.start(2), 
            static_cast<size_t>(current_timestep - 1), 
          }, {
            local_domain.size(0), 
            local_domain.size(1), 
            local_domain.size(2), 
            1
          }),
          ftk::ELEMENT_SCOPE_INTERVAL, 
          func3, nthreads);
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
    
    for (auto cp : results) {
      element_t e(4, 3);
      e.from_work_index(m, cp.tag, ordinal_core, ELEMENT_SCOPE_ORDINAL);
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
      fprintf(stderr, "interval_results#=%d\n", results.size());
      for (auto cp : results) {
        element_t e(4, 3);
        e.from_work_index(m, cp.tag, interval_core, ELEMENT_SCOPE_INTERVAL);
        discrete_critical_points[e] = cp;
      }
    }
#else
    assert(false);
#endif
  }
}

void critical_point_tracker_3d_regular::trace_connected_components()
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
      critical_point_traj_t traj; 
      for (int k = 0; k < linear_graphs[j].size(); k ++)
        traj.push_back(discrete_critical_points[linear_graphs[j][k]]);
      traced_critical_points.emplace_back(traj);
    }
  }
}

void critical_point_tracker_3d_regular::simplex_positions(
    const std::vector<std::vector<int>>& vertices, double X[4][4]) const
{
  for (int i = 0; i < 4; i ++)
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
}

void critical_point_tracker_3d_regular::simplex_vectors(
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

void critical_point_tracker_3d_regular::simplex_scalars(
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

void critical_point_tracker_3d_regular::simplex_jacobians(
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


bool critical_point_tracker_3d_regular::check_simplex(
    const simplicial_regular_mesh_element& e,
    critical_point_t& cp)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m);

  double v[4][3]; // vector values on vertices
  simplex_vectors(vertices, v);
  // ftk::print4x3("v", v);

  double mu[4]; // check intersection
  bool succ = ftk::inverse_lerp_s3v3(v, mu);
  if (!succ) return false;
  
  double X[4][4], x[4]; // position
  simplex_positions(vertices, X);
  lerp_s3v4(X, mu, x);
  cp.x[0] = x[0];
  cp.x[1] = x[1];
  cp.x[2] = x[2];
  cp.t = x[3];

  return true; // TODO
 
  if (scalar_field_source != SOURCE_NONE) {
    double values[3];
    simplex_scalars(vertices, values);
    cp.scalar[0] = lerp_s3(values, mu);
  }

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
} 

inline std::vector<critical_point_t> critical_point_tracker_3d_regular::get_critical_points() const
{
  std::vector<critical_point_t> results;
  for (const auto &kv : discrete_critical_points) 
    results.push_back(kv.second);
  return results;
}

}

#endif
