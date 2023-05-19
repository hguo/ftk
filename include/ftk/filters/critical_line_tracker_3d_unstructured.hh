#ifndef _FTK_CRITICAL_LINE_TRACKER_3D_UNSTRUCTURED_HH
#define _FTK_CRITICAL_LINE_TRACKER_3D_UNSTRUCTURED_HH

#include <ftk/config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/filters/critical_line_tracker.hh>
#include <ftk/filters/unstructured_3d_tracker.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_3d_mesh.hh>
#include <ftk/external/diy/serialization.hpp>

namespace ftk {

struct critical_line_tracker_3d_unstructured :
  public virtual critical_line_tracker, public virtual unstructured_3d_tracker
{
  critical_line_tracker_3d_unstructured(diy::mpi::communicator comm, std::shared_ptr<simplicial_unstructured_3d_mesh<>> m) : 
    critical_line_tracker(comm), unstructured_3d_tracker(comm, m), tracker(comm) {}
  ~critical_line_tracker_3d_unstructured() {}
  
  void initialize() {}
  void finalize();
  void reset() {}

  void update_timestep();

public:
  bool check_simplex(int, feature_point_t& cp);

  template <int n, typename T>
  void simplex_values(const int verts[n], 
      T x[n][4], 
      T UV[][2]) const;

public:
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_intersections_vtp() const;
#endif

protected:
  std::map<int, feature_point_t> intersections;
  std::set<int> related_cells;
};

////////////

template <int n, typename T>
inline void critical_line_tracker_3d_unstructured::simplex_values(
    const int verts[n], T X[n][4], T uv[n][2]) const
{
  for (int i = 0; i < n; i ++) {
    const int iv = m->flat_vertex_time(verts[i]) == current_timestep ? 0 : 1;
    const int k = m->flat_vertex_id(verts[i]);
    const auto &data = field_data_snapshots[iv];
    m->get_coords(verts[i], X[i]);

#if 0
    if (!data.scalar.empty())
      for (int j = 0; j < get_num_scalar_components(); j ++)
        f[i][j] = data.scalar(j, k);
#endif

    for (int j = 0; j < 2; j ++) 
      uv[i][j] = data.uv(j, k);
  }
}

inline bool critical_line_tracker_3d_unstructured::check_simplex(int i, feature_point_t& p)
{
  int tri[3];
  m->get_simplex(2, i, tri);

  double X[3][4], UV[3][2], uv[3][2];
  simplex_values<3, double>(tri, X, UV);

  // print3x2("UV", UV);

  // locate zero
  double mu[3], // barycentric coordinates
        cond; // condition number

  const long long factor = 2 << 20; // TODO
  long long UVf[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++) {
      UVf[i][j] = factor * UV[i][j];
      uv[i][j] = UV[i][j];
    }
  
  bool succ = robust_critical_point_in_simplex2(UVf, tri);
  if (!succ) return false;

  inverse_lerp_s2v2(uv, mu, &cond);
  // mu[0] = mu[1] = mu[2] = 0.3333;
  clamp_barycentric<3, double>(mu);

  // interpolation
  double x[4];
  lerp_s2v4(X, mu, x);

  // result
  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];
  p.t = x[3];
#if 0
  p.v[0] = mu[0]; // using v channel to store mu for now for the derived classes
  p.v[1] = mu[1];
  p.v[2] = mu[2];
  // p.scalar[0] = mu[0] * UV[0][2] + mu[1] * UV[1][2] + mu[2] * UV[2][2];
  // p.cond = cond;
#endif
  p.tag = i;
  p.ordinal = m->is_ordinal(2, i);
  p.timestep = current_timestep;

  // fprintf(stderr, "%f, %f, %f, %f\n", p[0], p[1], p[2], p[3]);
  return true;
}

inline void critical_line_tracker_3d_unstructured::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  
  auto func = [&](int i) {
    feature_point_t cp;
    if (check_simplex(i, cp)) {
      std::lock_guard<std::mutex> guard(mutex);
      intersections[i] = cp;
    }
  };

  m->element_for_ordinal(2, current_timestep, func, xl, nthreads, enable_set_affinity);
  if (field_data_snapshots.size() >= 2)
    m->element_for_interval(2, current_timestep, func, xl, nthreads, enable_set_affinity);
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> critical_line_tracker_3d_unstructured::get_intersections_vtp() const
{
  vtkSmartPointer<vtkPolyData> poly = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
 
  vtkIdType pid[1];
  for (const auto &kv : intersections) {
    const auto &cp = kv.second;
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]}; 
    // if (cpdims() == 2) p[2] = cp.t;
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  poly->SetPoints(points);
  poly->SetVerts(vertices);
  
  return poly;
}
#endif

inline void critical_line_tracker_3d_unstructured::finalize()
{
  if (comm.rank() == get_root_proc()) {
    fprintf(stderr, "#intersections=%zu, #related_cells=%zu\n", 
        intersections.size(), related_cells.size());

#if 0
    if (start_timestep == current_timestep) // single timestep
      build_vortex_lines();
    else // multi timesteps
      build_vortex_surfaces();
#endif
  }
}

}

#endif
