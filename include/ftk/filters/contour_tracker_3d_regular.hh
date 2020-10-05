#ifndef _FTK_CONTOUR_TRACKER_3D_REGULAR_HH
#define _FTK_CONTOUR_TRACKER_3D_REGULAR_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/filters/contour_tracker_regular.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy-ext/gather.hh>

#if FTK_HAVE_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#endif

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

namespace ftk {

struct contour_tracker_3d_regular : public contour_tracker_regular {
  contour_tracker_3d_regular() : contour_tracker_regular(3) {}
  virtual ~contour_tracker_3d_regular() {}

  int cpdims() const { return 3; }

  void initialize();
  void finalize();
  void reset();

  void update_timestep();

protected:
  void write_trajectories_vtk(const std::string& filename)  const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> get_trajectories_vtk() const;
#endif

protected:
  typedef simplicial_regular_mesh_element element_t;
  
protected:
  bool check_simplex(const element_t& s, feature_point_t& cp);
  void trace_intersections();
  void trace_connected_components();

  void simplex_coordinates(const std::vector<std::vector<int>>& vertices, double X[][4]) const;
  void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[]) const;
};


////////////////////
inline void contour_tracker_3d_regular::initialize()
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

inline void contour_tracker_3d_regular::finalize()
{
  diy::mpi::gather(comm, intersections, intersections, get_root_proc());
  return; 
}

inline void contour_tracker_3d_regular::reset()
{
  current_timestep = 0;

  field_data_snapshots.clear();
  intersections.clear();

  contour_tracker::reset();
}

inline bool contour_tracker_3d_regular::check_simplex(
    const element_t& e, feature_point_t& p)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex
  
  int indices[2];
  simplex_indices(vertices, indices);
  
  double f[2];
  simplex_scalars(vertices, f);
 
  const long long factor = 2 << 20;
  long long fi[2];
  for (int i = 0; i < 2; i ++) {
    f[i] = f[i] - threshold;
    fi[i] = f[i] * factor;
  }

  bool succ = robust_critical_point_in_simplex1(f, indices);
  if (!succ) return false;

  double mu[2];
  bool succ2 = inverse_lerp_s1v1(f, mu);
  // if (!succ2) return false;

  double X[2][4], x[4];
  simplex_coordinates(vertices, X);
  lerp_s1v4(X, mu, x);

  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];
  p.t = x[3];
  p.tag = e.to_integer(m);
  p.ordinal = e.is_ordinal(m);
  p.timestep = current_timestep;

  // std::cerr << e << std::endl;
  // print_matrix<double, 2, 4>("X", X);
  // fprintf(stderr, "%f, %f, %f, %f\n", p.x[0], p.x[1], p.x[2], p.t);

  return true;
}

inline void contour_tracker_3d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  auto func = [=](element_t e) {
    feature_point_t p;
    if (check_simplex(e, p)) {
      std::lock_guard<std::mutex> guard(mutex);
      // std::cerr << e << std::endl;

      intersections[e] = p;

#if 1
      auto tris = e.side_of(m);
      for (auto tri : tris) {
        auto tets = tri.side_of(m);
        for (auto tet : tets) {
          if (tet.valid(m)) {
            auto pents = tet.side_of(m);
            for (auto pent : pents) 
              if (pent.valid(m))
                related_cells.insert(pent); 
          }
        }
      }
#endif
    }
  };

  m.element_for_ordinal(1, current_timestep, func, nthreads);
  if (field_data_snapshots.size() >= 2) // interval
    m.element_for_interval(1, current_timestep, current_timestep+1, func, nthreads);
}

inline void contour_tracker_3d_regular::simplex_coordinates(
    const std::vector<std::vector<int>>& vertices, double X[][4]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
}

inline void contour_tracker_3d_regular::simplex_scalars(
    const std::vector<std::vector<int>>& vertices, double values[]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    values[i] = field_data_snapshots[iv].scalar(
        vertices[i][0],
        vertices[i][1],
        vertices[i][2]);
  }
}

#if FTK_HAVE_VTK
inline void contour_tracker_3d_regular::write_trajectories_vtk(const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = 
      vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData( get_trajectories_vtk() );
    writer->Write();
  }
}

vtkSmartPointer<vtkUnstructuredGrid> contour_tracker_3d_regular::get_trajectories_vtk() const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkDataArray> array_time = vtkDoubleArray::New();
 
  array_time->SetName("time");
  array_time->SetNumberOfComponents(1);
  array_time->SetNumberOfTuples(intersections.size());

  auto my_intersections = intersections; // get_intersections();
  unsigned long long i = 0;
  for (auto &kv : my_intersections) {
    double p[3] = {kv.second.x[0], kv.second.x[1], kv.second.x[2]}; // TODO: time
    kv.second.tag = i;
    points->InsertNextPoint(p);
    array_time->SetTuple1(i, kv.second.t);
    i ++;
  }
  grid->SetPoints(points);
  grid->GetPointData()->AddArray(array_time);

  auto add_tet = [&grid](vtkIdType ids[4]) {
    grid->InsertNextCell(VTK_TETRA, 4, ids);
  };

  for (const auto &e : related_cells) { // pentachoron
    int count = 0;
    vtkIdType ids[10]; // 10 edges

    std::set<element_t> unique_edges;
    for (auto tet : e.sides(m)) 
      for (auto tri : tet.sides(m))
        for (auto edge : tri.sides(m))
          unique_edges.insert(edge);

    for (auto edge : unique_edges)
      if (my_intersections.find(edge) != my_intersections.end())
        ids[count ++] = my_intersections[edge].tag;

    if (count == 4) {
      add_tet(ids);
    } else if (count == 6) {

    } else {
      // should not happen
    }
  }

  return grid;
}
#else
inline void contour_tracker_3d_regular::write_trajectories_vtk(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

}

#endif
