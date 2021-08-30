#ifndef _FTK_CONTOUR_TRACKER_2D_REGULAR_HH
#define _FTK_CONTOUR_TRACKER_2D_REGULAR_HH

#include <ftk/config.hh>
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
#include <ftk/utils/gather.hh>

#if FTK_HAVE_VTK
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkPolyDataNormals.h>
// #include <vtkPCANormalEstimation.h>
#endif

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

namespace ftk {

struct contour_tracker_2d_regular : public contour_tracker_regular {
  contour_tracker_2d_regular(diy::mpi::communicator comm) : 
    contour_tracker_regular(comm, 2),
    tracker(comm)
  {}
  virtual ~contour_tracker_2d_regular() {}

  // int cpdims() const { return 2; }

  void finalize();
  void reset();

  void update_timestep();

  const feature_surface_t& get_surfaces() const { return surfaces; }

public:
  void write_isovolume_vtu(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_isovolume_vtp() const;
  vtkSmartPointer<vtkPolyData> get_isovolume_vtu() const;
#endif

protected:
  typedef simplicial_regular_mesh_element element_t;

  feature_surface_t surfaces;

protected:
  bool check_simplex(const element_t& s, feature_point_t& cp);
  void trace_intersections();
  void trace_connected_components();

  void build_surfaces();

  virtual void simplex_coordinates(const std::vector<std::vector<int>>& vertices, double X[][3]) const;
  virtual void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[]) const;
};


////////////////////
inline void contour_tracker_2d_regular::build_surfaces()
{
  fprintf(stderr, "building spacetime isosurfaces...\n");
  
  std::mutex my_mutex;
  auto add_triangle = [&](int i0, int i1, int i2) {
    std::lock_guard<std::mutex> guard(my_mutex);
    surfaces.tris.push_back({i0, i1, i2});
  };
  
  auto my_intersections = intersections; // get_intersections();
  unsigned long long i = 0;
  for (auto &kv : my_intersections) {
    double p[3] = {kv.second.x[0], kv.second.x[1], kv.second.t};
    kv.second.tag = i ++;
    surfaces.pts.push_back(kv.second);
  }

  parallel_for<element_t>(related_cells, [&](const element_t &e) {
    int count = 0;
    unsigned long long ids[4];

    std::set<element_t> unique_edges;
    for (auto tri : e.sides(m)) {
      for (auto edge : tri.sides(m)) {
        unique_edges.insert(edge);
      }
    }
    
    for (auto edge : unique_edges) 
      if (my_intersections.find(edge) != my_intersections.end())
        ids[count ++] = my_intersections[edge].tag;

    if (count == 3) {
      add_triangle(ids[0], ids[1], ids[2]);
    } else if (count == 4) { // quad
      add_triangle(ids[0], ids[1], ids[2]);
      add_triangle(ids[1], ids[3], ids[2]);
    } else {
      // fprintf(stderr, "irregular count=%d\n", count); // WIP: triangulation
    }
  }, FTK_THREAD_PTHREAD, nthreads, enable_set_affinity);

  surfaces.relabel();
  fprintf(stderr, "#pts=%zu, #tri=%zu\n", surfaces.pts.size(), surfaces.tris.size());
}

inline void contour_tracker_2d_regular::finalize()
{
  diy::mpi::gather(comm, intersections, intersections, get_root_proc());
  diy::mpi::gather(comm, related_cells, related_cells, get_root_proc());

  if (comm.rank() == get_root_proc())
    build_surfaces();

#if 0
  // TODO: build surfaces
  for (const auto &e : related_cells) {
    auto sides = e.sides(m);
    int count = 0;

    std::set<element_t> unique_edges;
    for (auto tri : e.sides(m)) {
      for (auto edge : tri.sides(m)) {
        unique_edges.insert(edge);
      }
    }
    for (auto edge : unique_edges) 
      if (intersections.find(edge) != intersections.end())
        count ++;
    std::cerr << "count=" << count << ", " << e << std::endl;
  }
#endif
}

inline void contour_tracker_2d_regular::reset()
{
  surfaces.clear();
  contour_tracker_regular::reset();
}

inline bool contour_tracker_2d_regular::check_simplex(
    const element_t& e, feature_point_t& p)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex
  
  int indices[3];
  simplex_indices(vertices, indices);
  
  double f[2];
  simplex_scalars(vertices, f);
 
  const long long factor = 2 << 20;
  long long fi[2];
  for (int i = 0; i < 2; i ++) {
    f[i] = f[i] - threshold;
    fi[i] = f[i] * factor;
  }

  bool succ = robust_critical_point_in_simplex1(fi, indices);
  if (!succ) return false;

  double mu[2];
  bool succ2 = inverse_lerp_s1v1(f, mu);
  // if (!succ2) return false;

  double X[2][3], x[3];
  simplex_coordinates(vertices, X);
  lerp_s1v3(X, mu, x);

  p.x[0] = x[0];
  p.x[1] = x[1];
  p.t = x[2];
  p.tag = e.to_integer(m);
  p.ordinal = e.is_ordinal(m);
  p.timestep = current_timestep;

  // p.print(std::cerr, 2, scalar_components) << std::endl;
  // std::cerr << e << ", (" << f[0] << ", " << f[1] << "), mu=" << mu[0] << std::endl;

  return true;
}

inline void contour_tracker_2d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  auto func = [=](element_t e) {
    feature_point_t p;
    if (check_simplex(e, p)) {
      std::lock_guard<std::mutex> guard(mutex);
      // std::cerr << e << std::endl;

      intersections[e] = p;

      auto tris = e.side_of(m);
      for (auto tri : tris) {
        if (tri.valid(m)) {
          auto tets = tri.side_of(m);
          for (auto tet : tets) 
            if (tet.valid(m))
              related_cells.insert(tet); // tets.begin(), tets.end());
        }
      }
    }
  };

  m.element_for_ordinal(1, current_timestep, func, xl, nthreads, enable_set_affinity);
  if (field_data_snapshots.size() >= 2) // interval
    m.element_for_interval(1, current_timestep, current_timestep+1, func, xl, nthreads, enable_set_affinity);
}

inline void contour_tracker_2d_regular::simplex_coordinates(
    const std::vector<std::vector<int>>& vertices, double X[][3]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j];
}

inline void contour_tracker_2d_regular::simplex_scalars(
    const std::vector<std::vector<int>>& vertices, double values[]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][2] == current_timestep ? 0 : 1;
    values[i] = field_data_snapshots[iv].scalar(
        vertices[i][0] - local_array_domain.start(0), 
        vertices[i][1] - local_array_domain.start(1));
  }
}

#if FTK_HAVE_VTK
inline void contour_tracker_2d_regular::write_isovolume_vtu(const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_isovolume_vtp();

    // write_polydata(filename, poly);
#if 1
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputData(poly);
    normalGenerator->ConsistencyOn();
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputePointNormalsOff();
    // normalGenerator->ComputeCellNormalsOn();
    // normalGenerator->SetFlipNormals(true);
    normalGenerator->AutoOrientNormalsOn();
    normalGenerator->Update();

    write_polydata(filename, normalGenerator->GetOutput());
#endif

#if 0
    vtkSmartPointer<vtkPCANormalEstimation> normals = vtkSmartPointer<vtkPCANormalEstimation>::New();
    normals->SetSampleSize(10);
    normals->SetInputData(poly);
    normals->FlipNormalsOn();
    normals->SetNormalOrientationToGraphTraversal();
    // normals->SetNormalOrientationToAsComputed();
    // normals->SetNormalOrientationToPoint();
    normals->Update();
    poly->GetPointData()->SetNormals( normals->GetOutput()->GetPointData()->GetArray(0) );
    write_polydata(filename, poly);
#endif
  }
}

vtkSmartPointer<vtkPolyData> contour_tracker_2d_regular::get_isovolume_vtp() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();

  auto add_triangle = [&cells](int i0, int i1, int i2) {
    vtkSmartPointer<vtkTriangle> tri = vtkTriangle::New();
    tri->GetPointIds()->SetId(0, i0);
    tri->GetPointIds()->SetId(1, i1);
    tri->GetPointIds()->SetId(2, i2);
    cells->InsertNextCell(tri);
  };

  auto my_intersections = intersections; // get_intersections();
  unsigned long long i = 0;
  for (auto &kv : my_intersections) {
    double p[3] = {kv.second.x[0], kv.second.x[1], kv.second.t};
    kv.second.tag = i ++;
    points->InsertNextPoint(p);
  }

  for (const auto &e : related_cells) {
  // m.element_for(3, [&](element_t e) {
    int count = 0;
    unsigned long long ids[4];

    std::set<element_t> unique_edges;
    for (auto tri : e.sides(m)) {
      for (auto edge : tri.sides(m)) {
        unique_edges.insert(edge);
      }
    }
    
    for (auto edge : unique_edges) 
      if (my_intersections.find(edge) != my_intersections.end())
        ids[count ++] = my_intersections[edge].tag;

    if (count == 3) {
      add_triangle(ids[0], ids[1], ids[2]);
    } else if (count == 4) { // quad
#if 0
      double p[4][3];
      for (int i = 0; i < 4; i ++)
        points->GetPoint(ids[i], p[i]);

      double A[3][3];
      for (int i = 0; i < 3; i ++)
        for (int j = 0; j < 3; j ++)
          A[i][j] = p[j][i];
      // print3x3("A", A);
      double b[3] = {p[3][0], p[3][1], p[3][2]};
      // fprintf(stderr, "b=%f, %f, %f\n", b[0], b[1], b[2]);

      double mu[3];
      solve_linear3x3(A, b, mu);
      // fprintf(stderr, "mu=%f, %f, %f, sum=%f\n", mu[0], mu[1], mu[2], 
      //     mu[0] + mu[1] + mu[2]);
#endif

      add_triangle(ids[0], ids[1], ids[2]);
      add_triangle(ids[1], ids[3], ids[2]);
    } else {
      // fprintf(stderr, "shoot..\n");
    };
  }; //, 1); 

  polyData->SetPoints(points);
  polyData->SetPolys(cells);

  return polyData;
}
#else
inline void contour_tracker_2d_regular::write_isovolume_vtu(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

}

#endif
