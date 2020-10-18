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
#include <ftk/ndarray/writer.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/filters/contour_tracker_regular.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy-ext/gather.hh>

#if FTK_HAVE_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkPointSet.h>
#include <vtkDelaunay3D.h>
#include <vtkCellIterator.h>
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
  void build_isovolume();

protected:
  void write_isovolume_vtu(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> get_isovolume_vtu() const; // legacy
#endif

  // void write_sliced_vtp(const std::string& pattern) const;
  void write_sliced_vtu(const std::string& pattern) const;
  feature_surface_t get_sliced(int t) const;

protected:
  typedef simplicial_regular_mesh_element element_t;

  feature_volume_t isovolume;
  
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

inline void contour_tracker_3d_regular::build_isovolume()
{
  fprintf(stderr, "building isovolumes...\n");

  int i = 0;
  for (auto &kv : intersections) {
    kv.second.id = i ++;
    isovolume.pts.push_back(kv.second);
  }

  std::mutex my_mutex;
  auto add_tet = [&](int i0, int i1, int i2, int i3) {
    std::lock_guard<std::mutex> guard(my_mutex);
    isovolume.conn.push_back({i0, i1, i2, i3});
  };
  
  parallel_for<element_t>(related_cells, nthreads, [&](const element_t &e) {
    int count = 0;
    int ids[6]; // 6 intersected edges max

    std::set<element_t> unique_edges;
    for (auto tet : e.sides(m))
      for (auto tri : tet.sides(m))
        for (auto edge : tri.sides(m))
          unique_edges.insert(edge);

    for (auto edge : unique_edges)
      if (intersections.find(edge) != intersections.end())
        ids[count ++] = intersections[edge].id;

    if (count == 4) {
      add_tet(ids[0], ids[1], ids[2], ids[3]);
    } else if (count == 6) {
      // triangulation. // if the pentachoron has six intersected 2-edges, the isovolume must be a prism
      // (1) find two tets, each of which has a triangular isosurface patch
      // (2) *sort vertex of each patch  *already sorted
      // (3) staircase triangulation

      int triangles[2][3], triangle_count = 0;

      for (auto tet : e.sides(m)) {
        int my_count = 0;
        int my_ids[4];

        std::set<element_t> my_unique_edges;
        for (auto tri : tet.sides(m))
          for (auto edge : tri.sides(m))
            my_unique_edges.insert(edge);

        for (auto edge : my_unique_edges)
          if (intersections.find(edge) != intersections.end())
            my_ids[my_count ++] = intersections[edge].id;

        // fprintf(stderr, "my_count=%d\n", my_count);
        if (my_count == 3) { // triangle
          for (int i = 0; i < 3; i ++)
            triangles[triangle_count][i] = my_ids[i];
          triangle_count ++;
        }
      }
      // fprintf(stderr, "triangle_count=%d\n", triangle_count);
      assert(triangle_count == 2);

      // staircase triangulation
      add_tet(triangles[0][0], triangles[0][1], triangles[0][2], triangles[1][2]);
      add_tet(triangles[0][0], triangles[0][1], triangles[1][1], triangles[1][2]);
      add_tet(triangles[0][0], triangles[1][0], triangles[1][1], triangles[1][2]);
    }
  });

  isovolume.relabel();
  fprintf(stderr, "isovolumes built, #pts=%zu, #tet=%zu\n", isovolume.pts.size(), isovolume.conn.size());
}

inline void contour_tracker_3d_regular::finalize()
{
  diy::mpi::gather(comm, intersections, intersections, get_root_proc());
  diy::mpi::gather(comm, related_cells, related_cells, get_root_proc());

  if (comm.rank() == get_root_proc()) {
    build_isovolume();
    // feature_volume_set_t isovolumes;
    
#if 0 // this is not quite efficient; better to directly build feature volume
    std::set<element_t> elements;
    for (const auto &kv : intersections)
      elements.insert(kv.first);

    auto cc = extract_connected_components<element_t, std::set<element_t>>(
        [&](element_t e) {
          std::set<element_t> neighbors;
          for (auto tri : e.side_of(m))
            for (auto tet : tri.side_of(m))
              for (auto pen : tet.side_of(m))
                for (auto tet : pen.sides(m))
                  for (auto tri : tet.sides(m))
                    for (auto edge : tri.sides(m))
                      neighbors.insert(edge);
          return neighbors;
        }, elements);

    fprintf(stderr, "#cc=%zu\n", cc.size());
#endif
  }
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
      std::set<element_t> my_related_cells;

      auto tris = e.side_of(m);
      for (auto tri : tris) {
        auto tets = tri.side_of(m);
        for (auto tet : tets) {
          if (tet.valid(m)) {
            auto pents = tet.side_of(m);
            for (auto pent : pents) 
              if (pent.valid(m))
                my_related_cells.insert(pent); 
          }
        }
      }
     
      {
        std::lock_guard<std::mutex> guard(mutex);
        // std::cerr << e << std::endl;
        intersections[e] = p;
        related_cells.insert(my_related_cells.begin(), my_related_cells.end());
      }
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

feature_surface_t contour_tracker_3d_regular::get_sliced(int t) const
{
  return isovolume.slice([&](const feature_point_t& p) {
    if (p.timestep == t && p.ordinal) return true;
    else return false;
  });
}

#if FTK_HAVE_VTK
inline void contour_tracker_3d_regular::write_isovolume_vtu(const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = 
      vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData( isovolume.to_vtu() );
    writer->Write();
  }
}

inline void contour_tracker_3d_regular::write_sliced_vtu(const std::string& pattern) const 
{
  if (comm.rank() == get_root_proc()) {
    for (int i = 0; i < current_timestep; i ++) {
      const auto filename = ndarray_writer<double>::filename(pattern, i);
      
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = 
        vtkXMLUnstructuredGridWriter::New();
      writer->SetFileName(filename.c_str());
      writer->SetInputData(get_sliced(i).to_vtu());
      writer->Write();
    }
  }
}
#else
inline void contour_tracker_3d_regular::write_isovolume_vtu(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}

inline void contour_tracker_3d_regular::write_sliced_vtu(const std::string& pattern) const 
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

}

#endif
