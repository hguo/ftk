#ifndef _FTK_CONTOUR_TRACKER_3D_REGULAR_HH
#define _FTK_CONTOUR_TRACKER_3D_REGULAR_HH

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
#include <ftk/ndarray/writer.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/filters/contour_tracker_regular.hh>
#include <ftk/utils/gather.hh>

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

extern std::vector<ftk::feature_point_lite_t> 
extract_contour_3dt_cuda(
    int scope, 
    int current_timestep, 
    const ftk::lattice& domain,
    const ftk::lattice& core, 
    const ftk::lattice& ext, 
    double threshold,
    const double *F_c, 
    const double *F_n);

namespace ftk {

struct contour_tracker_3d_regular : public contour_tracker_regular {
  contour_tracker_3d_regular(diy::mpi::communicator comm) : contour_tracker_regular(comm, 3), tracker(comm) {}
  virtual ~contour_tracker_3d_regular() {}

  // int cpdims() const { return 3; }

  void finalize();
  void reset();

  void update_timestep();

public:
  const feature_volume_t& get_isovolume() { return isovolume; }

protected:
  void build_isovolume();

protected:
  void write_isovolume_vtu(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> get_isovolume_vtu() const; // legacy
#endif

  // void write_sliced_vtp(const std::string& pattern) const;
  void write_sliced_vtu(const std::string& pattern) const;

protected:
  typedef simplicial_regular_mesh_element element_t;

  feature_volume_t isovolume;
  
protected:
  bool check_simplex(const element_t& s, feature_point_t& cp);
  void trace_intersections();
  void trace_connected_components();

  void simplex_coordinates(const std::vector<std::vector<int>>& vertices, double X[][4]) const;
  void simplex_values(const std::vector<std::vector<int>>& vertices, double scalars[], double gradients[][3]) const;
};


////////////////////
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
  
  parallel_for<element_t>(related_cells, [&](const element_t &e) {
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
      // (2) sort vertices of each patch
      // (3) staircase triangulation

      int quadrilaterals[3][4], quadrilaterals_count = 0;
      int triangles[2][3], triangle_count = 0;

      for (auto tet : e.sides(m)) {
        int my_count = 0;
        int my_ids[4];

        std::set<element_t> my_unique_edges;
        for (auto tri : tet.sides(m))
          if (tri.valid(m)) 
            for (auto edge : tri.sides(m))
              if (edge.valid(m))
                my_unique_edges.insert(edge);

        for (auto edge : my_unique_edges)
          if (intersections.find(edge) != intersections.end())
            my_ids[my_count ++] = intersections[edge].id;

        // fprintf(stderr, "my_count=%d\n", my_count);
        if (my_count == 3) { // triangle
          for (int i = 0; i < 3; i ++)
            triangles[triangle_count][i] = my_ids[i];
          triangle_count ++;
        } else if (my_count == 4) { // quadrilateral
          for (int i = 0; i < 4; i ++) 
            quadrilaterals[quadrilaterals_count][i] = my_ids[i];
          quadrilaterals_count ++;
        }
      }
      // fprintf(stderr, "triangle_count=%d\n", triangle_count);
      assert(triangle_count == 2);
      assert(quadrilaterals_count == 3);

      // staircase triangulatio8n
      add_tet(triangles[0][0], triangles[0][1], triangles[0][2], triangles[1][2]);
      add_tet(triangles[0][0], triangles[0][1], triangles[1][1], triangles[1][2]);
      add_tet(triangles[0][0], triangles[1][0], triangles[1][1], triangles[1][2]);
    }
  }, FTK_THREAD_PTHREAD, nthreads, enable_set_affinity);

  isovolume.relabel();
  fprintf(stderr, "isovolumes built, #pts=%zu, #tet=%zu\n", isovolume.pts.size(), isovolume.conn.size());
}

inline void contour_tracker_3d_regular::finalize()
{
  double max_accumulated_kernel_time;
  diy::mpi::reduce(comm, accumulated_kernel_time, max_accumulated_kernel_time, get_root_proc(), diy::mpi::maximum<double>());
  if (comm.rank() == get_root_proc())
    fprintf(stderr, "max_accumulated_kernel_time=%f\n", accumulated_kernel_time);
  
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
  isovolume.clear();
  contour_tracker_regular::reset();
}

inline bool contour_tracker_3d_regular::check_simplex(
    const element_t& e, feature_point_t& p)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex
  
  int indices[2];
  simplex_indices(vertices, indices);
  
  double f[2], g[2][3];
  simplex_values(vertices, f, g);
 
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

  double grad[3];
  lerp_s1v3(g, mu, grad);
  
  double ff = lerp_s1(f, mu) + threshold;

  double X[2][4], x[4];
  simplex_coordinates(vertices, X);
  lerp_s1v4(X, mu, x);

  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];
  p.t = x[3];
  p.v[0] = grad[0];
  p.v[1] = grad[1];
  p.v[2] = grad[2];
  p.scalar[0] = ff;
  p.tag = e.to_integer(m);
  p.ordinal = e.is_ordinal(m);
  p.timestep = current_timestep;

  // std::cerr << p.ordinal << ", " 
  //   << p.timestep << ", "
  //   << e << std::endl;
  // print_matrix<double, 2, 4>("X", X);
  // fprintf(stderr, "%f, %f, %f, %f\n", p.x[0], p.x[1], p.x[2], p.t);

  return true;
}

inline void contour_tracker_3d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  
  typedef std::chrono::high_resolution_clock clock_type;
  auto t0 = clock_type::now();
  
  auto get_relatetd_cels = [&](element_t e) {
    std::set<element_t> my_related_cells;
    
    auto tris = e.side_of(m);
    for (auto tri : tris) {
      if (tri.valid(m)) {
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
    }

    return my_related_cells;
  };

  auto func = [=](element_t e) {
    feature_point_t p;
    if (check_simplex(e, p)) {
      std::set<element_t> my_related_cells = get_relatetd_cels(e);
     
      {
        std::lock_guard<std::mutex> guard(mutex);
        // std::cerr << e << std::endl;
        intersections[e] = p;
        related_cells.insert(my_related_cells.begin(), my_related_cells.end());
      }
    }
  };

  if (xl == FTK_XL_CUDA) {
#if FTK_HAVE_CUDA
    ftk::lattice domain4({
          domain.start(0), 
          domain.start(1), 
          domain.start(2), 
          0
        }, {
          domain.size(0),
          domain.size(1),
          domain.size(2),
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
        {field_data_snapshots[0].scalar.dim(0), 
         field_data_snapshots[0].scalar.dim(1),
         field_data_snapshots[0].scalar.dim(2)});

    // ordinal
    auto results = extract_contour_3dt_cuda(
        ELEMENT_SCOPE_ORDINAL, 
        current_timestep, 
        domain4,
        ordinal_core,
        ext,
        threshold,
        field_data_snapshots[0].scalar.data(),
        NULL
      );
  
    for (auto lcp : results) {
      feature_point_t cp(lcp);
      element_t e(4, 1);
      e.from_work_index(m, cp.tag, ordinal_core, ELEMENT_SCOPE_ORDINAL);
      cp.tag = e.to_integer(m);
      cp.ordinal = true;
      cp.timestep = current_timestep;

      intersections[e] = cp;
      std::set<element_t> my_related_cells = get_relatetd_cels(e);
      related_cells.insert(my_related_cells.begin(), my_related_cells.end());
    }

    if (field_data_snapshots.size() >= 2) { // interval
      fprintf(stderr, "processing interval %d, %d\n", current_timestep, current_timestep + 1);
      auto results = extract_contour_3dt_cuda(
          ELEMENT_SCOPE_INTERVAL, 
          current_timestep,
          domain4,
          interval_core,
          ext,
          threshold,
          field_data_snapshots[0].scalar.data(),
          field_data_snapshots[1].scalar.data()
        );
      
      fprintf(stderr, "interval_results#=%zu\n", results.size());
      for (auto lcp : results) {
        feature_point_t cp(lcp);
        element_t e(4, 1);
        e.from_work_index(m, cp.tag, interval_core, ELEMENT_SCOPE_INTERVAL);
        cp.tag = e.to_integer(m);
        cp.ordinal = false;
        cp.timestep = current_timestep;
      
        intersections[e] = cp;
        std::set<element_t> my_related_cells = get_relatetd_cels(e);
        related_cells.insert(my_related_cells.begin(), my_related_cells.end());
      }
    }
#else
    fatal("FTK not compiled with CUDA.");
#endif

  } else {
    element_for_ordinal(1, func);
    if (field_data_snapshots.size() >= 2) // interval
      element_for_interval(1, func);
  }
  
  auto t1 = clock_type::now();
  accumulated_kernel_time += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9;
}

inline void contour_tracker_3d_regular::simplex_coordinates(
    const std::vector<std::vector<int>>& vertices, double X[][4]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
}

inline void contour_tracker_3d_regular::simplex_values(
    const std::vector<std::vector<int>>& vertices, double scalars[], double grads[][3]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    scalars[i] = field_data_snapshots[iv].scalar(
        vertices[i][0] - local_array_domain.start(0), 
        vertices[i][1] - local_array_domain.start(1), 
        vertices[i][2] - local_array_domain.start(2));
    for (int j = 0; j < 3; j ++)
      grads[i][j] = field_data_snapshots[iv].gradient(
          j, vertices[i][0], vertices[i][1], vertices[i][2]);
  }
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
      const auto filename = series_filename(pattern, i);
      
      auto poly = isovolume.slice_time(i).to_vtp();
   
      vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
      normalGenerator->SetInputData(poly);
      normalGenerator->ConsistencyOn();
      normalGenerator->ComputePointNormalsOff();
      normalGenerator->ComputeCellNormalsOn();
      // normalGenerator->SetFlipNormals(true);
      // normalGenerator->AutoOrientNormalsOn();
      normalGenerator->Update();

      write_polydata(filename, normalGenerator->GetOutput());
    
#if 0
      vtkSmartPointer<vtkPLYWriter> writer = vtkPLYWriter::New();
      writer->SetFileName(filename.c_str());
      auto sliced = isovolume.slice_time(i);
      sliced.reorient();
      writer->SetInputData(sliced.to_vtp());
      writer->Write();
#endif
#if 0
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = 
        vtkXMLUnstructuredGridWriter::New();
      writer->SetFileName(filename.c_str());

      auto sliced = isovolume.slice_time(i);
      sliced.reorient();

      writer->SetInputData(sliced.to_vtu());
      writer->Write();
#endif
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
