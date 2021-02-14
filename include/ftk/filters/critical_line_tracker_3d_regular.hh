#ifndef _FTK_CRITICAL_LINE_TRACKER_3D_REGULAR_HH
#define _FTK_CRITICAL_LINE_TRACKER_3D_REGULAR_HH

#include <ftk/config.hh>
#include <ftk/filters/critical_line_tracker.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>

namespace ftk {

struct critical_line_tracker_3d_regular : public virtual critical_line_tracker, public virtual regular_tracker
{
  critical_line_tracker_3d_regular(diy::mpi::communicator comm) : critical_line_tracker(comm), regular_tracker(comm, 3), tracker(comm) {}
  virtual ~critical_line_tracker_3d_regular() {}

  void finalize();
  void reset();

  void update_timestep();

public:
  void build_vortex_surfaces();

  void read_surfaces(const std::string& filename, std::string format="auto");

  void write_intersections(const std::string& filename) const;
  void write_sliced(const std::string& pattern) const;
  void write_surfaces(const std::string& filename, std::string format="auto") const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_intersections_vtp() const;
#endif

protected:
  typedef simplicial_regular_mesh_element element_t;
  
  std::map<element_t, feature_point_t> intersections;
  std::set<element_t> related_cells;

  feature_surface_t surfaces;

protected:
  virtual bool check_simplex(const element_t& s, feature_point_t& cp) const;
  
  void simplex_values(
      const std::vector<std::vector<int>>& vertices,
      float X[][4],
      float UV[][2]) const;

  virtual std::vector<std::string> varnames() const { return {}; } // varnames for additional variables stored in scalar
};


////////////////////
inline void critical_line_tracker_3d_regular::finalize()
{
  double max_accumulated_kernel_time;
  diy::mpi::reduce(comm, accumulated_kernel_time, max_accumulated_kernel_time, get_root_proc(), diy::mpi::maximum<double>());
  if (comm.rank() == get_root_proc())
    fprintf(stderr, "max_accumulated_kernel_time=%f\n", accumulated_kernel_time);
  
  diy::mpi::gather(comm, intersections, intersections, get_root_proc());
  diy::mpi::gather(comm, related_cells, related_cells, get_root_proc());
  
  if (comm.rank() == get_root_proc()) {
    fprintf(stderr, "#intersecttions=%zu, #related_cells=%zu\n", 
        intersections.size(), related_cells.size());
    build_vortex_surfaces();
  }
}

inline void critical_line_tracker_3d_regular::build_vortex_surfaces()
{
  fprintf(stderr, "building vortex surfaces...\n");

  int i = 0;
  for (auto &kv : intersections) {
    kv.second.id = i ++;
    surfaces.pts.push_back(kv.second);
  }

  std::mutex my_mutex;
  auto add_tri = [&](int i0, int i1, int i2) {
    std::lock_guard<std::mutex> guard(my_mutex);
    surfaces.tris.push_back({i0, i1, i2});
  };

  parallel_for<element_t>(related_cells, [&](const element_t &e) {
    int count = 0;
    int ids[6];

    std::set<element_t> unique_tris;
    for (auto tet : e.sides(m))
      for (auto tri : tet.sides(m))
        unique_tris.insert(tri);

    for (auto tri : unique_tris)
      if (intersections.find(tri) != intersections.end())
        ids[count ++] = intersections[tri].id;

    if (count == 3) {
      add_tri(ids[0], ids[1], ids[2]);
    } else if (count == 4) {
      // std::lock_guard<std::mutex> guard(my_mutex);
      // surfaces.quads.push_back({ids[0], ids[1], ids[2], ids[3]});
      add_tri(ids[0], ids[1], ids[2]);
      add_tri(ids[0], ids[1], ids[3]);
      add_tri(ids[0], ids[2], ids[3]);
      add_tri(ids[1], ids[2], ids[3]);
    } else if (count == 5) {
      // std::lock_guard<std::mutex> guard(my_mutex);
      // surfaces.pentagons.push_back({ids[0], ids[1], ids[2], ids[3], ids[4]});
      add_tri(ids[0], ids[1], ids[2]);
      add_tri(ids[0], ids[1], ids[3]);
      add_tri(ids[0], ids[1], ids[4]);
      add_tri(ids[0], ids[2], ids[3]);
      add_tri(ids[0], ids[2], ids[4]);
      add_tri(ids[0], ids[3], ids[4]);
      add_tri(ids[1], ids[2], ids[3]);
      add_tri(ids[1], ids[2], ids[4]);
      add_tri(ids[2], ids[3], ids[4]);
    } else {
      fprintf(stderr, "irregular count=%d\n", count); // WIP: triangulation
    }
  }, FTK_THREAD_PTHREAD, nthreads, enable_set_affinity);

  surfaces.relabel();
  fprintf(stderr, "#pts=%zu, #tri=%zu\n", surfaces.pts.size(), surfaces.tris.size());
}

inline void critical_line_tracker_3d_regular::reset()
{
  current_timestep = 0;

  field_data_snapshots.clear();
  intersections.clear();

  critical_line_tracker::reset();
}

inline bool critical_line_tracker_3d_regular::check_simplex(
    const element_t& e, feature_point_t& p) const
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex

  float X[3][4], UV[3][2], uv[3][2];
  simplex_values(vertices, X, UV);

  // locate zero
  float mu[3], // barycentric coordinates
        cond; // condition number

  const long long factor = 2 << 20; // TODO
  long long UVf[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++) {
      UVf[i][j] = factor * UV[i][j];
      uv[i][j] = UV[i][j];
    }
  
  int indices[3];
  simplex_indices(vertices, indices);
  
  bool succ = robust_critical_point_in_simplex2(UVf, indices);
  if (!succ) return false;

  inverse_lerp_s2v2(uv, mu, &cond);
  // mu[0] = mu[1] = mu[2] = 0.3333;
  clamp_barycentric<3, float>(mu);

  // interpolation
  float x[4];
  lerp_s2v4(X, mu, x);

  // result
  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];
  p.t = x[3];
  p.v[0] = mu[0]; // using v channel to store mu for now for the derived classes
  p.v[1] = mu[1];
  p.v[2] = mu[2];
  // p.scalar[0] = mu[0] * UV[0][2] + mu[1] * UV[1][2] + mu[2] * UV[2][2];
  p.cond = cond;
  p.tag = e.to_integer(m);
  p.ordinal = e.is_ordinal(m);
  p.timestep = current_timestep;

  // fprintf(stderr, "%f, %f, %f, %f\n", p[0], p[1], p[2], p[3]);

  return true;
}

inline void critical_line_tracker_3d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  // std::cerr << field_data_snapshots[0].uv.shape() << std::endl;

  typedef std::chrono::high_resolution_clock clock_type;
  auto t0 = clock_type::now();

  auto get_relatetd_cels = [&](element_t e) {
    std::set<element_t> my_related_cells;
    
    auto tets = e.side_of(m);
    for (auto tet : tets) {
      if (tet.valid(m)) {
        auto pents = tet.side_of(m);
        for (auto pent : pents)
          if (pent.valid(m))
            my_related_cells.insert(pent);
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
        intersections[e] = p;
        related_cells.insert(my_related_cells.begin(), my_related_cells.end());
      }
    }
  };

  element_for_ordinal(2, func);
  if (field_data_snapshots.size() >= 2) 
    element_for_interval(2, func);
  
  auto t1 = clock_type::now();
  accumulated_kernel_time += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9;
}
  
inline void critical_line_tracker_3d_regular::simplex_values(
      const std::vector<std::vector<int>>& vertices,
      float X[][4],
      float UV[][2]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    const auto &f = field_data_snapshots[iv];

    const auto idx = f.uv.index(std::vector<size_t>({
          0,
          vertices[i][0] - local_array_domain.start(0), 
          vertices[i][1] - local_array_domain.start(1), 
          vertices[i][2] - local_array_domain.start(2)}));
    UV[i][0] = f.uv[idx];
    UV[i][1] = f.uv[idx+1];
    // UV[i][2] = f.uv[idx+2];
    
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
  }
}

inline void critical_line_tracker_3d_regular::write_surfaces(const std::string& filename, std::string format) const 
{
  if (comm.rank() == get_root_proc()) {
    surfaces.save(filename, format);
  }
}

inline void critical_line_tracker_3d_regular::read_surfaces(const std::string& filename, std::string format)
{
  if (comm.rank() == get_root_proc()) {
    surfaces.load(filename, format);
    fprintf(stderr, "readed surfaces #pts=%zu, #tris=%zu\n", surfaces.pts.size(), surfaces.tris.size());
  }
}

inline void critical_line_tracker_3d_regular::write_intersections(const std::string& filename) const
{
#if FTK_HAVE_VTK
  if (comm.rank() == get_root_proc())
    write_polydata(filename, get_intersections_vtp());
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

inline void critical_line_tracker_3d_regular::write_sliced(const std::string& pattern) const
{
#if FTK_HAVE_VTK
  if (comm.rank() == get_root_proc()) {
    for (int i = 0; i < end_timestep; i ++) {
      auto sliced = surfaces.slice_time(i);
      fprintf(stderr, "sliced timestep %d, #curves=%zu\n", i, sliced.size());

      // auto poly = sliced.to_vtp(3, {"residue"}); // std::vector<std::string>());
      auto poly = sliced.to_vtp(3, varnames()); // std::vector<std::string>());
      
      const auto filename = series_filename(pattern, i);
      write_polydata(filename, poly);
    }
  }
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> critical_line_tracker_3d_regular::get_intersections_vtp() const
{
  vtkSmartPointer<vtkPolyData> poly = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
 
  vtkIdType pid[1];
  for (const auto &kv : intersections) {
    const auto &cp = kv.second;
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]}; 
    if (cpdims() == 2) p[2] = cp.t;
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  poly->SetPoints(points);
  poly->SetVerts(vertices);
  
  return poly;
}
#endif

}

#endif
