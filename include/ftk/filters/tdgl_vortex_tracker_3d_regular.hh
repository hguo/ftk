#ifndef _FTK_TDGL_TRACKER_3D_REGULAR_HH
#define _FTK_TDGL_TRACKER_3D_REGULAR_HH

#include <ftk/config.hh>
#include <ftk/filters/tdgl_vortex_tracker.hh>
#include <ftk/filters/regular_tracker.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/write_polydata.hh>
#include <ftk/ndarray/writer.hh>
#include <ftk/numeric/fmod.hh>

extern std::vector<ftk::feature_point_lite_t> 
extract_tdgl_vortex_3dt_cuda(
    int scope, 
    int current_timestep, 
    const ftk::lattice& domain,
    const ftk::lattice& core, 
    const ftk::lattice& ext, 
    const ftk::tdgl_metadata_t &h_c,
    const ftk::tdgl_metadata_t &h_n,
    const float *rho_c, 
    const float *rho_l,
    const float *phi_c, 
    const float *phi_l);

namespace ftk {

struct tdgl_vortex_tracker_3d_regular : public virtual tdgl_vortex_tracker, public virtual regular_tracker
{
  tdgl_vortex_tracker_3d_regular(diy::mpi::communicator comm) : tdgl_vortex_tracker(comm), regular_tracker(comm, 3), tracker(comm) {}
  virtual ~tdgl_vortex_tracker_3d_regular() {}

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
  bool check_simplex(const element_t& s, feature_point_t& cp);
  
  void simplex_values(
      const std::vector<std::vector<int>>& vertices,
      float X[][4],
      float A[][3], 
      float rho[], float phi[],
      float re[], float im[]);

  void magnetic_potential(const tdgl_metadata_t& m, const float X[3], float A[3]) const;

  static float line_integral(float X0[], float X1[], float A0[], float A1[]);

  // template <typename T> inline static T mod2pi(T x) { T y = fmod(x, 2*M_PI); if (y<0) y+= 2*M_PI; return y; }
  // template <typename T> static T mod2pi1(T x) { return mod2pi(x + M_PI) - M_PI; }
};


////////////////////
inline void tdgl_vortex_tracker_3d_regular::finalize()
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

inline void tdgl_vortex_tracker_3d_regular::build_vortex_surfaces()
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

inline void tdgl_vortex_tracker_3d_regular::reset()
{
  current_timestep = 0;

  field_data_snapshots.clear();
  intersections.clear();

  tdgl_vortex_tracker::reset();
}

inline bool tdgl_vortex_tracker_3d_regular::check_simplex(
    const element_t& e, feature_point_t& p)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex

  float X[3][4], // coordinates
        A[3][3]; // averaged magnetic potential
  float rho[3], phi[3], re[3], im[3]; // values
  simplex_values(vertices, X, A, rho, phi, re, im);
  // fprintf(stderr, "rho=%f, %f, %f, phi=%f, %f, %f\n", rho[0], rho[1], rho[2], phi[0], phi[1], phi[2]);

  // compute contour integral
  float delta[3], phase_shift = 0;
  for (int i = 0; i < 3; i ++) { // ignoring quasi periodical boundary conditions
    int j = (i+1) % 3;
    float li = line_integral(X[i], X[j], A[i], A[j]);
    delta[i] = mod2pi1( phi[j] - phi[i] - li ); // gauge transformation
    phase_shift -= delta[i];
  }

  // check contour integral
  float critera = phase_shift / (2 * M_PI);
  if (std::abs(critera) < 0.5) return false; // ignoring chiralities

  // fprintf(stderr, "delta=%f, %f, %f\n", delta[0], delta[1], delta[2]);

  // guage transformation
  float psi[3][2]; // in re/im
  for (int i = 0; i < 3; i ++) {
    if (i != 0) phi[i] = phi[i-1] + delta[i-1];
    psi[i][0] = rho[i] * cos(phi[i]);
    psi[i][1] = rho[i] * sin(phi[i]);
  }

  // locate zero
  float mu[3], // barycentric coordinates
        cond; // condition number
  inverse_lerp_s2v2(psi, mu, &cond);
  // mu[0] = mu[1] = mu[2] = 0.3333;

  // interpolation
  float x[4];
  lerp_s2v4(X, mu, x);

  // result
  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];
  p.t = x[3];
  // p.cond = cond;
  p.tag = e.to_integer(m);
  p.ordinal = e.is_ordinal(m);
  p.timestep = current_timestep;

  // fprintf(stderr, "%f, %f, %f, %f\n", p[0], p[1], p[2], p[3]);

  return true;
}

inline void tdgl_vortex_tracker_3d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  
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
        {field_data_snapshots[0].rho.dim(0), 
         field_data_snapshots[0].rho.dim(1),
         field_data_snapshots[0].rho.dim(2)});

    // ordinal
    auto results = extract_tdgl_vortex_3dt_cuda(
        ELEMENT_SCOPE_ORDINAL, 
        current_timestep, 
        domain4,
        ordinal_core,
        ext,
        field_data_snapshots[0].meta,
        field_data_snapshots[0].meta,
        field_data_snapshots[0].rho.data(),
        NULL, // gradV[0].data(),
        field_data_snapshots[0].phi.data(),
        NULL // scalar[0].data(),
      );
  
    for (auto lcp : results) {
      feature_point_t cp(lcp);
      element_t e(4, 2);
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
      auto results = extract_tdgl_vortex_3dt_cuda(
          ELEMENT_SCOPE_INTERVAL, 
          current_timestep,
          domain4,
          interval_core,
          ext,
          field_data_snapshots[0].meta,
          field_data_snapshots[1].meta,
          field_data_snapshots[0].rho.data(),
          field_data_snapshots[1].rho.data(),
          field_data_snapshots[0].phi.data(),
          field_data_snapshots[1].phi.data()
        );
      
      fprintf(stderr, "interval_results#=%zu\n", results.size());
      for (auto lcp : results) {
        feature_point_t cp(lcp);
        element_t e(4, 2);
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
    element_for_ordinal(2, func);
    if (field_data_snapshots.size() >= 2) 
      element_for_interval(2, func);
  }
  
  auto t1 = clock_type::now();
  accumulated_kernel_time += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9;
}
  
inline void tdgl_vortex_tracker_3d_regular::simplex_values(
      const std::vector<std::vector<int>>& vertices,
      float X[][4],
      float A[][3],
      float rho[], float phi[],
      float re[], float im[])
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    const auto &f = field_data_snapshots[iv];

    const auto idx = f.rho.index(std::vector<size_t>({
          vertices[i][0] - local_array_domain.start(0), 
          vertices[i][1] - local_array_domain.start(1), 
          vertices[i][2] - local_array_domain.start(2)}));
    rho[i] = f.rho[idx];
    phi[i] = f.phi[idx];
    re[i] = f.re[idx];
    im[i] = f.im[idx];
    
    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j] * f.meta.cell_lengths[j] + f.meta.origins[j];
    X[i][3] = vertices[i][3];
      
    magnetic_potential(f.meta, X[i], A[i]);
  }
}
  
inline float tdgl_vortex_tracker_3d_regular::line_integral(float X0[], float X1[], float A0[], float A1[])
{
  float dX[3] = {X1[0] - X0[0], X1[1] - X0[1], X1[2] - X0[2]};
  float A[3]  = {A0[0] + A1[0], A0[1] + A1[1], A0[2] + A1[2]};
  return 0.5 * inner_product3(A, dX);
}

inline void tdgl_vortex_tracker_3d_regular::magnetic_potential(const tdgl_metadata_t& m, const float X[3], float A[3]) const
{
  if (m.B[1] > 0) {
    A[0] = -m.Kex;
    A[1] = X[0] * m.B[2];
    A[2] = -X[0] * m.B[1];
  } else {
    A[0] = -X[1] * m.B[2] - m.Kex;
    A[1] = 0;
    A[2] = X[1] * m.B[0];
  }
}

inline void tdgl_vortex_tracker_3d_regular::write_surfaces(const std::string& filename, std::string format) const 
{
  if (comm.rank() == get_root_proc()) {
    surfaces.save(filename, format);
  }
}

inline void tdgl_vortex_tracker_3d_regular::read_surfaces(const std::string& filename, std::string format)
{
  if (comm.rank() == get_root_proc()) {
    surfaces.load(filename, format);
    fprintf(stderr, "readed surfaces #pts=%zu, #tris=%zu\n", surfaces.pts.size(), surfaces.tris.size());
  }
}

  
#if FTK_HAVE_VTK
inline void tdgl_vortex_tracker_3d_regular::write_intersections(const std::string& filename) const
{
  if (comm.rank() == get_root_proc())
    write_polydata(filename, get_intersections_vtp());
}

inline void tdgl_vortex_tracker_3d_regular::write_sliced(const std::string& pattern) const
{
  if (comm.rank() == get_root_proc()) {
    for (int i = 0; i < end_timestep; i ++) {
      auto sliced = surfaces.slice_time(i);
      fprintf(stderr, "sliced timestep %d, #curves=%zu\n", i, sliced.size());

      auto poly = sliced.to_vtp(std::vector<std::string>());
      
      const auto filename = series_filename(pattern, i);
      write_polydata(filename, poly);
    }
  }
}

inline vtkSmartPointer<vtkPolyData> tdgl_vortex_tracker_3d_regular::get_intersections_vtp() const
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
#else
inline void tdgl_vortex_tracker_3d_regular::write_intersections(const std::string& filename) const
{
  fatal("FTK not compiled with VTK.");
}

inline void tdgl_vortex_tracker_3d_regular::write_sliced(const std::string& pattern) const
{
  fatal("FTK not compiled with VTK.");
}
#endif

}

#endif
