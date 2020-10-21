#ifndef _FTK_TDGL_TRACKER_3D_REGULAR_HH
#define _FTK_TDGL_TRACKER_3D_REGULAR_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/tdgl_vortex_tracker.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/geometry/points2vtk.hh>

#if FTK_HAVE_VTK
#include <vtkPolyDataNormals.h>
#include <vtkPLYWriter.h>
#endif

namespace ftk {

struct tdgl_vortex_tracker_3d_regular : public tdgl_vortex_tracker
{
  tdgl_vortex_tracker_3d_regular() : m(4) {}
  virtual ~tdgl_vortex_tracker_3d_regular() {}
  
  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void initialize();
  void finalize();
  void reset();

  void update_timestep();

public:
  void build_vortex_surfaces();

  void write_intersections_vtp(const std::string& filename) const;
  void write_sliced_vtp(const std::string& pattern) const {}
  void write_surfaces_vtp(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_intersections_vtp() const;
#endif

protected:
  simplicial_regular_mesh m;
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

  template <typename T> inline static T mod2pi(T x) { T y = fmod(x, 2*M_PI); if (y<0) y+= 2*M_PI; return y; }
  template <typename T> static T mod2pi1(T x) { return mod2pi(x + M_PI) - M_PI; }

protected: // config
  lattice domain, array_domain;
};


////////////////////
inline void tdgl_vortex_tracker_3d_regular::initialize()
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

inline void tdgl_vortex_tracker_3d_regular::finalize()
{
  diy::mpi::gather(comm, intersections, intersections, get_root_proc());
  diy::mpi::gather(comm, related_cells, related_cells, get_root_proc());
  
  if (comm.rank() == get_root_proc()) {
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

  parallel_for<element_t>(related_cells, nthreads, [&](const element_t &e) {
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
      std::lock_guard<std::mutex> guard(my_mutex);
      surfaces.quads.push_back({ids[0], ids[1], ids[2], ids[3]});
      // add_tri(ids[0], ids[1], ids[2]);
      // add_tri(ids[1], ids[3], ids[2]);
    }
    // fprintf(stderr, "count=%d\n", count); // WIP: triangulation
  });

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
  p.cond = cond;
  p.tag = e.to_integer(m);
  p.timestep = current_timestep;

  // fprintf(stderr, "%f, %f, %f, %f\n", p[0], p[1], p[2], p[3]);

  return true;
}

inline void tdgl_vortex_tracker_3d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  auto func = [=](element_t e) {
    feature_point_t p;
    if (check_simplex(e, p)) {
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

      {
        std::lock_guard<std::mutex> guard(mutex);
        intersections[e] = p;
        related_cells.insert(my_related_cells.begin(), my_related_cells.end());
      }
    }
  };

  m.element_for_ordinal(2, current_timestep, func, nthreads);
  if (field_data_snapshots.size() >= 2) // interval
    m.element_for_interval(2, current_timestep, current_timestep+1, func, nthreads);
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

    const size_t idx = f.rho.index(std::vector<int>({vertices[i][0], vertices[i][1], vertices[i][2]}));
    rho[i] = f.rho[idx];
    phi[i] = f.phi[idx];
    re[i] = f.re[idx];
    im[i] = f.im[idx];
    // rho[i] = f.rho(vertices[i][0], vertices[i][1], vertices[i][2]);
    // phi[i] = f.phi(vertices[i][0], vertices[i][1], vertices[i][2]);
    // re[i] = f.re(vertices[i][0], vertices[i][1], vertices[i][2]);
    // im[i] = f.im(vertices[i][0], vertices[i][1], vertices[i][2]);
    
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j] * f.meta.cell_lengths[j] + f.meta.origins[j];
      
    for (int j = 0; j < 4; j ++)
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
  
#if FTK_HAVE_VTK
inline void tdgl_vortex_tracker_3d_regular::write_intersections_vtp(const std::string& filename) const
{
  if (comm.rank() == get_root_proc())
    write_vtp(filename, get_intersections_vtp());
}

inline void tdgl_vortex_tracker_3d_regular::write_surfaces_vtp(const std::string& filename) const 
{
  if (comm.rank() == get_root_proc()) {
    auto poly = surfaces.to_vtp();
      
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputData(poly);
    normalGenerator->ConsistencyOn();
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOn();
    // normalGenerator->SetFlipNormals(true);
    normalGenerator->AutoOrientNormalsOn();
    normalGenerator->Update();

    write_vtp(filename, poly);
#if 0
    vtkSmartPointer<vtkPLYWriter> writer = vtkPLYWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
#endif
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
    if (cpdims() == 2) p[2] = cp.t;
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  poly->SetPoints(points);
  poly->SetVerts(vertices);
  
  return poly;
}
#else
inline void tdgl_vortex_tracker_3d_regular::write_intersections_vtp(const std::string& filename) const
{
  fatal("FTK not compiled with VTK.");
}
#endif

}

#endif
