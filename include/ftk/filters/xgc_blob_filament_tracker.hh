#ifndef _FTK_XGC_BLOB_FILAMENT_TRACKER_HH
#define _FTK_XGC_BLOB_FILAMENT_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/xgc_tracker.hh>

#if FTK_HAVE_TBB
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_unordered_set.h>
#endif

#if FTK_HAVE_CUDA
#include "xgc_filament_tracker.cuh"
#endif

namespace ftk {
  
struct xgc_blob_filament_tracker : public xgc_tracker {
  xgc_blob_filament_tracker(diy::mpi::communicator comm,
      std::shared_ptr<simplicial_xgc_3d_mesh<>> mx);
  virtual ~xgc_blob_filament_tracker();

  // xgc_blob_filament_tracker(diy::mpi::communicator comm, 
  //     std::shared_ptr<simplicial_unstructured_2d_mesh<>> m2, 
  //     int nphi_, int iphi_);

  int cpdims() const { return 3; }
 
  void initialize();
  void reset() {
    field_data_snapshots.clear();
  }
  void update() {}
  void finalize();

public:
  void update_timestep();
  
  void push_field_data_snapshot(
      const ndarray<double> &scalar, 
      const ndarray<double> &vector,
      const ndarray<double> &jacobian);
  void push_field_data_snapshot(
      const ndarray<double> &scalar);

protected:
  bool check_simplex(int, feature_point_t& cp);
  void check_penta(int);
  void add_penta_tri(int, int, int);
 
  template <int n, typename T>
  void simplex_values(
      const int i,
      int verts[n],
      int t[n],
      int p[n],
      T rzpt[n][4], 
      T f[n], // scalars
      T v[n][2], // vectors
      T j[n][2][2]); // jacobians

public:
  void build_critical_surfaces();
  void post_process_surfaces();
  void post_process_curves(feature_curve_set_t& curves) const;

public:
  void write_intersections_binary(const std::string& filename) const;
  void write_intersections(const std::string& filename, std::string format="") const;
  void read_intersections_binary(const std::string& filename);
  
  void write_surfaces(const std::string& filename, std::string format="", bool torus=true) const;
  void read_surfaces(const std::string& filename, std::string format="");

  void write_sliced(const std::string& filename, std::string format="", bool torus=true) const;

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_intersections_vtp(bool torus = true) const;
  vtkSmartPointer<vtkPolyData> get_critical_line_vtp(bool torus = true) const;
  vtkSmartPointer<vtkPolyData> get_critical_surfaces_vtp(bool torus = true) const;

  vtkSmartPointer<vtkPolyData> transform_vtp_coordinates(vtkSmartPointer<vtkPolyData>) const;
#endif

protected:
#if FTK_HAVE_TBB
  tbb::concurrent_hash_map<int, feature_point_t> intersections;
  tbb::concurrent_unordered_set<int> related_cells;
#else
  std::map<int, feature_point_t> intersections;
  std::set<int> related_cells;
#endif
  feature_surface_t surfaces;

protected:
#if FTK_HAVE_CUDA
  xft_ctx_t *ctx = NULL;
#endif
};

/////

xgc_blob_filament_tracker::xgc_blob_filament_tracker(
    diy::mpi::communicator comm, 
    std::shared_ptr<simplicial_xgc_3d_mesh<>> mx) :
  xgc_tracker(comm, mx)
{
}

xgc_blob_filament_tracker::~xgc_blob_filament_tracker()
{
#if FTK_HAVE_CUDA
  if (xl == FTK_XL_CUDA)
    xft_destroy_ctx(&ctx);
#endif
}

inline void xgc_blob_filament_tracker::initialize()
{
#if FTK_HAVE_CUDA
  if (xl == FTK_XL_CUDA) {
    xft_create_ctx(&ctx);
    xft_load_mesh(ctx, 
        m3->get_nphi(), m3->get_iphi(), m3->get_vphi(), 
        m2->n(0), m2->n(1), m2->n(2), 
        m2->get_coords().data(), 
        m2->get_edges().data(), 
        m2->get_triangles().data());
    xft_load_interpolants(ctx, 
        m3->get_interpolants());
    xft_load_smoothing_kernel(ctx, 
        m2->get_smoothing_kernel_size(),
        m2->get_smoothing_kernel());
  }
#endif
}

inline void xgc_blob_filament_tracker::push_field_data_snapshot(
      const ndarray<double> &scalar, 
      const ndarray<double> &vector,
      const ndarray<double> &jacobian)
{
#if FTK_HAVE_CUDA
  if (xl == FTK_XL_CUDA) 
    xft_load_data(ctx, scalar.data(), vector.data(), jacobian.data());
#endif
  xgc_tracker::push_field_data_snapshot(scalar, vector, jacobian);
}

inline void xgc_blob_filament_tracker::push_field_data_snapshot(
      const ndarray<double> &scalar)
{
#if FTK_HAVE_CUDA
  if (xl == FTK_XL_CUDA) {
    xft_load_scalar_data(ctx, scalar.data());
    xgc_tracker::push_field_data_snapshot(ndarray<double>(), ndarray<double>(), ndarray<double>()); // empty
  } else 
    xgc_tracker::push_field_data_snapshot(scalar);
#else
  xgc_tracker::push_field_data_snapshot(scalar);
#endif
}

inline void xgc_blob_filament_tracker::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
 
  auto get_related_cels = [&](int i) {
    std::set<int> my_related_cells;
    
    auto tets = m4->side_of(2, i);
    for (auto tet : tets) {
      if (1) { // TODO: if valid tet
        auto pents = m4->side_of(3, tet);
        for (auto pent : pents)
          if (1) // TODO if valid pent
            my_related_cells.insert(pent);
      }
    }

    return my_related_cells;
  };

  auto func = [&](int i) {
    feature_point_t cp;
    if (check_simplex(i, cp)) {
      std::set<int> my_related_cells = get_related_cels(i);

      {
#if !FTK_HAVE_TBB
        std::lock_guard<std::mutex> guard(mutex);
#endif

        intersections.insert({i, cp});
        related_cells.insert(my_related_cells.begin(), my_related_cells.end());
      }
    }
  };

  if (xl == FTK_XL_CUDA) {
#if FTK_HAVE_CUDA
    // ordinal
    xft_execute(ctx, 1 /* ordinal */, current_timestep);
    std::vector<feature_point_lite_t> results(ctx->hcps, ctx->hcps + ctx->hncps);
    for (auto lcp : results) {
      feature_point_t cp(lcp);
      cp.tag += current_timestep * m4->n(2);
      cp.ordinal = true;
      cp.timestep = current_timestep;

      intersections.insert({cp.tag, cp});

      std::set<int> related = get_related_cels(cp.tag);
      related_cells.insert( related.begin(), related.end() );
    }

    // interval
    if (field_data_snapshots.size() >= 2) {
      xft_execute(ctx, 2 /* interval */, current_timestep);
      // fprintf(stderr, "** current_timestep=%d, gpu done interval.\n", current_timestep);
      std::vector<feature_point_lite_t> results(ctx->hcps, ctx->hcps + ctx->hncps);
      for (auto lcp : results) {
        feature_point_t cp(lcp);
        cp.tag += current_timestep * m4->n(2);
        cp.ordinal = false;
        cp.timestep = current_timestep; // - 1;

        intersections.insert({cp.tag, cp});

        std::set<int> related = get_related_cels(cp.tag);
        related_cells.insert( related.begin(), related.end() );
      }
    }
#else
    fatal("FTK not compiled with CUDA.");
#endif
  } else {
    // TODO: strange that m4->element_for has proformance problem...
    // m4->element_for_ordinal(2, current_timestep, func, xl, nthreads, false); // enable_set_affinity);
    const auto nd = m4->n(2), no = m4->n_ordinal(2), ni = m4->n_interval(2);
    parallel_for(no, [&](int i) {func(i + current_timestep * nd);}, xl, nthreads, enable_set_affinity);

    if (field_data_snapshots.size() >= 2)
      parallel_for(ni, [&](int i) {func(i + no + current_timestep * nd);}, xl, nthreads, enable_set_affinity);
      //   m4->element_for_interval(2, current_timestep, func, xl, nthreads, false); // enable_set_affinity);
  }
}

inline bool xgc_blob_filament_tracker::check_simplex(int i, feature_point_t& cp)
{
  int tri[3], t[3], p[3];
  double rzpt[3][4], f[3], v[3][2], j[3][2][2];
  long long vf[3][2];
  const long long factor = 1 << 15; // WIP
 
  simplex_values<3, double>(i, tri, t, p, rzpt, f, v, j);

  // print3x2("rz", rz);
  // print3x2("v", v);

  for (int k = 0; k < 3; k ++) 
    for (int l = 0; l < 2; l ++) 
      vf[k][l] = factor * v[k][l];

  bool succ = ftk::robust_critical_point_in_simplex2(vf, tri);
  if (!succ) return false;

  double mu[3], x[4];
  bool succ2 = ftk::inverse_lerp_s2v2(v, mu);
  ftk::clamp_barycentric<3>(mu);

  ftk::lerp_s2v4(rzpt, mu, x);
  for (int k = 0; k < 3; k ++)
    cp.x[k] = x[k];
  cp.t = x[3];

  cp.scalar[0] = f[0] * mu[0] + f[1] * mu[1] + f[2] * mu[2];

  double h[2][2];
  ftk::lerp_s2m2x2(j, mu, h);
  cp.type = ftk::critical_point_type_2d(h, true);

  cp.tag = i;
  cp.ordinal = m4->is_ordinal(2, i); 
  cp.timestep = current_timestep;

  // fprintf(stderr, "succ, mu=%f, %f, %f, x=%f, %f, %f, %f, timestep=%d, type=%d, t=%d, %d, %d\n", 
  //     mu[0], mu[1], mu[2], 
  //     x[0], x[1], x[2], x[3], cp.timestep, cp.type, 
  //     t[0], t[1], t[2]);
  
  return true;
}

void xgc_blob_filament_tracker::add_penta_tri(int i0, int i1, int i2)
{
  std::lock_guard<std::mutex> guard(mutex);
  // tri_count ++;
  // fprintf(stderr, "pushing %d, %d, %d, count=%zu\n", i0, i1, i2, surfaces.tris.size()); // , tri_count);
  surfaces.tris.push_back({i0, i1, i2});
}

void xgc_blob_filament_tracker::check_penta(int e)
{
  int penta[5], t[5];
  m4->get_simplex(4, e, penta);
  for (int i = 0; i < 5; i ++) {
    t[i] = m4->flat_vertex_time(penta[i]);
    if (t[i] < 0 || t[i] >= end_timestep)
      return;
  }

  int count = 0;
  int ids[20]; // some large buffer
  
  std::set<int> unique_tris;
  for (auto tet : m4->sides(4, e))
    for (auto tri : m4->sides(3, tet))
      unique_tris.insert(tri);
  
  // sanity check
#if 0
  if (unique_tris.size() != 10) { 
    fprintf(stderr, "penta %d, penta_type=%d\n", e, m4->simplex_type(4, e));
    for (auto tet : m4->sides(4, e)) {
      fprintf(stderr, "--tet %d, tet_type=%d\n", tet, m4->simplex_type(3, tet));
      for (auto tri : m4->sides(3, tet))
        fprintf(stderr, "----tri %d, tri_type=%d\n", tri, m4->simplex_type(2, tri));
    }
    fprintf(stderr, "#unique_tris=%zu\n", unique_tris.size());
    for (auto tri : unique_tris) 
      fprintf(stderr, "--tri %d\n", tri);
  }
#endif
  assert( unique_tris.size() == 10 );

  for (auto tri : unique_tris) {
#if FTK_HAVE_TBB
    tbb::concurrent_hash_map<int, feature_point_t>::const_accessor it;
    bool found = intersections.find(it, tri);
    if (found)
      ids[count ++] = it->second.id;
#else
    if (intersections.find(tri) != intersections.end())
      ids[count ++] = intersections[tri].id;
#endif
  }

  if (count == 0) return;
  else if (count == 3) {
    // std::lock_guard<std::mutex> guard(my_mutex);
    add_penta_tri(ids[0], ids[1], ids[2]);
  } else if (count == 4) {
    // std::lock_guard<std::mutex> guard(my_mutex);
    // surfaces.quads.push_back({ids[0], ids[1], ids[2], ids[3]});
    add_penta_tri(ids[0], ids[1], ids[2]);
    add_penta_tri(ids[0], ids[1], ids[3]);
    add_penta_tri(ids[0], ids[2], ids[3]);
    add_penta_tri(ids[1], ids[2], ids[3]);
  } else if (count == 5) {
    // std::lock_guard<std::mutex> guard(my_mutex);
    // surfaces.pentagons.push_back({ids[0], ids[1], ids[2], ids[3], ids[4]});
    add_penta_tri(ids[0], ids[1], ids[2]);
    add_penta_tri(ids[0], ids[1], ids[3]);
    add_penta_tri(ids[0], ids[1], ids[4]);
    add_penta_tri(ids[0], ids[2], ids[3]);
    add_penta_tri(ids[0], ids[2], ids[4]);
    add_penta_tri(ids[0], ids[3], ids[4]);
    add_penta_tri(ids[1], ids[2], ids[3]);
    add_penta_tri(ids[1], ids[2], ids[4]);
    add_penta_tri(ids[2], ids[3], ids[4]);
  } else {
    fprintf(stderr, "irregular count=%d, penta=%d: %d, %d, %d, %d, %d, t=%d, %d, %d, %d, %d\n", 
        count, e, 
        penta[0], penta[1], penta[2], penta[3], penta[4], 
        t[0], t[1], t[2], t[3], t[4]); // WIP: triangulation
  }
}

void xgc_blob_filament_tracker::finalize()
{
  diy::mpi::gather(comm, intersections, intersections, get_root_proc()); // TODO
  diy::mpi::gather(comm, related_cells, related_cells, get_root_proc());
  
  if (is_root_proc()) {
    fprintf(stderr, "#intersections=%zu, #related_cells=%zu\n", 
        intersections.size(), related_cells.size());
    build_critical_surfaces();
  }
}

template<int n, typename T>
void xgc_blob_filament_tracker::simplex_values(
    const int i, 
    int verts[n],
    int t[n], // time
    int p[n], // poloidal
    T rzpt[n][4], 
    T f[n],
    T v[n][2],
    T j[n][2][2])
{
  const int m2n0 = m2->n(0), m3n0 = m3->n(0);
  const int nphi = m3->get_nphi(), 
            iphi = m3->get_iphi(),
            vphi = m3->get_vphi();
  const int n0ivphi = m2n0 * vphi * nphi;

  m4->get_simplex(n-1, i, verts);
  for (int k = 0; k < n; k ++) {
    t[k] = verts[k] / m4->n(0); // m4->flat_vertex_time(verts[k]);
    const int v3 = verts[k] % m3->n(0); // vertex id in m3 mesh
    const int v2 = v3 % m2->n(0); // vertex id in m2 mesh
    p[k] = v3 / m2->n(0); // poloidal plane coords

    // if (m3->is_poloidal(p[k])) fprintf(stderr, "non-poloidal virtual plane\n");
    m3->get_coords_rzp(v3, rzpt[k]);
    rzpt[k][3] = t[k];

    const int iv = (t[k] == current_timestep) ? 0 : 1; 
    const field_data_snapshot_t &data = field_data_snapshots[iv];

#if 1
    m3->interpolate(data.scalar, data.vector, data.jacobian, v3, &f[k], v[k], j[k]);
#else
    const int dp0 = p[k] / iphi / vphi, // plane id in the density array
              dp1 = dp0 + 1; // (dp0 + 1) % nphi;
    const int vc0 = v2 + dp0 * m2n0;
   
    // if (m3->is_poloidal(p[k])) {
    // if (p[k] % vphi == 0) { // poloidal
    if (1) { 
      f[k] = data.scalar[vc0];

      v[k][0] = data.vector[vc0*2];
      v[k][1] = data.vector[vc0*2+1];

      j[k][0][0] = data.jacobian[vc0*4];
      j[k][0][1] = data.jacobian[vc0*4+1];
      j[k][1][0] = data.jacobian[vc0*4+2];
      j[k][1][1] = data.jacobian[vc0*4+3];
    } else { // virtual polodal plane, using linear interpolation
      int next = m2->nextnode(v2);
      const int vc1 = next + (dp1 % nphi) * m2n0;
      const T beta = T(p[k]) / iphi / vphi - dp0, 
              alpha = T(1) - beta;
      
      f[k] = alpha * data.scalar[vc0] + beta * data.scalar[vc1];

      v[k][0] = alpha * data.vector[vc0*2] + beta * data.vector[vc1*2];
      v[k][1] = alpha * data.vector[vc0*2+1] + beta * data.vector[vc1*2+1];

      j[k][0][0] = alpha * data.jacobian[vc0*4] + beta * data.jacobian[vc1*4];
      j[k][0][1] = alpha * data.jacobian[vc0*4+1] + beta * data.jacobian[vc1*4+1];
      j[k][1][0] = alpha * data.jacobian[vc0*4+2] + beta * data.jacobian[vc1*4+2];
      j[k][1][1] = alpha * data.jacobian[vc0*4+3] + beta * data.jacobian[vc1*4+3];

#if 0
      T rzp0[3], rzp1[3];
      m3->get_coords_rzp(v2 + dp0 * m2n0 * iphi * vphi, rzp0);
      m3->get_coords_rzp(next + dp1 * m2n0 * iphi * vphi, rzp1);
      if (dp1 == nphi) { // the next plane is in the next period
        // rzp1[2] += nphi * iphi * vphi;
        // fprintf(stderr, "alpha=%f, beta=%f, dp0=%d, dp1=%d, p0=%f, p1=%f\n", alpha, beta, 
        //     dp0, dp1, rzp0[2], rzp1[2]);
      }

      for (int j = 0; j < 3; j ++)
        rzpt[k][j] = alpha * rzp0[j] + beta * rzp1[j];
#endif
    }
#endif
  }

#if 1
  bool b0 = false, b1 = false;
  for (int k = 0; k < n; k ++) {
    if (p[k] == 0)
      b0 = true;
    else if (p[k] == m3->np()-1) 
      b1 = true;
  }
  if (b0 && b1) { // periodical
    // fprintf(stderr, "periodical.\n");
    for (int k = 0; k < n; k ++)
      if (p[k] == 0)
        // zpt[k][2] += m3->np()-1;
        rzpt[k][2] += m3->np();
  }
#endif
}

inline void xgc_blob_filament_tracker::build_critical_surfaces()
{
  fprintf(stderr, "building critical surfaces...\n");

  int i = 0;
  for (auto &kv : intersections) {
    kv.second.id = i ++;
    surfaces.pts.push_back(kv.second);
  }

  int tri_count = 0;
#if 0 // for all 4-simplicies
  for (int timestep = 0; timestep < current_timestep; timestep ++) {
    fprintf(stderr, "pass II, timestep=%d\n", timestep);
    m4->element_for_ordinal(4, timestep, 
        std::bind(&xgc_blob_filament_tracker::check_penta, this, std::placeholders::_1), 
        xl, nthreads, enable_set_affinity);
  }
#else // for all related 4-simplicies
  parallel_for<int>(related_cells, [&](const int e) {
    check_penta(e);
  }, FTK_XL_NONE, nthreads, enable_set_affinity);
#endif

  surfaces.relabel();
  fprintf(stderr, "#pts=%zu, #tri=%zu, #tri_count=%d\n", 
      surfaces.pts.size(), surfaces.tris.size(), tri_count);
}

#if 0 // legacy code, to be removed later
inline /*static*/ std::shared_ptr<xgc_blob_filament_tracker>
xgc_blob_filament_tracker::from_augmented_mesh_file(
    diy::mpi::communicator comm, const std::string &filename)
{
  FILE *fp = fopen(filename.c_str(), "rb");
  assert(fp);
  diy::detail::FileBuffer bb(fp);

  std::shared_ptr<xgc_blob_filament_tracker> tracker(
      new xgc_blob_filament_tracker(comm));

  diy::load(bb, tracker->nphi);
  diy::load(bb, tracker->iphi);
  
  tracker->m2.reset(new simplicial_unstructured_2d_mesh<>);
  diy::load(bb, *tracker->m2);

  tracker->m3.reset(new simplicial_unstructured_3d_mesh<>);
  diy::load(bb, *tracker->m3);

  tracker->m4.reset(new simplicial_unstructured_extruded_3d_mesh<>(*tracker->m3));
  fclose(fp);
  return tracker;
}

inline void xgc_blob_filament_tracker::to_augmented_mesh_file(const std::string& filename)
{
  if (!is_root_proc()) return;

  FILE *fp = fopen(filename.c_str(), "wb");
  assert(fp);
  diy::detail::FileBuffer bb(fp);

  diy::save(bb, nphi);
  diy::save(bb, iphi);
  // diy::save(bb, *m2); // WIP
  // diy::save(bb, *m3); // WIP

  fclose(fp);
}
#endif

inline void xgc_blob_filament_tracker::write_intersections_binary(const std::string& filename) const
{
  if (is_root_proc())
    diy::serializeToFile(intersections, filename); // TODO: TBB
}

inline void xgc_blob_filament_tracker::read_intersections_binary(const std::string& filename)
{
  if (is_root_proc())
    diy::unserializeFromFile(filename, intersections); // TODO: TBB
}

inline void xgc_blob_filament_tracker::write_surfaces(const std::string& filename, std::string format, bool torus) const 
{
  if (!is_root_proc()) return;

  if (ends_with_lower(filename, "vtp")) {
#if FTK_HAVE_VTK
    fprintf(stderr, "writing surfaces...\n");
    write_polydata(filename, transform_vtp_coordinates(surfaces.to_vtp()));
#else
    fatal("FTK not compiled with VTK.");
#endif
  } else {
    surfaces.save(filename, format);
  }
}

inline void xgc_blob_filament_tracker::read_surfaces(const std::string& filename, std::string format)
{
  if (!is_root_proc()) return;

  surfaces.load(filename, format);
  fprintf(stderr, "readed surfaces #pts=%zu, #tris=%zu\n", surfaces.pts.size(), surfaces.tris.size());
}

inline void xgc_blob_filament_tracker::post_process_curves(feature_curve_set_t& curves) const
{
  fprintf(stderr, "before post process: #%zu\n", curves.size());

  const int nphi = m3->get_nphi(), 
            iphi = m3->get_iphi(),
            vphi = m3->get_vphi();
  const int np = m3->np();

  // unwrap and update statistics
  curves.foreach([&](feature_curve_t& curve) {
    curve.unwrap<2>(np);
    curve.update_statistics(); 
  });

  // discard non-poloidal points
  curves.foreach([&](feature_curve_t& c) {
    // c.discard_high_cond();
    c.discard([&](const feature_point_t& p) {
      return !m3->is_poloidal(2, p.tag);
      // return !m3->is_actual_poloidal(2, p.tag); // remove virtual planes as well
    });
  });
  
  // remove trajectories with no points
  curves.filter([&](const feature_curve_t& c) {
    return c.size() > 0;
  });

  // apply only when virtual poloidal planes are used
  if (vphi > 1) {
    // discard loops between non-virtual poloidal planes
#if 1
    curves.filter([&](const feature_curve_t& c) {
      const double delta_phi = c.bbmax[2] - c.bbmin[2];
      // fprintf(stderr, "delta_phi=%f\n", delta_phi);
      return delta_phi > vphi;
    });
#endif

    // TODO: rotate loops to make the first point in an actual poloidal plane
#if 0
    // force phi being monotonous between non-virtual poloidal planes
    curves.foreach([&](feature_curve_t& c) { // TODO: consider loops
      auto poloidals = c.select([&](const feature_point_t& p) { return m3->is_actual_poloidal(2, p.tag); });
      for (int i = 0; i < poloidals.size()-1; i ++) {
        if (c[poloidals[i]].x[2] < c[poloidals[i+1]].x[2]) { // ascending phi
          for (int j = poloidals[i]; j < poloidals[i+1]; j ++) 
            if (c[j].x[2] > c[j+1].x[2]) 
              c[j].x[2] = c[j+1].x[2];
        } else if (c[poloidals[i]].x[2] > c[poloidals[i+1]].x[2]) { // descending phi
          for (int j = poloidals[i]; j < poloidals[i+1]; j ++) 
            if (c[j].x[2] < c[j+1].x[2]) 
              c[j].x[2] = c[j+1].x[2];
        } else { // same phi, flatten things up
          for (int j = poloidals[i]+1; j < poloidals[i+1]; j ++) 
            c[j].x[2] = c[poloidals[i]].x[2];
        }
      }
      // c.loop = false;
    });
#endif
  }
  
  // remove trajectories with no points (again)
  curves.filter([&](const feature_curve_t& c) {
    return c.size() > 0;
  });

  // if the curve is loop and has inconsistent type, rotate the loop
  curves.foreach([&](feature_curve_t& c) { c.rotate(); });

  // smooth types
  const int half_window_size = 2;
  curves.foreach([&](feature_curve_t& c) {
    if (c.size() < half_window_size*2+1) return;

    std::map<int, unsigned int> pending_changes;
    for (int i = half_window_size; i < c.size() - half_window_size; i ++) {
      unsigned int local_consistent_type = c.at(i - half_window_size).type;
      for (int j = i - half_window_size; j <= i + half_window_size; j ++)  {
        if (j == i) continue;
        else if (local_consistent_type != c.at(j).type) {
          local_consistent_type = 0; // type is inconsistent within the window
          break;
        }
      }
        
      if (local_consistent_type != 0 && c.at(i).type != local_consistent_type)
        pending_changes[i] = local_consistent_type;
    }
    for (const auto &kv : pending_changes)
      c.at(kv.first).type = kv.second;
  });

  // print for debugging
#if 0
  curves.foreach([&](const feature_curve_t& c) {
    for (int i = 0; i < c.size(); i ++)
      fprintf(stderr, "r=%f, z=%f, p=%f, type=%d\n", 
          c[i].x[0], c[i].x[1], c[i].x[2], c[i].type);
  });
#endif

  // split into consistent type curves
  curves.split_all();
#if 0
  curves.foreach([&](feature_curve_t& c) {
    c.discard([&](const feature_point_t& p) {
      return !m3->is_actual_poloidal(2, p.tag); // remove virtual planes as well
    });
  });
#endif

  fprintf(stderr, "after post process: #%zu\n", curves.size());
}

inline void xgc_blob_filament_tracker::write_sliced(const std::string& pattern, std::string format, bool torus) const
{
  if (!is_root_proc()) return;

#if FTK_HAVE_VTK
  for (int i = 0; i < end_timestep; i ++) { // TODO
    auto sliced = surfaces.slice_time(i);
    post_process_curves(sliced);
    fprintf(stderr, "sliced timestep %d, #curves=%zu\n", i, sliced.size());

    auto poly = sliced.to_vtp(3, std::vector<std::string>());
    
    const auto filename = ndarray_writer<double>::filename(pattern, i);
    write_polydata(filename, transform_vtp_coordinates(poly));
    // write_polydata(filename, poly);
  }
#else
  fatal("FTK not compiled with VTK.");
#endif
}

inline void xgc_blob_filament_tracker::write_intersections(const std::string& filename, std::string format) const
{
#if FTK_HAVE_VTK
  if (comm.rank() == get_root_proc()) {
    auto poly = get_intersections_vtp();
    write_polydata(filename, poly);
  }
#else
  fatal("FTK not compiled with VTK.");
#endif
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> xgc_blob_filament_tracker::get_intersections_vtp(bool torus) const
{
  const int np = m3->np();

  vtkSmartPointer<vtkPolyData> poly = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  
  vtkIdType pid[1];
  
  // const auto intersections = get_intersections();
  for (const auto &kv : intersections) {
    const auto &cp = kv.second;
    // double p[3] = {cp.x[0], cp.x[1], cp.x[2]}; // rzp coords
    const double phi = cp.x[2] * 2 * M_PI / np;
    const double p[3] = {
      cp.x[0] * cos(phi), 
      cp.x[0] * sin(phi), 
      cp.x[1]
    };
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  poly->SetPoints(points);
  poly->SetVerts(vertices);

  return poly;
}
  
vtkSmartPointer<vtkPolyData> xgc_blob_filament_tracker::transform_vtp_coordinates(vtkSmartPointer<vtkPolyData> poly) const
{
  const int np = m3->np();
  
  vtkSmartPointer<vtkPoints> pts = poly->GetPoints();
  for (auto i = 0; i < pts->GetNumberOfPoints(); i ++) {
    double rzp[3], xyz[3];
    pts->GetPoint(i, rzp);

    const double phi = rzp[2] * 2 * M_PI / np;
    xyz[0] = rzp[0] * cos(phi);
    xyz[1] = rzp[0] * sin(phi);
    xyz[2] = rzp[1];

    pts->SetPoint(i, xyz);
  }
  
  poly->SetPoints(pts);
  return poly;
}

inline vtkSmartPointer<vtkPolyData> xgc_blob_filament_tracker::get_critical_surfaces_vtp(bool torus) const
{
  return surfaces.to_vtp();
}
#endif

} // namespace ftk

#endif
