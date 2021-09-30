#ifndef _FTK_XGC_BLOB_FILAMENT_TRACKER_HH
#define _FTK_XGC_BLOB_FILAMENT_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/filters/xgc_tracker.hh>
#include <ftk/features/feature_line.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/clamp.hh>
#include <ftk/ndarray/writer.hh>

#if FTK_HAVE_TBB
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_unordered_set.h>
#endif

#if FTK_HAVE_CUDA
#include "xgc_blob_filament_tracker.cuh"
#endif

namespace ftk {
  
struct xgc_blob_filament_tracker : public xgc_tracker {
  xgc_blob_filament_tracker(diy::mpi::communicator comm,
      std::shared_ptr<simplicial_xgc_3d_mesh<>> mx);
  virtual ~xgc_blob_filament_tracker();

  // xgc_blob_filament_tracker(diy::mpi::communicator comm, 
  //     std::shared_ptr<simplicial_unstructured_2d_mesh<>> m2, 
  //     int nphi_, int iphi_);

  // int cpdims() const { return 3; }
 
  void initialize();
  void reset() {
    field_data_snapshots.clear();
  }
  void update() {}
  void finalize();

public:
  void update_timestep();
  
  void push_field_data_snapshot(std::shared_ptr<ndarray_group> g) { xgc_tracker::push_field_data_snapshot(g); }
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
      T psin[n], // normalized psi
      T f[n], // scalars
      T v[n][2], // vectors
      T j[n][2][2], // jacobians
      T er[n]); 

  std::set<int> get_related_penta_cells(size_t tri) const;
  void add_lite_feature_points(const std::vector<feature_point_lite_t>& pts, bool ordinal);

public:
  void start_lite_feature_point_consumer_threads(int);
  void join_lite_feature_point_consumer_threads();

public:
  void build_critical_line();

public:
  void build_critical_surfaces();
  void post_process_surfaces();
  void post_process_curves(feature_curve_set_t& curves) const;
  void post_analysis_curves(feature_curve_set_t& curves) const; // compute magnetic line alignment, etc.

  void set_enable_post_processing(bool b) { enable_post_processing = b; }
  void set_ordinal_only(bool b) { ordinal_only = b; }
  void set_ordinal_output_pattern(const std::string &str) { ordinal_output_pattern = str; }

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

  const feature_surface_t& get_surfaces() const { return surfaces; }

protected:
#if FTK_HAVE_TBB
  tbb::concurrent_hash_map<int, feature_point_t> intersections;
  tbb::concurrent_unordered_set<int> related_cells;
#else
  std::map<int, feature_point_t> intersections;
  std::set<int> related_cells;
#endif
  feature_surface_t surfaces;

  bool enable_post_processing = true;
  bool ordinal_only = false;

  std::string ordinal_output_pattern;

protected:
  double vector_field_maxabs = 0;
  double vector_field_resolution = std::numeric_limits<double>::max(); // min abs nonzero value of vector field.  for robust cp detection w/o gmp
  double vector_field_scaling_factor = 1;
  
  void update_vector_field_scaling_factor(int minbits=8, int maxbits=16);

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

inline void xgc_blob_filament_tracker::update_vector_field_scaling_factor(int minbits, int maxbits)
{
  // vector_field_resolution = std::numeric_limits<double>::max();
  for (const auto &s : field_data_snapshots) {
    // vector_field_resolution = std::min(vector_field_resolution, s.vector.resolution());
    vector_field_maxabs = std::max(vector_field_maxabs, s.vector.maxabs());
    // std::cerr << s.vector << std::endl;
  }

  vector_field_scaling_factor = std::exp2(-std::ceil(std::log2(vector_field_maxabs)) + 60); // 20 bits

  fprintf(stderr, "maxabs=%f, factor=%f\n", vector_field_maxabs, vector_field_scaling_factor);

#if 0
  int nbits = std::ceil(std::log2(1.0 / vector_field_resolution));
  nbits = std::max(minbits, std::min(nbits, maxbits));

  vector_field_scaling_factor = 1 << nbits;
 
  std::cerr 
    << "maxabs=" << vector_field_maxabs
    << ", resolution=" << vector_field_resolution 
    << ", factor=" << vector_field_scaling_factor 
    << ", nbits=" << nbits << std::endl;
#endif
}

inline void xgc_blob_filament_tracker::initialize()
{
#if FTK_HAVE_CUDA
  if (xl == FTK_XL_CUDA) {
    int device = device_ids.empty() ? 0 : device_ids[0];

    xft_create_ctx(&ctx, device, device_buffer_size_in_mb);
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

    ndarray<double> psin = m2->get_psifield();
    psin /= m2->get_units().psi_x;
    xft_load_psin(ctx, psin.data());
  }
#endif

  // initialize roi
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

inline std::set<int> xgc_blob_filament_tracker::get_related_penta_cells(size_t i) const
{
  std::set<int> my_related_cells;

  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m4 = 
    use_roi ? this->mr4 : this->m4;

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
}

inline void xgc_blob_filament_tracker::start_lite_feature_point_consumer_threads(int nt)
{
  for (int i = 0; i < nt; i ++) {
    // worker_threads.push_back(std::thread([=]() {
    // }));
  }
}

inline void xgc_blob_filament_tracker::add_lite_feature_points(const std::vector<feature_point_lite_t>& pts, bool ordinal) 
{
  // tbb::parallel_for(tbb::blocked_range<size_t>(0, pts.size()), 
  //   [=](const tbb::blocked_range<size_t>& r) {
  //     for (size_t i = r.begin(); i != r.end(); ++ i) {
#if FTK_HAVE_TBB
  fprintf(stderr, "deriving and inserting related cells...\n");
  parallel_for(pts.size(), [=](int i) {
    feature_point_t cp(pts[i]);
    cp.tag += current_timestep * m4->n(2);
    cp.ordinal = ordinal;
    cp.timestep = current_timestep;
    
    std::set<int> related = get_related_penta_cells(cp.tag);

    // fprintf(stderr, "adding tag %llu, ordinal=%d, current_timestep=%d, #related=%zu\n", 
    //     cp.tag, cp.ordinal, cp.timestep, related.size());

    // both intersections and related_cells are concurrent tbb containers
    intersections.insert({cp.tag, cp});
    related_cells.insert( related.begin(), related.end() );
  });
  fprintf(stderr, "done.\n");
#else
  for (auto lcp : pts) {
    feature_point_t cp(lcp);
    cp.tag += current_timestep * m4->n(2);
    cp.ordinal = ordinal;
    cp.timestep = current_timestep;

    intersections.insert({cp.tag, cp});

    if (!ordinal_only) {
      std::set<int> related = get_related_penta_cells(cp.tag);
      related_cells.insert( related.begin(), related.end() );
    }
  }
#endif
};


inline void xgc_blob_filament_tracker::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  update_vector_field_scaling_factor();
 
  auto func = [=](int i) {
    feature_point_t cp;
    if (check_simplex(i, cp)) {
      std::set<int> my_related_cells = get_related_penta_cells(i);

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
    if (current_timestep == end_timestep) 
      xft_swap(ctx); // swap buffers in order to correctly access the very last timestep
    xft_execute(ctx, 1 /* ordinal */, current_timestep);
    std::vector<feature_point_lite_t> results(ctx->hcps, ctx->hcps + ctx->hncps);
    add_lite_feature_points(results, true);

    // interval
    if (field_data_snapshots.size() >= 2) {
      xft_execute(ctx, 2 /* interval */, current_timestep);
      // fprintf(stderr, "** current_timestep=%d, gpu done interval.\n", current_timestep);
      std::vector<feature_point_lite_t> results(ctx->hcps, ctx->hcps + ctx->hncps);
      add_lite_feature_points(results, false);
    }
#else
    fatal("FTK not compiled with CUDA.");
#endif
  } else {
    std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m4 = 
      use_roi ? this->mr4 : this->m4;

    // TODO: strange that m4->element_for has proformance problem...
    // m4->element_for_ordinal(2, current_timestep, func, xl, nthreads, false); // enable_set_affinity);
    const auto nd = m4->n(2), no = m4->n_ordinal(2), ni = m4->n_interval(2);
    // object::parallel_for(no, [=](int i) {func(i + current_timestep * nd);}, FTK_THREAD_PTHREAD, 8, true);
    
    if (ordinal_only) 
      object::parallel_for(no, [=](int i) {func(i);});
    else 
      object::parallel_for(no, [=](int i) {func(i + current_timestep * nd);});

    if (ordinal_only) {
      build_critical_line();
      intersections.clear();
    } else if (field_data_snapshots.size() >= 2)
      object::parallel_for(ni, [=](int i) {func(i + no + current_timestep * nd);});
      //   m4->element_for_interval(2, current_timestep, func, xl, nthreads, false); // enable_set_affinity);
  }
}

inline bool xgc_blob_filament_tracker::check_simplex(int i, feature_point_t& cp)
{
  int tri[3], t[3], p[3];
  double rzpt[3][4], psin[3], f[3], v[3][2], j[3][2][2], er[3];
  __int128 vf[3][2];
  const long long factor = 1 << 15; // WIP
 
  simplex_values<3, double>(i, tri, t, p, rzpt, psin, f, v, j, er);

  // print3x2("rz", rz);
  // print3x2("v", v);

  for (int k = 0; k < 3; k ++) 
    for (int l = 0; l < 2; l ++) 
      // vf[k][l] = vector_field_scaling_factor * v[k][l];
      vf[k][l] = vector_field_scaling_factor * v[k][l];

  // if (vf[0][0] + vf[0][1] + vf[1][0] + vf[1][1] + vf[2][0] + vf[2][1] != 0)
  //   fprintf(stderr, "vf=%lld, %lld; %lld, %lld; %lld, %lld\n", 
  //       vf[0][0], vf[0][1], vf[1][0], vf[1][1], vf[2][0], vf[2][1]);

  bool succ = ftk::robust_critical_point_in_simplex2(vf, tri);
  if (!succ) return false;

  double mu[3], x[4];
  bool succ2 = ftk::inverse_lerp_s2v2(v, mu);
  ftk::clamp_barycentric<3>(mu);

  ftk::lerp_s2v4(rzpt, mu, x);
  for (int k = 0; k < 3; k ++)
    cp.x[k] = x[k];
  cp.t = x[3];

  cp.scalar[0] = lerp_s2(f, mu); // dneOverne0
  cp.scalar[1] = lerp_s2(psin, mu); // psin
  cp.scalar[2] = m2->theta(cp.x[0], cp.x[1]); // theta /// scalar[2] will be the offset to the magnetic lines in post analysis
  cp.scalar[4] = lerp_s2(er, mu);
  // fprintf(stderr, "er=%f, %f, %f\n", er[0], er[1], er[2]);
  
  double h[2][2];
  ftk::lerp_s2m2x2(j, mu, h);
  h[1][0] = h[0][1] = 0.5 * (h[1][0] + h[0][1]);
  cp.type = ftk::critical_point_type_2d(h, true);

  cp.tag = i;
  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m4 = 
    use_roi ? this->mr4 : this->m4;
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
  auto type0 = surfaces.pts[i0].type,
       type1 = surfaces.pts[i1].type,
       type2 = surfaces.pts[i2].type;

  // if (type0 == type1 && type1 == type2 && type0 == CRITICAL_POINT_2D_MAXIMUM)
  //   surfaces.tris.push_back({i0, i1, i2});
  
  {
    std::lock_guard<std::mutex> guard(mutex);
    surfaces.tris.push_back({i0, i1, i2});
  }
}

void xgc_blob_filament_tracker::check_penta(int e)
{
  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m4 = 
    use_roi ? this->mr4 : this->m4;
  
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
    if (!ordinal_only)
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
    T psin[n], // normalized psi
    T f[n],
    T v[n][2],
    T j[n][2][2], 
    T er[n])
{
  const int m2n0 = m2->n(0), m3n0 = m3->n(0);
  const int nphi = m3->get_nphi(), 
            iphi = m3->get_iphi(),
            vphi = m3->get_vphi();
  const int n0ivphi = m2n0 * vphi * nphi;

  if (use_roi)
    mr4->get_simplex(n-1, i, verts);
  else 
    m4->get_simplex(n-1, i, verts);

#if 0
  if (use_roi) {
    mr4->get_simplex(n-1, i, verts);
    for (int k = 0; k < n; k ++) {
      const int rv3 = verts[k] % mr3n0;
      const int p = rv3 / mr2n0;
      const int rv2 = rv3 % mr2n0;
      verts[i] = roi_inverse_node_map[verts[i]]; // FIXME: wrong
    }
  } else 
    m4->get_simplex(n-1, i, verts);
#endif

  for (int k = 0; k < n; k ++) {
    int v3, v2;

    if (use_roi) {
      t[k] = verts[k] / mr4->n(0);
      v3 = verts[k] % mr3->n(0); // needs translation
      v2 = v3 % mr2->n(0); // needs translation
      p[k] = v3 / mr2->n(0);

      v2 = roi_inverse_node_map[ v2 ];
      v3 = m2n0 * p[k] + v2;

      // verts[k] = t[k] * m4->n(0) + v3;
    } else {
      t[k] = verts[k] / m4->n(0); // m4->flat_vertex_time(verts[k]);
      v3 = verts[k] % m3->n(0); // vertex id in m3 mesh
      v2 = v3 % m2->n(0); // vertex id in m2 mesh
      p[k] = v3 / m2->n(0); // poloidal plane coords
    }
    psin[k] = m2->psin(v2); // normalized psi

    // if (m3->is_poloidal(p[k])) fprintf(stderr, "non-poloidal virtual plane\n");
    m3->get_coords_rzp(v3, rzpt[k]);
    // rzpt[k][3] = t[k];
    rzpt[k][2] = rzpt[k][3]; 
    rzpt[k][3] = t[k];

    const int iv = (t[k] == current_timestep) ? 0 : 1; 
    const field_data_snapshot_t &data = field_data_snapshots[iv];

    // m3->interpolate_central_difference(data.scalar, data.vector, data.jacobian, v3, &f[k], v[k], j[k]);
    m3->interpolate(data.scalar, data.vector, data.jacobian, v3, &f[k], v[k], j[k]);
    if (!data.Er.empty())
      er[k] = m3->interpolate(data.Er, v3);
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

inline void xgc_blob_filament_tracker::build_critical_line()
{
#if FTK_HAVE_TBB // TODO
#else
  feature_line_t line;
  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m4 = 
    use_roi ? this->mr4 : this->m4;
  
  int i = 0; 
  // auto lower = intersections.lower_bound(m4->n(2) * current_timestep), 
  //      upper = intersections.upper_bound(m4->n(2) * current_timestep + m4->n_ordinal(2));
  auto lower = intersections.lower_bound(0),
       upper = intersections.upper_bound(m4->n_ordinal(2));

  for (auto it = lower; it != upper; it ++) {
    it->second.id = i ++;
    line.pts.push_back(it->second);
  }

  m4->element_for_ordinal(3, 0, // current_timestep,
    [&](const int tetid) {
      auto sides = m4->sides(3, tetid);
      int count = 0;
      int ids[4]; // some large number;
      for (auto i : sides) {
        if (intersections.find(i) != intersections.end()) {
          ids[count ++] = intersections[i].id;
        }
      }

      if (count == 0) return; 
      else if (count == 2) {
        std::lock_guard<std::mutex> guard(mutex);
        line.edges.push_back(std::array<int, 2>({ids[0], ids[1]}));
      }
      else fprintf(stderr, "irregular count=%d\n", count);
    }, FTK_THREAD_PTHREAD, nthreads, true);

  line.relabel();
  feature_curve_set_t curves = line.to_curve_set();
  fprintf(stderr, "critical line built, #pts=%zu, #edge=%zu, #curves=%zu\n", 
      line.pts.size(), line.edges.size(), curves.size());

  post_process_curves( curves );
  post_analysis_curves( curves );

#if FTK_HAVE_VTK
  auto poly = curves.to_vtp({
      "dneOverne0", "psin", "theta", "offset", "Er"});
  const auto filename = series_filename(ordinal_output_pattern, current_timestep);
  write_polydata(filename, transform_vtp_coordinates(poly));
#endif
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

#if 0 // for all 4-simplicies
  for (int timestep = 0; timestep < current_timestep; timestep ++) {
    fprintf(stderr, "pass II, timestep=%d\n", timestep);
    m4->element_for_ordinal(4, timestep, 
        std::bind(&xgc_blob_filament_tracker::check_penta, this, std::placeholders::_1), 
        xl, nthreads, enable_set_affinity);
  }
#else // for all related 4-simplicies
  parallel_for<int>(related_cells, [=](const int e) {
    check_penta(e);
  }, FTK_XL_NONE, nthreads, true);
#endif

  surfaces.relabel();
  fprintf(stderr, "#pts=%zu, #tri=%zu\n",
      surfaces.pts.size(), surfaces.tris.size());

  post_process_surfaces();
}

inline void xgc_blob_filament_tracker::post_process_surfaces()
{
  // an alternative way to calculate radial velocities:
  // for each ordinal+poloidal point, find all ordinal+poloidal point in 
  // the next timestep at the same poloidal plane (and with the same type), 
  // and then find the nearest.
  // for now, maybe find one with bfs is ok.
  
  fprintf(stderr, "post processing critical surfaces...\n");

  std::shared_ptr<simplicial_xgc_3d_mesh<>> m3 = 
    use_roi ? this->mr3 : this->m3;
  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m4 = 
    use_roi ? this->mr4 : this->m4;
  
  std::vector<std::set<int>> neighbors(surfaces.pts.size());
  for (int i = 0; i < surfaces.tris.size(); i ++) {
    const auto tri = surfaces.tris[i];

    for (int j = 0; j < 3; j ++)
      for (int k = 0; k < 3; k ++)
        if (j != k)
          neighbors[tri[j]].insert(tri[k]);
  }

  parallel_for(surfaces.pts.size(), [=](int i) {
  // for (int i = 0; i < surfaces.pts.size(); i ++) {
    // compute only for poloidal+ordinal elements
    auto &p = surfaces.pts[i];
    const int tp = p.timestep;
    const int pp = m3->get_poloidal(2, p.tag % m4->n(2));

    if (m3->is_poloidal(2, p.tag % m4->n(2)) && p.ordinal) {
      // fprintf(stderr, "investigating %d\n", i);
      std::set<int> visited;
      std::queue<int> Q;
      Q.push(i);
      visited.insert(i);

      while (!Q.empty()) {
        int current = Q.front();
        Q.pop();

        const auto &q = surfaces.pts[current];
        if (current != i &&
            q.ordinal && 
            q.timestep == tp + 1 && 
            // q.type == p.type &&
            m3->is_poloidal(2, q.tag % m4->n(2)) && 
            m3->get_poloidal(2, q.tag % m4->n(2)) == pp)
        {
          p.v[0] = q.scalar[1] - p.scalar[1]; // dpsin
          p.v[1] = q.scalar[2] - p.scalar[2]; // dtheta
          // p.v[0] = q.x[0] - p.x[0];
          // p.v[1] = q.x[1] - p.x[1];
#if 1
          fprintf(stderr, "gotcha: i=%d, j=%d, xi=%f, %f, %f, %f, xj=%f, %f, %f, %f, v=%f, %f\n", 
              i, current,
              p.x[0], p.x[1], p.x[2], p.t, 
              q.x[0], q.x[1], q.x[2], q.t, 
              p.v[0], p.v[1]);
#endif
          break;
        }

        for (auto j : neighbors[current]) {
          if (visited.find(j) == visited.end()
              && surfaces.pts[j].timestep >= tp
              && surfaces.pts[j].timestep <= tp + 1)
          {
            visited.insert(j);
            Q.push(j);
          }
        }
      }
    }
  });


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
    diy::serializeToFile(intersections, filename); 
}

inline void xgc_blob_filament_tracker::read_intersections_binary(const std::string& filename)
{
  if (is_root_proc()) {
    diy::unserializeFromFile(filename, intersections); 
    for (const auto &kv : intersections) {
      std::set<int> related = get_related_penta_cells(kv.second.tag);
      related_cells.insert( related.begin(), related.end() );
    }
  }
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

#if 0
  // simplify
  fprintf(stderr, "simplifying critical surfaces...");
  surfaces.discard([&](const feature_point_t& p) {
    if (p.type != CRITICAL_POINT_2D_MAXIMUM) return true;
    // if (!p.ordinal) return true;
    return false;
  });
  fprintf(stderr, "after simplification: #pts=%zu, #tri=%zu\n",
      surfaces.pts.size(), surfaces.tris.size());
#endif
}

inline void xgc_blob_filament_tracker::post_analysis_curves(feature_curve_set_t& curves) const
{
  // calculate the deviation from magnetic lines, assuming that 
  // curves are already sampled at (virtual) poloidal planes
  curves.foreach([&](feature_curve_t& c) {
    // fprintf(stderr, "curve:\n");
    double accumulated_offset = 0.0;
    for (int i = 0; i < c.size(); i ++) {
      double offset = 0;
      if (i > 0) {
        const feature_point_t &prev = c[i-1], &curr = c[i];
        double rzp[3] = {prev.x[0], prev.x[1], prev.x[2] / m3->np() * M_PI * 2};
        m2->magnetic_map(rzp, curr.x[2] / m3->np() * M_PI * 2);
        double rzp1[3] = {curr.x[0], curr.x[1], curr.x[2] / m3->np() * M_PI * 2};

        offset = vector_dist_2norm_3(rzp, rzp1);
        accumulated_offset += offset;

        // fprintf(stderr, "----projected_rzp=%f, %f, %f, actual_rzp=%f, %f, %f, offset=%f\n", 
        //     rzp[0], rzp[1], rzp[2], rzp1[0], rzp1[1], rzp1[2], offset);
      }

      c[i].scalar[3] = offset;
      // fprintf(stderr, "---r=%f, z=%f, p=%f, t=%f, offset=%f\n", 
      //     c[i].x[0], c[i].x[1], c[i].x[2], c[i].t, offset);
    }

    // fprintf(stderr, "----average_offset=%f\n", accumulated_offset / (c.size() - 1));
  });
}

inline void xgc_blob_filament_tracker::post_process_curves(feature_curve_set_t& curves) const
{
  fprintf(stderr, "before post process: #%zu\n", curves.size());

  std::shared_ptr<simplicial_xgc_3d_mesh<>> m3 = 
    use_roi ? this->mr3 : this->m3;
  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m4 = 
    use_roi ? this->mr4 : this->m4;

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
      return !m3->is_poloidal(2, p.tag % m4->n(2));
      // return !m3->is_actual_poloidal(2, p.tag % m4->n(2)); // remove virtual planes as well
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
    return c.size() > 1;
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
    
  // remove saddle curves
  curves.filter([&](const feature_curve_t& c) {
    // return c.consistent_type == CRITICAL_POINT_2D_MAXIMUM ||
    //        c.consistent_type == CRITICAL_POINT_2D_MINIMUM;
    return c.consistent_type == CRITICAL_POINT_2D_MAXIMUM;
  });
  
  // remove trajectories with no points (again..)
  curves.filter([&](const feature_curve_t& c) {
    return c.size() > 1;
  });

  fprintf(stderr, "after post process: #%zu\n", curves.size());
}

inline void xgc_blob_filament_tracker::write_sliced(const std::string& pattern, std::string format, bool torus) const
{
  if (!is_root_proc()) return;

#if FTK_HAVE_VTK
  for (int i = 0; i < end_timestep; i ++) { // TODO
    fprintf(stderr, "slicing timestep %d\n", i);
    auto sliced = surfaces.slice_time(i);
    if (enable_post_processing) {
      post_process_curves(sliced);
      post_analysis_curves(sliced);
    }
    fprintf(stderr, "sliced timestep %d, #curves=%zu\n", i, sliced.size());

    // auto poly = sliced.to_vtp(3, std::vector<std::string>());
    auto poly = sliced.to_vtp({
        "dneOverne0", "psin", "theta", "offset", "Er"});
    
    const auto filename = series_filename(pattern, i);
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
