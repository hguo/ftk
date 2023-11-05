#ifndef _FTK_PARTICLE_TRACER_MPAS_OCEAN_HH
#define _FTK_PARTICLE_TRACER_MPAS_OCEAN_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/particle_tracer.hh>
#include <ftk/filters/mpas_ocean_tracker.hh>
#include <ftk/utils/gather.hh>
#include <ftk/numeric/wachspress_interpolation.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/rad.hh>

#if FTK_HAVE_CUDA
#include "mpas_ocean_particle_tracker.cuh"
#endif

namespace ftk {

struct particle_tracer_mpas_ocean : public particle_tracer, public mpas_ocean_tracker
{
  particle_tracer_mpas_ocean(diy::mpi::communicator comm, 
      std::shared_ptr<mpas_mesh<>> m) : 
    particle_tracer(comm, 3), 
    mpas_ocean_tracker(comm, m), 
    tracker(comm) 
  {
    // this->integrator = PARTICLE_TRACER_INTEGRATOR_SPHERICAL_RK1;
    // this->integrator = PARTICLE_TRACER_INTEGRATOR_SPHERICAL_RK1_WITH_VERTICAL_VELOCITY;
    this->integrator = PARTICLE_TRACER_INTEGRATOR_SPHERICAL_RK4_WITH_VERTICAL_VELOCITY;
  }

  virtual ~particle_tracer_mpas_ocean() {
#if FTK_HAVE_CUDA
    if (xl == FTK_XL_CUDA)
      mop_destroy_ctx(&ctx);
#endif
  }

  void initialize();
  void initialize_particles_at_grid_points(std::vector<int> strides);
  void initialize_particles_latlonz(
      const int nlat, const double lat0, const double lat1,
      const int nlon, const double lon0, const double lon1,
      const int nz,   const double z0,   const double z1);

  static constexpr double earth_radius = 6371229.0;

  void push_field_data_snapshot(std::shared_ptr<ndarray_group> g);
  void prepare_timestep();
  void update_timestep();

protected:
  bool eval_v(int t, const double* x, double *v, int *hint);
  bool eval_v_vertical(int t, const double* x, double *v, int *hint);
  bool eval_v_with_vertical_velocity(int t, const double* x, double *v, int *hint);

  int nch() const { return 7; }
  std::vector<std::string> scalar_names() const { return {"vertVelocity", "salinity", "temperature"}; }

  static double deltaT(std::tm t0, std::tm t1);

protected:
  // std::shared_ptr<ndarray<double>> V[2]; // inherited from particle_tracer
  std::shared_ptr<ndarray<double>> zTop[2];
  std::shared_ptr<ndarray<double>> vertVelocityTop[2];
  std::shared_ptr<ndarray<double>> salinity[2];
  std::shared_ptr<ndarray<double>> temperature[2];
  std::tm timestamp[2];

#if FTK_HAVE_CUDA
  mop_ctx_t *ctx;
#endif
  void load_particles_cuda();
};

////
inline void particle_tracer_mpas_ocean::initialize()
{
  if (xl == FTK_XL_CUDA) {
#if FTK_HAVE_CUDA
    // fprintf(stderr, "loading mesh to gpu...\n");
    mop_create_ctx(&ctx);
    mop_load_mesh(ctx, 
        m->n_cells(),
        m->n_edges(),
        m->n_vertices(),
        m->n_layers(),
        m->max_edges_on_cell(),
        2, // just temperature and salinity for now
        m->xyzCells.data(),
        m->xyzVertices.data(),
        m->nEdgesOnCell.data(),
        m->cellsOnCell.data(),
        m->cellsOnEdge.data(),
        m->cellsOnVertex.data(),
        m->edgesOnCell.data(),
        m->verticesOnCell.data());

    // std::cerr << m->coeffsReconstruct.shape() << std::endl;
    mop_load_e2c_interpolants(ctx, m->coeffsReconstruct.data());
#endif
  }
}

inline void particle_tracer_mpas_ocean::initialize_particles_latlonz(
      const int nlat, const double lat0, const double lat1,
      const int nlon, const double lon0, const double lon1,
      const int nz,   const double z0,   const double z1)
{
  fprintf(stderr, "initializing geo particles: nlat=%d, nlon=%d, nz=%d, lat0=%f, lat1=%f, lon0=%f, lon1=%f, z0=%f, z1=%f\n", 
      nlat, nlon, nz, lat0, lat1, lon0, lon1, z0, z1);

  const double dz   = nz == 1 ? 0.0 : (z1 - z0) / (nz - 1), 
               dlat = nlat == 1 ? 0.0 : (lat1 - lat0) / (nlat - 1),
               dlon = nlon == 1 ? 0.0 : (lon1 - lon0) / (nlon - 1);
  int hint = 0;

  for (int i = 0; i < nlat; i ++) {
    const double lat = deg2rad(i * dlat + lat0);
    const double slat = std::sin(lat), 
                 clat = std::cos(lat);

    for (int j = 0; j < nlon; j ++) {
      const double lon = deg2rad(j * dlon + lon0);
      const double clon = std::cos(lon),
                   slon = std::sin(lon);

      for (int k = 0; k < nz; k ++) {
        const double z = k * dz + z0;
        const double r = earth_radius + z;
       
        feature_curve_t curve;
        curve.id = i*nlat*nz + j*nz + k;

        feature_point_t p;
        p.x[0] = r * clon * clat;
        p.x[1] = r * slon * clat;
        p.x[2] = r * slat;

        p.id = curve.id;
        const int ci = m->locate_cell_i({p.x[0], p.x[1], p.x[2]}, hint);
        p.tag = m->i2cid(ci);
       
#if 0
        fprintf(stderr, "lat=%f, lon=%f, z=%f, x=%f, %f, %f, tag=%llu\n", 
            rad2deg(lat), rad2deg(lon), z, 
            p.x[0], p.x[1], p.x[2], p.tag);
#endif

        if (ci >= 0) {
          curve.push_back(p);
          trajectories.add(curve);
        }
      }
    }
  }

  fprintf(stderr, "#traj=%zu\n", trajectories.size());

  if (xl == FTK_XL_CUDA)
    load_particles_cuda();
}

inline void particle_tracer_mpas_ocean::initialize_particles_at_grid_points(std::vector<int> strides)
{
  int stride_horizontal = 1, stride_vertical = 1;

  if (strides.size() >= 1)
    stride_horizontal = strides[0];

  if (strides.size() >= 2)
    stride_vertical = strides[1];

  // fprintf(stderr, "strides=%d, %d\n", stride_horizontal, stride_vertical);
  // fprintf(stderr, "ncells=%zu, nlayers=%zu\n", m->n_cells(), m->n_layers());

  int nparticles = 0;

  for (auto i = 0; i < m->n_cells(); i += stride_horizontal) {
    double x0[3];
    for (auto k = 0; k < 3; k ++)
      x0[k] = m->xyzCells(k, i);
    const double R = vector_2norm<3>(x0);

    for (auto j = 0; j < m->n_layers(); j += stride_vertical) {
      const double thickness = m->accRestingThickness(j, i);
      const double r = R - thickness;
      for (auto k = 0; k < 3; k ++)
        x0[k] = x0[k] * r / R;

      feature_curve_t curve;
      curve.id = nparticles;
      
      feature_point_t p;
      for (auto k = 0; k < 3; k ++)
        p.x[k] = x0[k]; // m->xyzCells(k, i); 
      p.id = nparticles;
      p.tag = i + 1; // cell id for now for GPU code

      // fprintf(stderr, "cell=%zu, layer=%zu, thickness=%f, x0=%f, %f, %f\n", 
      //     i, j, thickness, x0[0], x0[1], x0[2]);
      curve.push_back(p);
      trajectories.add(curve);

      nparticles ++;
    }
  }
  fprintf(stderr, "#trajectories=%zu\n", trajectories.size());

  if (xl == FTK_XL_CUDA) {
    // fprintf(stderr, "loading particles to gpu...\n");
    load_particles_cuda();
  }
}

inline void particle_tracer_mpas_ocean::load_particles_cuda()
{
#if FTK_HAVE_CUDA
  std::vector<feature_point_lite_t> particles;
  particles.reserve(trajectories.size());
  for (const auto &traj : trajectories)
    particles.push_back(traj.second[0].to_lite());

  mop_load_particles(ctx, particles.size(), particles.data());
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_CUDA);
#endif
}

inline void particle_tracer_mpas_ocean::update_timestep()
{
  prepare_timestep();
  
  typedef std::chrono::high_resolution_clock clock_type;
    
  const int ndays = current_delta_t / 86400.;
  nsteps_per_interval = nsteps_per_day * ndays;

  if (checkpoint_days) 
    nsteps_per_checkpoint = nsteps_per_day * checkpoint_days;
  else if (checkpoint_months)
    nsteps_per_checkpoint = nsteps_per_day * 30 * checkpoint_months;
  else 
    nsteps_per_checkpoint = nsteps_per_interval;
  
  if (xl == FTK_XL_CUDA) {
#if FTK_HAVE_CUDA
    if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d, delta_t=%f\n", 
        current_timestep, current_delta_t);

    current_t = current_timestep;

    bool streamlines = false;
    if (this->ntimesteps == 1) {
      streamlines = true;
      fprintf(stderr, "tracing mpas streamlines...\n");
    }

    // if ((!streamlines) && this->snapshots.size() < 2)
    //   return; // nothing can be done

#if 0
    const int years = 1;
    const double T = 86400.0 * 365 * years;
    const int nsteps = nsteps_per_day * 365 * years;
    const int nsubsteps = nsteps_per_day * 10; // 10 days
#endif

    for (int istep = 0; istep < nsteps_per_interval; istep += nsteps_per_checkpoint) {
      auto t1 = clock_type::now();
      mop_execute(ctx, current_delta_t, nsteps_per_interval, nsteps_per_checkpoint);

      auto t2 = clock_type::now();
      fprintf(stderr, "checkpoint=%d, t_comp=%f, istep=%d, nsteps_per_interval=%d, nsteps_per_checkpoint=%d\n",
          istep / nsteps_per_day,
          std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() * 1e-9, 
          istep, nsteps_per_interval, nsteps_per_checkpoint);

      // push resulting particles to the trajectory
      for (int i = 0; i < ctx->nparticles; i ++) {
        feature_point_lite_t &p = ctx->hparts[i];
        feature_point_t pt(p);

        pt.v[0] = p.scalar[0];
        pt.v[1] = p.scalar[1];
        pt.v[2] = p.scalar[2];
        pt.scalar[0] = p.scalar[3]; // vertVel
        pt.scalar[1] = p.scalar[4]; // salinity
        pt.scalar[2] = p.scalar[5]; // temperature
        pt.id = i;

        auto &traj = trajectories.find(i)->second;
        traj.push_back(pt);
      }
      
      auto t3 = clock_type::now();
      fprintf(stderr, "t_unload=%f\n", 
          std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() * 1e-9);
    }
#endif
  } else 
    particle_tracer::update_timestep();
}

inline void particle_tracer_mpas_ocean::push_field_data_snapshot(std::shared_ptr<ndarray_group> g)
{
  if (xl == FTK_XL_CUDA) {
#if FTK_HAVE_CUDA
    auto t0 = clock_type::now();
    
    // const auto V = g->get_ptr<double>("velocity");
    const auto V = g->get_ptr<double>("normalVelocity");
    const auto zTop = g->get_ptr<double>("zTop");
    const auto vertVelocityTop = g->get_ptr<double>("vertVelocityTop");
    const auto salinity = g->get_ptr<double>("salinity");
    const auto temperature = g->get_ptr<double>("temperature");

    ndarray<double> attrs;
    attrs.reshape(2, m->n_vertices(), m->n_layers());
    for (auto i = 0; i < m->n_layers(); i ++) {
      for (auto j = 0; j < m->n_vertices(); j ++) {
        // data(0, j, i) = V[0]->at(0, j, i);
        // data(1, j, i) = V[0]->at(1, j, i);
        // data(2, j, i) = V[0]->at(2, j, i);
        // data(3, j, i) = zTop[0]->at(0, j, i);
        // data(4, j, i) = vertVelocityTop[0]->at(0, j, i);
        attrs(0, j, i) = salinity->at(0, j, i);
        attrs(1, j, i) = temperature->at(0, j, i);
      }
    }
    mop_load_data_with_normal_velocity(ctx,
        (double)current_timestep,
        V->data(), 
        vertVelocityTop->data(),
        zTop->data(),
        attrs.data());
    
    auto t1 = clock_type::now();
    fprintf(stderr, "t_load=%f\n", 
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9);

#endif
  } 
    
  particle_tracer::push_field_data_snapshot(g);
}

inline double particle_tracer_mpas_ocean::deltaT(std::tm tm0, std::tm tm1)
{
  if (tm0.tm_year < 2000) // a workaround for year before 1970
    tm0.tm_year += 2000;

  if (tm1.tm_year < 2000)
    tm1.tm_year += 2000;

  std::time_t t0 = std::mktime(&tm0),
              t1 = std::mktime(&tm1);

  return std::difftime(t1, t0);
}

inline void particle_tracer_mpas_ocean::prepare_timestep()
{
  V[0] = snapshots[0]->get_ptr<double>("velocity");
  V[1] = snapshots.size() > 1 ? snapshots[1]->get_ptr<double>("velocity") : nullptr;
  
  zTop[0] = snapshots[0]->get_ptr<double>("zTop");
  zTop[1] = snapshots.size() > 1 ? snapshots[1]->get_ptr<double>("zTop") : nullptr;
  
  vertVelocityTop[0] = snapshots[0]->get_ptr<double>("vertVelocityTop");
  vertVelocityTop[1] = snapshots.size() > 1 ? snapshots[1]->get_ptr<double>("vertVelocityTop") : nullptr;

  salinity[0] = snapshots[0]->get_ptr<double>("salinity");
  salinity[1] = snapshots.size() > 1 ? snapshots[1]->get_ptr<double>("salinity") : nullptr;
  
  temperature[0] = snapshots[0]->get_ptr<double>("temperature");
  temperature[1] = snapshots.size() > 1 ? snapshots[1]->get_ptr<double>("temperature") : nullptr;
 
  // timestamp
  for (int i = 0; i < snapshots.size(); i ++)
  {
    const auto xtime = snapshots[i]->get_ptr<char>("xtime");
    const std::string xtimestr(xtime->data(), xtime->size());
    std::istringstream ss(xtimestr);
    ss >> std::get_time(&timestamp[i], "%Y-%m-%d_%H:%M:%S");
  }

  if (snapshots.size() > 1) {
    // fprintf(stderr, "difftime=%f\n", deltaT(timestamp[0], timestamp[1]));
    set_delta_t(deltaT(timestamp[0], timestamp[1]));
  }
}

inline bool particle_tracer_mpas_ocean::eval_v_vertical(int t, const double *x, double *v, int *hint)
{
  static const int max_nverts = 10, max_nlayers = 100;
 
  int cell_i = m->locate_cell_i(x);
  assert(cell_i < m->n_cells());

  int verts_i[max_nverts]; //  = {-1};
  const int nverts = m->verts_i_on_cell_i(cell_i, verts_i);
  const int nlayers = m->n_layers();

  double Xv[max_nverts][3]; // coordinates of vertices
  m->verts_i_coords(nverts, verts_i, Xv);

  double omega[max_nverts] = {0};
  wachspress_weights(nverts, Xv, x, omega);
  // with the wachpress weights, first interpolate the thickness
  // of each layer, and then find the proper layer

  double tops[max_nlayers] = {0};
  for (int l = 0; l < nlayers; l ++)
    for (int i = 0; i < nverts; i ++)
      tops[l] += zTop[t]->at(0, l, verts_i[i]);

  // locate depth layer
  const double z = vector_2norm<3>(x) - earth_radius;
  int layer = -1;
  for (int l = 0; l < nlayers; l ++)
    if (z <= tops[l] && z > tops[l+1]) {
      layer = l;
      break;
    }

  if (layer < 0 || layer == nlayers - 1) {
    // fprintf(stderr, "vertical layer not found, %f, %f, %f\n", x[0], x[1], x[2]);
    return false; 
  } else {
    // fprintf(stderr, "ilayer=%d\n", layer);
  }

  double Vu[max_nverts][3], Vl[max_nverts][3]; // upper/lower velocities on vertices
  for (int i = 0; i < nverts; i ++)
    for (int k = 0; k < 3; k ++) {
      Vu[i][k] = V[t]->at(k, layer, verts_i[i]);
      Vl[i][k] = V[t]->at(k, layer+1, verts_i[i]);
    }
  
  double vu[3] = {0}, vl[3] = {0};
  for (int i = 0; i < nverts; i ++)
    for (int k = 0; k < 3; k ++) {
      vu[k] += omega[i] * Vu[i][k];
      vl[k] += omega[i] * Vl[i][k];
    }

  const double alpha = (z - tops[layer+1]) / (tops[layer] - tops[layer+1]);

  for (int k = 0; k < 3; k ++)
    v[k] = alpha * vu[k] + (1.0 - alpha) * vl[k];

  return true;
}

inline bool particle_tracer_mpas_ocean::eval_v(
    int t, const double *x, double *v, int *hint)
{
  return eval_v_with_vertical_velocity(t, x, v, hint);
  // return eval_v_vertical(t, x, v);

  // 1. find the cell and number of vertices
  // 2. apply wachpress interpolation
  
  int cell_i = m->locate_cell_i(x);
  assert(cell_i < m->n_cells());

  static const int max_nverts = 10;
  int verts_i[max_nverts]; //  = {-1};
  const int nverts = m->verts_i_on_cell_i(cell_i, verts_i);

  // fprintf (stderr, "nverts=%d, verts=%d, %d, %d, %d, %d, %d, %d\n", 
  //     nverts, verts_i[0], verts_i[1], verts_i[2], verts_i[3], verts_i[4], verts_i[5], verts_i[6]);

  double Xv[max_nverts][3]; // coordinates of vertices
  m->verts_i_coords(nverts, verts_i, Xv);

  // for (int i = 0; i < nverts; i ++) 
  //   fprintf(stderr, "x%d=%f, %f, %f\n", i, Xv[i][0], Xv[i][1], Xv[i][2]);

  // std::cerr << V->shape() << std::endl;

  double Vv[max_nverts][3]; // velocities on vertices
  for (int i = 0; i < nverts; i ++)
    for (int k = 0; k < 3; k ++)
      Vv[i][k] = V[t]->at(k, 0, verts_i[i]);
  
  // for (int i = 0; i < nverts; i ++) 
  //   fprintf(stderr, "v%d=%f, %f, %f\n", i, Vv[i][0], Vv[i][1], Vv[i][2]);

  double omega[max_nverts] = {0};
  wachspress_weights(nverts, Xv, x, omega);

  // fprintf(stderr, "omega=%f, %f, %f, %f, %f, %f, %f\n", 
  //     omega[0], omega[1], omega[2], omega[3], omega[4], omega[5], omega[6]);

  v[0] = v[1] = v[2] = 0.0;
  for (int i = 0; i < nverts; i ++)
    for (int k = 0; k < 3; k ++)
      v[k] += omega[i] * Vv[i][k]; // meters per second

  // v[3] = 1.0;

  // fprintf(stderr, "x=%f, %f, %f, %f, cell_i=%d, v=%f, %f, %f\n", 
  //     x[0], x[1], x[2], x[3], cell_i, v[0], v[1], v[2]);

#if 0
  for (size_t k = 0; k < 3; k ++)
    v[k] = V->get(k, 0, cellId) * 1e5;
#endif
  
  // fprintf(stderr, "x=%f, %f, %f, %f, cell_i=%zu, v=%f, %f, %f\n", 
  //     x[0], x[1], x[2], x[3], cell_i, v[0], v[1], v[2]);

  return true;
}

inline bool particle_tracer_mpas_ocean::eval_v_with_vertical_velocity(int t, const double *x, double *v, int *hint)
{
  static const int max_nverts = 10, max_nlayers = 100;
  
  // fprintf(stderr, "t%d, x=%f, %f, %f\n", 
  //     t, x[0], x[1], x[2]);

  int cell_i;
  if (hint) {
    cell_i = m->locate_cell_i(x, hint[0]);
    hint[0]= cell_i;
  }
  else 
    cell_i = m->locate_cell_i(x);

  if (cell_i < 0) // not found
    return false;

  int verts_i[max_nverts]; //  = {-1};
  const int nverts = m->verts_i_on_cell_i(cell_i, verts_i);
  const int nlayers = m->n_layers();

  double Xv[max_nverts][3]; // coordinates of vertices
  m->verts_i_coords(nverts, verts_i, Xv);
 
#if 0 // no need to check point-in-cell anymore.  the locate_cell function will check it
  if (!m->point_in_cell(nverts, Xv, x)) {
    fprintf(stderr, "point x=%f, %f, %f not in cell %d\n", x[0], x[1], x[2], cell_i);
    return false;
  }
#endif

  double omega[max_nverts] = {0};
  wachspress_weights(nverts, Xv, x, omega);
  // with the wachpress weights, first interpolate the thickness
  // of each layer, and then find the proper layer

  // double tops[max_nlayers] = {0};
  bool succ = false;
  const double z = vector_2norm<3>(x) - earth_radius;
  
  int layer = hint[1];
  double upper = 0.0, lower = 0.0;
  int dir; // search dir: 0=up, 1=down

  if (layer >= 0) { // try if the point is still in the layer
    for (int i = 0; i < nverts; i ++) {
      upper += omega[i] * zTop[t]->at(0, layer, verts_i[i]);
      lower += omega[i] * zTop[t]->at(0, layer+1, verts_i[i]);
    }

    if (z > upper)
      dir = 0; // upward
    else if (z <= lower)
      dir = 1; // downward
    else 
      succ = true;
  } else {
    layer = 0;
    dir = 1; // downward
  }

  if (!succ) { // search along the dir
    if (dir == 1) { // downward
      upper = lower;
      for (layer = layer + 1 ; layer < nlayers-1; layer ++) {
        lower = 0.0;
        for (int k = 0; k < nverts; k ++)
          lower += omega[k] * zTop[t]->at(0, layer+1, verts_i[k]);

        if (z <= upper && z > lower) {
          succ = true;
          break;
        } else 
          upper = lower;
      }
    } else { // upward
      lower = upper;
      for (layer = layer - 1; layer >= 0; layer --) {
        upper = 0.0;
        for (int k = 0; k < nverts; k ++)
          upper += omega[k] * zTop[t]->at(0, layer, verts_i[k]);

        if (z <= upper && z > lower) {
          succ = true;
          break;
        } else 
          lower = upper;
      }
    }
  }

  if (!succ) 
    return false;

  hint[1] = layer;

  double Vu[max_nverts][3], Vl[max_nverts][3]; // upper/lower velocities on vertices
  double VVu[max_nverts], VVl[max_nverts]; // upper/lower vertical velocities on vertices
  double Su[max_nverts], Sl[max_nverts]; // salinity
  double Tu[max_nverts], Tl[max_nverts]; // temperature
  // double Su[max_nverts], Sl[max_nverts]; // salinity
  for (int i = 0; i < nverts; i ++) {
    for (int k = 0; k < 3; k ++) {
      Vu[i][k]  = V[t]->at(k, layer, verts_i[i]);
      Vl[i][k]  = V[t]->at(k, layer+1, verts_i[i]);
    }
    VVu[i] = vertVelocityTop[t]->at(0, layer, verts_i[i]);
    VVl[i] = vertVelocityTop[t]->at(0, layer+1, verts_i[i]);
    Su[i] = salinity[t]->at(0, layer, verts_i[i]);
    Sl[i] = salinity[t]->at(0, layer+1, verts_i[i]);
    Tu[i] = temperature[t]->at(0, layer, verts_i[i]);
    Tl[i] = temperature[t]->at(0, layer+1, verts_i[i]);
  }

  double vu[3] = {0}, vl[3] = {0};
  double vvu = 0.0, vvl = 0.0;
  double su = 0.0, sl = 0.0;
  double tu = 0.0, tl = 0.0;
  for (int i = 0; i < nverts; i ++) {
    for (int k = 0; k < 3; k ++) {
      vu[k]  += omega[i] * Vu[i][k];
      vl[k]  += omega[i] * Vl[i][k];
    }
    vvu += omega[i] * VVu[i];
    vvl += omega[i] * VVl[i];
    su += omega[i] * Su[i];
    sl += omega[i] * Sl[i];
    tu += omega[i] * Tu[i];
    tl += omega[i] * Tl[i];
  }

  const double alpha = (z - lower) / (upper - lower);
  const double beta = 1.0 - alpha;
  assert(alpha >= 0.0 && alpha <= 1.0);
  
  // fprintf(stderr, "x=%f, %f, %f, cell=%d, layer=%d, alpha=%f, beta=%f\n", 
  //     x[0], x[1], x[2], cell_i, layer, 
  //     alpha, beta);

  for (int k = 0; k < 3; k ++)
    v[k] = alpha * vu[k] + beta * vl[k];

  // assuming the 4th component is vertical
  for (int k = 0; k < 3; k ++) {
    v[4] = alpha * vvu + beta * vvl;
    v[5] = alpha * su  + beta * sl;
    v[6] = alpha * tu  + beta * tl;
  }

#if 0
  fprintf(stderr, "t%d, cell_i=%d, ilayer=%d x=%f, %f, %f, vx=%f, vy=%f, vz=%f, vv=%f, salinity=%f, temperature=%f\n", 
      t, cell_i, layer, x[0], x[1], x[2], v[0], v[1], v[2], v[4], v[5], v[6]);
#endif

  return true;
}

}

#endif
