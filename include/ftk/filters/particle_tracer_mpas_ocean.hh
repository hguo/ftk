#ifndef _FTK_PARTICLE_TRACER_MPAS_OCEAN_HH
#define _FTK_PARTICLE_TRACER_MPAS_OCEAN_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/particle_tracer.hh>
#include <ftk/filters/mpas_ocean_tracker.hh>
#include <ftk/utils/gather.hh>
#include <ftk/numeric/wachspress_interpolation.hh>
#include <ftk/numeric/vector_norm.hh>

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
    if (xl == FTK_XL_CUDA)
      mop_destroy_ctx(&ctx);
  }

  void initialize();
  void initialize_particles_at_grid_points(std::vector<int> strides);

  static constexpr double earth_radius = 6371229.0;

  void prepare_timestep();

protected:
  bool eval_v(int t, const double* x, double *v, int *hint);
  bool eval_v_vertical(int t, const double* x, double *v, int *hint);
  bool eval_v_with_vertical_velocity(int t, const double* x, double *v, int *hint);

  int nch() const { return 7; }
  std::vector<std::string> scalar_names() const { return {"vertVelocity", "salinity", "temperature"}; }

protected:
  // std::shared_ptr<ndarray<double>> V[2]; // inherited from particle_tracer
  std::shared_ptr<ndarray<double>> zTop[2];
  std::shared_ptr<ndarray<double>> vertVelocityTop[2];
  std::shared_ptr<ndarray<double>> salinity[2];
  std::shared_ptr<ndarray<double>> temperature[2];

  mop_ctx_t *ctx;
};

////
inline void particle_tracer_mpas_ocean::initialize()
{
  if (xl == FTK_XL_CUDA) {
    fprintf(stderr, "loading mesh to gpu...\n");
    mop_create_ctx(&ctx);
    mop_load_mesh(ctx, 
        m->n_cells(),
        m->n_layers(),
        m->n_vertices(),
        m->max_edges_on_cell(),
        nch(),
        m->xyzCells.data(),
        m->xyzVertices.data(),
        m->nEdgesOnCell.data(),
        m->cellsOnCell.data(),
        m->verticesOnCell.data());
  }
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

      // fprintf(stderr, "cell=%zu, layer=%zu, thickness=%f, x0=%f, %f, %f\n", 
      //     i, j, thickness, x0[0], x0[1], x0[2]);
      curve.push_back(p);
      trajectories.add(curve);

      nparticles ++;
    }
  }
  
  fprintf(stderr, "#trajectories=%zu\n", trajectories.size());
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
