#ifndef _FTK_MPAS_STREAM_HH
#define _FTK_MPAS_STREAM_HH

#include <ftk/object.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/ndarray_group.hh>

namespace ftk {

using nlohmann::json;

struct mpas_stream : public object {
  mpas_stream(const std::string path_, diy::mpi::communicator comm_=MPI_COMM_WORLD) :
    path(path_), comm(comm_) {}
  
  void initialize();
  void set_callback(std::function<void(int, std::shared_ptr<ndarray_group>)> f) { callback = f; }

  std::shared_ptr<mpas_mesh<>> mesh() { return m; }

  void set_ntimesteps(int n) { ntimesteps = n; }
  bool advance_timestep();

public:
  diy::mpi::communicator comm;
  std::string path;
  
  std::shared_ptr<mpas_mesh<>> m;
  std::function<void(int, std::shared_ptr<ndarray_group>)> callback;
  
  int ncid;
  size_t start_timestep = 0, current_timestep = 0, ntimesteps = 0;
};

//////
void mpas_stream::initialize()
{
  m.reset(new mpas_mesh<>);
  m->read_netcdf(path);
  m->initialize();
  m->initialize_c2v_interpolants();

#if FTK_HAVE_NETCDF
#if NC_HAS_PARALLEL
  int rtn = nc_open_par(path.c_str(), NC_NOWRITE, comm, MPI_INFO_NULL, &ncid);
  if (rtn != NC_NOERR)
    NC_SAFE_CALL( nc_open(path.c_str(), NC_NOWRITE, &ncid) );
#else
  NC_SAFE_CALL( nc_open(path.c_str(), NC_NOWRITE, &ncid) );
#endif

  int dim_unlimited;
  NC_SAFE_CALL( nc_inq_unlimdim(ncid, &dim_unlimited) );
  NC_SAFE_CALL( nc_inq_dimlen(ncid, dim_unlimited, &ntimesteps) );

  fprintf(stderr, "mpas_ntimesteps=%zu\n", ntimesteps);

#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

bool mpas_stream::advance_timestep()
{
  if (current_timestep >= ntimesteps)
    return false;

  // fprintf(stderr, "current_timestep=%zu\n", current_timestep);
  std::shared_ptr<ndarray_group> g(new ndarray_group);

  {
    const size_t st[3] = {current_timestep, 0, 0}, 
                 sz[3] = {1, m->n_cells(), m->n_layers()};
    
    ndarray<double> velocityX, velocityY, velocityZ; 
    velocityX.read_netcdf(ncid, "velocityX", st, sz);
    velocityY.read_netcdf(ncid, "velocityY", st, sz);
    velocityZ.read_netcdf(ncid, "velocityZ", st, sz);
    // fprintf(stderr, "vectors read.\n");
  
    ndarray<double> velocity = ndarray<double>::concat({velocityX, velocityY, velocityZ});
    g->set("velocity", m->interpolate_c2v(velocity)); // vertexwise velocity
    // fprintf(stderr, "vectors interpolated.\n");
  
    ndarray<double> salinity; 
    salinity.read_netcdf(ncid, "salinity", st, sz);
    salinity.make_multicomponents();
    g->set("salinity", m->interpolate_c2v(salinity));

    ndarray<double> layerThickness;
    layerThickness.read_netcdf(ncid, "layerThickness", st, sz);
    layerThickness.make_multicomponents();
    g->set("layerThickness", m->interpolate_c2v(layerThickness));

    ndarray<double> zTop;
    zTop.read_netcdf(ncid, "zTop", st, sz);
    zTop.make_multicomponents();
    g->set("zTop", m->interpolate_c2v(zTop));
    
    // fprintf(stderr, "scalars read.\n");
  }

  {
    const size_t st[3] = {current_timestep, 0, 0}, 
                 sz[3] = {1, m->n_cells(), m->n_layers()+1};
    
    ndarray<double> vertVelocityTop;
    vertVelocityTop.read_netcdf(ncid, "vertVelocityTop", st, sz);
    vertVelocityTop.make_multicomponents();
    g->set("vertVelocityTop", m->interpolate_c2v(vertVelocityTop));
  }
  
  // fprintf(stderr, "callback..\n");
  callback(current_timestep, g);

  current_timestep ++;
  return true;
}

}

#endif
