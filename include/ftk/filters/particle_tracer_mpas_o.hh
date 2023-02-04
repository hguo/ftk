#ifndef _FTK_PARTICLE_TRACER_MPAS_O_HH
#define _FTK_PARTICLE_TRACER_MPAS_O_HH

#include <ftk/config.hh>
#include <ftk/filters/particle_tracer.hh>
#include <ftk/mesh/simplicial_mpas_2d_mesh.hh>

namespace ftk {

struct particle_tracer_mpas_o : public tracker
{
  particle_tracer_mpas_o(diy::mpi::communicator comm, 
      std::shared_ptr<simplicial_mpas_2d_mesh<>> m_) : tracker(comm), m(m_) {}
  virtual ~particle_tracer_mpas_o() {}

  void initialize();
  void update() {}
  void finalize() {}
  bool advance_timestep();

  void update_timestep();

  void push_field_data_snapshot(const ndarray<double>& uvw);
  bool pop_field_data_snapshot();

protected:
  bool eval_v(const double *x, double *v) const;

protected:
  std::shared_ptr<simplicial_mpas_2d_mesh<>> m;

protected:
  struct field_data_snapshot_t {
    ndarray<double> uvw;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;
};

//////////

inline void particle_tracer_mpas_o::push_field_data_snapshot(const ndarray<double>& arr)
{
  field_data_snapshot_t snapshot;
  snapshot.uvw = arr;

  field_data_snapshots.emplace_back(snapshot);
}

inline bool particle_tracer_mpas_o::pop_field_data_snapshot()
{
  if (field_data_snapshots.size() > 0) {
    field_data_snapshots.pop_front();
    return true;
  } else return false;
}

inline bool particle_tracer_mpas_o::eval_v(const double *x, double *v) const
{
  return false;
}

inline void particle_tracer_mpas_o::initialize()
{
  // initialize a seed per vertex
  for (size_t i = 0; i < m->n(0); i ++) {
    double coords[3];
    m->get_coords(i, coords);
    // fprintf(stderr, "i=%zu, %f, %f, %f\n", i,
    //     coords[0], coords[1], coords[2]);
  }
}

inline bool particle_tracer_mpas_o::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}

inline void particle_tracer_mpas_o::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
}

}

#endif
