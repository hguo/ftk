#ifndef _FTK_XGC_STREAM_HH
#define _FTK_XGC_STREAM_HH

#include <ftk/object.hh>
#include <ftk/ndarray.hh>
#include <ftk/external/json.hh>
#include <ftk/ndarray/ndarray_group.hh>
#include <ftk/mesh/simplicial_xgc_2d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_3d_mesh.hh>

namespace ftk {
using nlohmann::json;
  
struct xgc_stream : public object {
  static std::shared_ptr<xgc_stream> new_xgc_stream(const std::string& path, diy::mpi::communicator comm = MPI_COMM_WORLD);
  xgc_stream(const std::string& path_, diy::mpi::communicator comm = MPI_COMM_WORLD) : path(path_), object(comm) {}
  
  void set_start_timestep(int s) { start_timestep = s; }
  void set_ntimesteps(int n) { ntimesteps = n; }
  void set_vphi(int v) { vphi = v; }
  
  void initialize();
  void probe_nphi_iphi();

  void set_smoothing_kernel_size(double s) { smoothing_kernel_size = s; }
  void set_smoothing_kernel_filename(const std::string f) { smoothing_kernel_filename = f; }
  void set_interpolant_filename(const std::string f) { interpolant_filename = f; }

  void set_callback(std::function<void(int, std::shared_ptr<ndarray_group>)> f) { callback = f; }

  virtual bool advance_timestep() = 0;

  virtual std::string postfix() const = 0; 
  std::string mesh_filename() const { return path + "/xgc.mesh" + postfix(); }
  std::string oneddiag_filename() const { return path + "/xgc.oneddiag" + postfix(); }
  std::string bfield_filename() const { return path + "/xgc.bfield" + postfix(); }
  std::string units_filename() const { return path + "/units.m"; }
  std::string filename(int t) const { return series_filename(path + "/xgc.3d.%05d" + postfix(), t); }

  std::shared_ptr<simplicial_xgc_2d_mesh<>> get_m2() { return m2; }
  std::shared_ptr<simplicial_xgc_3d_mesh<>> get_m3() { return m3; }
  std::shared_ptr<simplicial_xgc_3d_mesh<>> get_mx3() { return mx3; }

protected:
  const std::string path;
  std::string smoothing_kernel_filename, interpolant_filename;

  int nphi = 1, iphi = 1, vphi = 1;
  int start_timestep = 1, current_timestep = 1, ntimesteps = 1;

  double smoothing_kernel_size = 0.03;

  std::shared_ptr<simplicial_xgc_2d_mesh<>> m2;
  std::shared_ptr<simplicial_xgc_3d_mesh<>> m3, mx3;

  std::function<void(int, std::shared_ptr<ndarray_group>)> callback;
};

/////
inline void xgc_stream::initialize()
{
  m2 = simplicial_xgc_2d_mesh<>::from_xgc_mesh_file(mesh_filename(), this->comm);
  m2->initialize_point_locator();

  m2->read_bfield( bfield_filename() );
  m2->read_units_m( units_filename() );
  // m2->read_oneddiag( oneddiag_filename() );

  probe_nphi_iphi();
  m3.reset( new ftk::simplicial_xgc_3d_mesh<>(m2, nphi, iphi) );

  if (vphi == 1) 
    mx3 = m3;
  else 
    mx3.reset( new ftk::simplicial_xgc_3d_mesh<>(m2, nphi, iphi, vphi) );

  // smoothing kernels
  if (file_exists( smoothing_kernel_filename ))
    m2->read_smoothing_kernel( smoothing_kernel_filename );
  else {
    m2->build_smoothing_kernel( smoothing_kernel_size );
    if (!smoothing_kernel_filename.empty())
      m2->write_smoothing_kernel( smoothing_kernel_filename );
  }

  // interpolants
  if (file_exists(interpolant_filename)) 
    mx3->read_interpolants( interpolant_filename );
  else {
    mx3->initialize_interpolants();
    if (!interpolant_filename.empty())
      mx3->write_interpolants( interpolant_filename );
  }
  
  current_timestep = start_timestep;
}

inline void xgc_stream::probe_nphi_iphi()
{
  const auto filename0 = filename(1);

  const auto array_nphi = ndarray<int>::from_file(filename0, "nphi");
  const auto array_iphi = ndarray<int>::from_file(filename0, "iphi");

  nphi = array_nphi[0];
  iphi = std::max(1, array_iphi[0]);
}

} // namespace ftk

#include "xgc_stream_h5.hpp"

namespace ftk {

std::shared_ptr<xgc_stream> xgc_stream::new_xgc_stream(const std::string& path, diy::mpi::communicator comm)
{
  std::shared_ptr<xgc_stream> s;
  if (file_exists( path + "/xgc.mesh.h5" )) s.reset(new xgc_stream_h5(path, comm)); 
  // TODO: adios

  return s;
}

} // namespace ftk

#endif
