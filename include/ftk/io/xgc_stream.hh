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

// step: time index in the simulation
// timestep: index of all available steps

struct xgc_stream : public object {
  static std::shared_ptr<xgc_stream> new_xgc_stream(const std::string& path, diy::mpi::communicator comm = MPI_COMM_WORLD);
  xgc_stream(const std::string& path_, diy::mpi::communicator comm = MPI_COMM_WORLD) : path(path_), object(comm) {}
  
  void set_start_timestep(int s) { start_timestep = s; }
  void set_ntimesteps(int n) { ntimesteps = n; }
  void set_vphi(int v) { vphi = v; }

  void set_enable_initialize_smoothing_kernel(bool b) { enable_initialize_smoothing_kernel = b; }
  void set_enable_initialize_interpolants(bool b) { enable_initialize_interpolants = b; }

  void initialize();
  void probe_nphi_iphi();

  void set_smoothing_kernel_size(double s) { smoothing_kernel_size = s; }
  void set_smoothing_kernel_filename(const std::string f) { smoothing_kernel_filename = f; }
  void set_interpolant_filename(const std::string f) { interpolant_filename = f; }

  void set_callback(std::function<void(int, std::shared_ptr<ndarray_group>)> f) { callback = f; }

  virtual std::shared_ptr<ndarray_group> request_step(int i) = 0;

  virtual bool read_oneddiag() = 0;
  virtual bool advance_timestep() = 0;

  virtual bool read_units();

  virtual std::string postfix() const = 0; 
  std::string mesh_filename() const { return path + "/xgc.mesh" + postfix(); }
  std::string oneddiag_filename() const { return path + "/xgc.oneddiag" + postfix(); }
  std::string bfield_filename() const { return path + "/xgc.bfield" + postfix(); }
  std::string units_filename() const { return path + "/units.m"; }
  std::string filename_step(int step) const { return series_filename(path + "/xgc.3d.%05d" + postfix(), step); }
  std::string filename(int t) const { 
    // fprintf(stderr, "%zu\n", filenames.size());
    if (filenames.empty()) return series_filename(path + "/xgc.3d.%05d" + postfix(), std::max(t, 1));
    else return filenames[t];
  }

  static int filename2step(const std::string&);

  std::shared_ptr<simplicial_xgc_2d_mesh<>> get_m2() { return m2; }
  std::shared_ptr<simplicial_xgc_3d_mesh<>> get_m3() { return m3; }
  std::shared_ptr<simplicial_xgc_3d_mesh<>> get_mx3() { return mx3; }

  const std::set<int>& get_available_steps() const { return available_steps; }

protected:
  bool enable_initialize_smoothing_kernel = true, 
       enable_initialize_interpolants = true;

protected:
  ndarray<int> steps;
  ndarray<double> time, Te1d, psi_mks, 
    e_gc_density_avg, 
    e_perp_temperature_avg, 
    e_parallel_mean_en_avg; // 1d

  std::set<int> available_steps;
  std::map<int, int> step2istep; // in oneddiag

protected:
  const std::string path;
  std::vector<std::string> filenames;

  std::string smoothing_kernel_filename, interpolant_filename;

  int nphi = 1, iphi = 1, vphi = 1;
  int start_timestep = 1, current_timestep = 1, ntimesteps = 0;

  double smoothing_kernel_size = 0.03;

  std::shared_ptr<simplicial_xgc_2d_mesh<>> m2;
  std::shared_ptr<simplicial_xgc_3d_mesh<>> m3, mx3;

  std::function<void(int, std::shared_ptr<ndarray_group>)> callback;
};

/////
inline int xgc_stream::filename2step(const std::string& str)
{
  size_t last = str.find_last_of("xgc.3d.");
  std::string sub = str.substr(last-5, 5);
  
  int i = -1;
  try {
    i = std::stoi(sub);
  } catch (...) {
    return -1;
  }

  return i;
}

inline bool xgc_stream::read_units()
{
  return m2->read_units_m( units_filename() );
}

inline void xgc_stream::initialize()
{
  const bool has_oneddiag = read_oneddiag();
    
  filenames = glob(path + "/xgc.3d.*");
  for (int i = 0; i < filenames.size(); i ++) {
    int s = filename2step(filenames[i]);
    if (s >= 0)
      available_steps.insert(s);
  }

  if (ntimesteps > 0)
    ntimesteps = std::min(size_t(ntimesteps), filenames.size());
  else 
    ntimesteps = filenames.size();
  
  m2 = simplicial_xgc_2d_mesh<>::from_xgc_mesh_file(mesh_filename(), this->comm);
  m2->initialize_point_locator();

  m2->read_bfield( bfield_filename() );
  
  read_units();

  // m2->read_oneddiag( oneddiag_filename() );

  probe_nphi_iphi();
  fprintf(stderr, "nphi=%d, iphi=%d, vphi=%d\n", nphi, iphi, vphi);
  m3.reset( new ftk::simplicial_xgc_3d_mesh<>(m2, nphi, iphi) );

  if (vphi == 1) 
    mx3 = m3;
  else 
    mx3.reset( new ftk::simplicial_xgc_3d_mesh<>(m2, nphi, iphi, vphi) );

  // smoothing kernels
  if (enable_initialize_smoothing_kernel) {
    if (file_exists( smoothing_kernel_filename ))
      m2->read_smoothing_kernel( smoothing_kernel_filename );
    else {
      m2->build_smoothing_kernel_cached( smoothing_kernel_size );
      if (!smoothing_kernel_filename.empty())
        m2->write_smoothing_kernel( smoothing_kernel_filename );
    }
  }

  // interpolants
  if (vphi > 1 && enable_initialize_interpolants) {
    if (file_exists(interpolant_filename)) 
      mx3->read_interpolants( interpolant_filename );
    else {
      mx3->initialize_interpolants_cached();
      if (!interpolant_filename.empty())
        mx3->write_interpolants( interpolant_filename );
    }
  }
  
  current_timestep = start_timestep;
}

inline void xgc_stream::probe_nphi_iphi()
{
  // try if we can find the information in xgc input file
  std::ifstream ifs(path + "/input");
  if (ifs.is_open()) {
    std::string str;
    std::string varname, dummy;
    double value;
    while (std::getline(ifs, str)) {
      if (str.find("&") == 0) // section begins
        continue;
      else if (str.find("/") == 0) // section ends
        continue;
      else if (str.find("!") == 0) // comment at the beginning
        continue; 
      else if (str.find("=") != std::string::npos) { // kv pair
        // fprintf(stderr, "----%s\n", str.c_str());
        if (str.find("!") != std::string::npos)
          str = str.substr(0, str.find("!")); // remove comment in the end of the line

        str.erase(std::remove(str.begin(), str.end(), ' '), str.end()); // remove spaces
        auto strs = split(str, "=");

        try {
          const auto var = strs[0];
          double val = std::stod(strs[1]); // note that the val may not be numeric

          if (var == "sml_wedge_n") {
            fprintf(stderr, "sml_wedge_n=%f\n", val);
            iphi = val;
          } else if (var == "sml_nphi_total") {
            fprintf(stderr, "sml_nphi_total=%f\n", val);
            nphi = val;
          }
        } catch (...) {
          continue;
        }

      }
    }
  } else {
    iphi = m2->get_units().sml_wedge_n;
  
    // probe nphi from the first availabe timestep
    const auto filename0 = filename(0);

    const auto array_nphi = ndarray<int>::from_file(filename0, "nphi");
    // const auto array_iphi = ndarray<int>::from_file(filename0, "iphi"); // iphi in xgc outputs are zero nowadays

    nphi = array_nphi[0];
    // iphi = std::max(1, array_iphi[0]);
  }
}

} // namespace ftk

#include "xgc_stream_h5.hpp"
#include "xgc_stream_adios2.hpp"

namespace ftk {

std::shared_ptr<xgc_stream> xgc_stream::new_xgc_stream(const std::string& path, diy::mpi::communicator comm)
{
  std::shared_ptr<xgc_stream> s;
  if (file_exists( path + "/xgc.mesh.h5" )) s.reset(new xgc_stream_h5(path, comm));
  else if (is_directory( path + "/xgc.mesh.bp" )) s.reset(new xgc_stream_adios2(path, comm));
  else fatal("cannot find xgc.mesh file");

  return s;
}

} // namespace ftk

#endif
