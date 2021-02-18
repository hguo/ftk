#ifndef _FTK_TRACKER_JSON_INTERFACE_HH
#define _FTK_TRACKER_JSON_INTERFACE_HH

#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/filters/critical_point_tracker_2d_unstructured.hh>
#include <ftk/filters/critical_point_tracker_3d_unstructured.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_2d_mesh.hh>
#include <ftk/ndarray/stream.hh>
#include <ftk/ndarray/writer.hh>
#include <ftk/io/util.hh>

namespace ftk {

struct json_interface : public object {
  json_interface(diy::mpi::communicator comm) : object(comm) {}

  // json options:
  // - feature, string, required: possible values include
  //      critical_point|cp: 2D/3D critical points
  //      critical_line|cl: 3D critical line (intersections of two isosurfaces)
  //      isosurface|iso: 3D isosurface
  //      tdgl_vortex|tdgl: magnetic flux vortices in 3D time-dependent Ginzburg-Landau (TDGL) simulations
  //      sujudi_haimes|sh: Sujudi-Haimes vortices in 3D vector fields
  //      levy_degani_seginer|lds: Levy-Degani-Seginer vortices in 3D vector fields
  //      connected_component|cc: 2D/3D connected components
  //      xgc_blob_filament_2d: 2D XGC blob filaments defined by local extrema
  //      xgc_blob_filament: 3D XGX blob filaments defined by local extrema
  //      xgc_blob_threshold: 3D XGC blob filaments defined by levelsets
  // - input, json, see ndarray/stream.hh for details
  // // - writeback, json, optional: write input data back to files; see ndarray/writer.hh for details
  // - archive, json, optional:
  //    - discrete, string, optional: file name to load/store discrete feature points w/o tracking
  //    - traced, string, optional: file name to load/store traced features
  // - output, json, required:
  //    - type, string, by default "traced": intersections, traced, sliced, or intercepted
  //    - format, string, by default "auto": auto, text, json, vtp, vtu, ply
  //    - pattern, string, required: e.g. "surface.vtp", "sliced-%04d.vtp"
  // - threshold, number, by default 0: threshold for some trackers, e.g. contour trackers
  // - accelerator, string, by default "none": none, cuda, or hipsycl
  // - thread_model, string by default "pthread": pthread, openmp, or tbb
  // - nblocks, int, by default 0: number of blocks; 0 will be replaced by the number of processes
  // - nthreads, int, by default 0: number of threads; 0 will replaced by the max available number of CPUs
  // - enable_streaming, bool, by default false
  // - enable_discarding_interval_points, bool, by default false
  // - enable_fast_detection, bool, by default true
  // - enable_deriving_velocities, bool, by default false
  // - enable_post_processing, bool, by default true
  // - xgc, json, optional: XGC-specific options
  //    - format, string, by default auto: auto, h5, or bp
  //    - path, string, optional: XGC data path, which contains xgc.mesh, xgc.bfield, units.m, 
  //      and xgc.3d data
  //    - mesh, string, optional: path to xgc.mesh
  //    - bfield, string, optional: path to xgc.bfield
  //    - oneddiag, string, optional: path to xgc.onediag
  //    - units.m, string, optional: path to units.m
  //    - vphi, integer, by default 1: number of virtual poloidal planes
  //    - interpolant, string, optional: path to interpolant file
  //    - smoothing_kernel_file, string, optional: path to smoothing kernel file; will (over)write 
  //      the file if the file does not exist or the file has a different smoothing kernel size
  //    - smoothing_kernel_size, number, optional: smoothing kernel size
  void configure(const json& j); // configure and validate json input

protected:
  json j;
  std::shared_ptr<tracker> tr;
  std::shared_ptr<ndarray_stream> st;

  int feature_type = 0;
};

///////
void json_interface::configure(const json& j)
{
  if (j.contains("feature")) {
    if (j["feature"].is_string()) {
      feature_type = tracker::str2tracker(j["feature"]);
      if (feature_type == 0) fatal("invalid feature");
      // otherwise ok
    }
  } else fatal("missing feature");

  if (j.contains("input")) {
    st.reset(new ndarray_stream<>);
    st->set_input_source_json(j_input);
  } else fatal("missing input");

  if (j.contains("archive")) {
    if (j["archive"].is_object()) {
      json ja = j["archive"];
      if (ja.contains("discrete")) {
        if (ja["discrete"].is_string()) { // OK
        } else fatal("invalid discrete archival file name");
      }

      if (ja.contains("traced")) {
        if (ja["traced"].is_string()) { // OK
        } else fatal("invalid traced archival file name");
      }
    }
  }

  if (j.contains("output")) {
    if (j["output"].is_object()) {
      json jo = j["output"];
      if (jo.contains("type")) {
        if (jo["type"].is_string()) { // OK
          const std::string t = jo["type"];
          // TODO: check if type is valid
        } else fatal("invalid output type");
      } else jo["type"] = "auto";
    }
  }

  if (j.contains("threshold")) {
    if (j["threshold"].is_number()) { // OK
    } else fatal("invalid threshold");
  } else j["threshold"] = 0.0;

  if (j.contains("accelerator")) {
    if (j["accelerator"].is_string()) { // OK
      const std::string a = j["accelerator"];

  }
}

}

#endif
