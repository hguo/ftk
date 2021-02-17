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

struct json_interface {
  // json options:
  // - feature_type, string, required: 
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
  // - input, json, see ndarray/stream.hh for spec
  // - output_type, string, by default "traced": intersections, traced, sliced, or intercepted
  // - output_format, string, by default "text": text, json, vtp, vtu, ply
  // - output_pattern, string, required
  // - accelerator, string, by default "none": none, cuda, or hipsycl
  // - nthreads, int, by default 0: number of threads; 0 will lead to the max available number of CPUs
  void configure(const json& j);
};

}

#endif
