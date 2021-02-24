#ifndef _FTK_TEST_CONSTANTS_HH
#define _FTK_TEST_CONSTANTS_HH

#include <ftk/external/json.hh>
#include <ftk/filters/json_interface.hh>

using nlohmann::json;

const int woven_width = 31, woven_height = 37;
const int merger_width = 31, merger_height = 37;
const int tornado_width = 31, tornado_height = 29, tornado_depth = 37;

// for streams
const json js_woven_synthetic = {
  {"type", "synthetic"},
  {"name", "woven"},
  {"dimensions", {woven_width, woven_height}}
};

const json js_woven_synthetic_perturbation = {
  {"type", "synthetic"},
  {"name", "woven"},
  {"dimensions", {128, 128}}
};

const json js_merger_2d_synthetic = {
  {"type", "synthetic"},
  {"name", "merger_2d"}
};

const json js_double_gyre_synthetic = {
  {"type", "synthetic"},
  {"name", "double_gyre"}
};

const json js_moving_extremum_2d_synthetic = {
  {"type", "synthetic"},
  {"name", "moving_extremum_2d"}
};

const json js_moving_extremum_3d_synthetic = {
  {"type", "synthetic"},
  {"name", "moving_extremum_3d"}
};

const json js_moving_ramp_3d_synthetic = {
  {"type", "synthetic"},
  {"name", "moving_ramp_3d"}
};

const json js_moving_extremum_2d_vti = {
  {"format", "vti"},
  {"filenames", "moving_extremum_2d-*.vti"},
};

const json js_volcano_2d_synthetic = {
  {"type", "synthetic"},
  {"name", "volcano_2d"}
};

const json js_tornado_synthetic = {
  {"type", "synthetic"},
  {"name", "tornado"},
  {"dimensions", {tornado_width, tornado_height, tornado_depth}}
};

const json js_woven_float64 = {
  {"type", "file"},
  {"format", "float64"},
  {"filenames", "woven-nonperturbation*.bin"},
  {"dimensions", {woven_width, woven_height}}
};

const json js_woven_nc_unlimited_time = {
  {"type", "file"},
  {"format", "nc"},
  {"filenames", "woven-unlimited-time*.nc"},
  {"variables", {"scalar"}},
};

const json js_woven_nc_no_time = {
  {"type", "file"},
  {"format", "nc"},
  {"filenames", "woven-no-time-*.nc"},
  {"variables", {"scalar"}}
};

const json js_woven_vti {
  {"format", "vti"},
  {"filenames", "woven-*.vti"}
};

const json js_double_gyre_vti {
  {"format", "vti"},
  {"filenames", "double-gyre-*.vti"}
};

// for writers
const json jw_woven_float64 = {
  {"nd", 2},
  {"format", "float64"},
  {"filename", "woven-nonperturbation-%04d.bin"}
};

const json jw_woven_float64_perturbation = {
  {"nd", 2},
  {"format", "float64"},
  {"filename", "woven-perturbation-%04d.bin"}
};

const json jw_tornado_float32 = {
  {"nd", 3},
  {"format", "float32"},
  {"filename", "tornado-%04d.bin"}
};

const json jw_moving_extremum_3d_float32 = {
  {"nd", 3},
  {"format", "float32"},
  {"filename", "moving_extremum_3d-%04d.bin"}
};

const json jw_woven_nc_unlimited_time = {
  {"nd", 2},
  {"format", "nc"},
  {"filename", "woven-unlimited-time-%04d.nc"},
  {"variable", "scalar"}
};

const json jw_tornado_nc = {
  {"nd", 3},
  {"format", "nc"},
  {"filename", "tornado-%04d.nc"},
  {"variable", {"u", "v", "w"}}
};

const json jw_woven_nc_no_time = {
  {"nd", 2},
  {"format", "nc"},
  {"unlimited_time", false},
  {"filename", "woven-no-time-%04d.nc"},
  {"variable", "scalar"}
};

const json jw_woven_vti {
  {"nd", 2},
  {"format", "vti"},
  {"filename", "woven-%04d.vti"},
  {"variable", "scalar"}
};

const json jw_double_gyre_vti {
  {"nd", 2},
  {"format", "vti"},
  {"filename", "double-gyre-%04d.vti"},
  {"variable", "vector"}
};

const json jw_volcano_vti {
  {"nd", 2},
  {"format", "vti"},
  {"filename", "volcano-%04d.vti"},
  {"variable", "scalar"}
};

const json jw_moving_extremum_2d_vti = {
  {"nd", 2},
  {"format", "vti"},
  {"filename", "moving_extremum_2d-%04d.vti"},
  {"variable", "scalar"}
};

const json jw_moving_ramp_3d_vti = {
  {"nd", 3},
  {"format", "vti"},
  {"filename", "moving_ramp_3d_-%04d.vti"},
  {"variable", "scalar"}
};

static std::tuple<size_t, size_t> track_cp2d(const json jstream, const json jconfig = json())
{
  ftk::ndarray_stream<> stream;
  stream.configure(jstream);

  ftk::json_interface consumer;
  consumer.configure(jconfig);
  consumer.consume(stream);
  consumer.post_process();
  consumer.write();
    
  auto tracker = std::dynamic_pointer_cast<ftk::critical_point_tracker_2d_regular>( consumer.get_tracker() );
  auto trajs = tracker->get_traced_critical_points();
  auto points = tracker->get_discrete_critical_points();
  return {trajs.size(), points.size()};
}

#endif
