#ifndef _FTK_TEST_CONSTANTS_HH
#define _FTK_TEST_CONSTANTS_HH

#include <ftk/external/json.hh>

using nlohmann::json;

const int woven_width = 31, woven_height = 37;
const int tornado_width = 31, tornado_height = 29, tornado_depth = 37;

// for streams
const json js_woven_synthetic = {
  {"type", "synthetic"},
  {"name", "woven"},
  {"width", woven_width},
  {"height", woven_height}
};

const json js_tornado_synthetic = {
  {"type", "synthetic"},
  {"name", "tornado"},
  {"width", tornado_width},
  {"height", tornado_height},
  {"depth", tornado_depth}
};

const json js_woven_float64 = {
  {"type", "file"},
  {"format", "float64"},
  {"filenames", "woven-*.bin"},
  {"width", woven_width},
  {"height", woven_height}
};

const json js_woven_nc_unlimited_time = {
  {"type", "file"},
  {"format", "nc"},
  {"filenames", "woven-unlimited-time*.nc"},
  {"variable", "scalar"},
  {"nd", 2}
};

const json js_woven_nc_no_time = {
  {"type", "file"},
  {"format", "nc"},
  {"filenames", "woven-no-time-*.nc"},
  {"variable", "scalar"}
};

const json js_woven_vti {
  {"format", "vti"},
  {"filenames", "woven-*.vti"}
};

// for writers
const json jw_woven_float64 = {
  {"nd", 2},
  {"format", "float64"},
  {"filename", "woven-%04d.bin"}
};

const json jw_tornado_float32 = {
  {"nd", 3},
  {"format", "float32"},
  {"filename", "tornado-%04d.bin"}
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

#endif
