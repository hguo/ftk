#ifndef _FTK_TEST_CONSTANTS_HH
#define _FTK_TEST_CONSTANTS_HH

#include <ftk/external/json.hh>

using nlohmann::json;

// for streams
const json js_woven_synthetic = {
  {"type", "synthetic"},
  {"name", "woven"}
};

const json js_woven_float64 = {
  {"type", "file"},
  {"format", "float64"},
  {"filenames", "woven-*.bin"},
  {"width", 32},
  {"height", 32}
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

// for writers
const json jw_woven_float64 = {
  {"nd", 2},
  {"format", "float64"},
  {"filename", "woven-%04d.bin"}
};

const json jw_woven_nc_unlimited_time = {
  {"nd", 2},
  {"format", "nc"},
  {"filename", "woven-unlimited-time-%04d.nc"},
  {"variable", "scalar"}
};

const json jw_woven_nc_no_time = {
  {"nd", 2},
  {"format", "nc"},
  {"unlimited_time", false},
  {"filename", "woven-no-time-%04d.nc"},
  {"variable", "scalar"}
};

#endif
