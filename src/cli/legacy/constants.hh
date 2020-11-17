#ifndef _FTK_CLI_CONSTANTS_HH
#define _FTK_CLI_CONSTANTS_HH

#include <ftk/ndarray/stream.hh>
#include <ftk/utils/string.hh>
#include <ftk/external/cxxopts.hpp>

#define COMMON_OPTS_INPUTS() \
    ("i,input", "Input file name pattern: a single file or a series of file, e.g. 'scalar.raw', 'cm1out_000*.nc'", \
     cxxopts::value<std::string>()) \
    ("f,input-format", "Input file format (auto|float32|float64|nc|h5|vti)", cxxopts::value<std::string>()) \
    ("synthetic", "Use a synthetic case (woven|double_gyre|merger) as inputs", cxxopts::value<std::string>()) \
    ("w,width", "Width", cxxopts::value<size_t>()) \
    ("h,height", "Height", cxxopts::value<size_t>()) \
    ("d,depth", "Depth.  Valid only for 3D data", cxxopts::value<size_t>()) \
    ("n,timesteps", "Number of timesteps", cxxopts::value<size_t>()) \
    ("var", "Variable name(s), e.g. `scalar', `u,v,w'.  Valid only for NetCDF, HDF5, and VTK.", cxxopts::value<std::string>()) \
    ("temporal-smoothing-kernel", "Temporal smoothing kernel bandwidth", cxxopts::value<double>()) \
    ("temporal-smoothing-kernel-size", "Temporal smoothing kernel size", cxxopts::value<size_t>()) \
    ("spatial-smoothing-kernel", "Spatial smoothing kernel bandwidth", cxxopts::value<double>()) \
    ("spatial-smoothing-kernel-size", "Spatial smoothing kernel size", cxxopts::value<size_t>()) \
    ("perturbation", "Gaussian perturbation sigma", cxxopts::value<double>())
    // ("dim", "Spatial dimensionality of data (auto|2|3)", cxxopts::value<std::string>())
    // ("nblocks-per-proc", "Number of blocks per proc", cxxopts::value<int>()) \

static const std::string 
        str_auto("auto"),
        str_none("none"),
        str_zero("0"),
        str_two("2"),
        str_three("3"),
        str_float32("float32"),
        str_float64("float64"),
        str_netcdf("nc"),
        str_hdf5("h5"),
        str_vti("vti"),
        str_vtp("vtp"),
        str_scalar("scalar"),
        str_vector("vector"),
        str_text("text"),
        str_cuda("cuda");

static const std::string
        str_ext_vti(".vti"), // vtkImageData
        str_ext_vtp(".vtp"), // vtkPolyData
        str_ext_netcdf(".nc"),
        str_ext_hdf5(".h5");

static const std::string
        str_critical_point_type_min("min"),
        str_critical_point_type_max("max"),
        str_critical_point_type_saddle("saddle");

static const std::set<std::string>
        set_valid_accelerator({str_none, str_cuda}),
        set_valid_input_format({str_auto, str_float32, str_float64, str_netcdf, str_hdf5, str_vti}),
        set_valid_input_dimension({str_auto, str_two, str_three});

static inline nlohmann::json args_to_json(cxxopts::ParseResult& results)
{
  using nlohmann::json;
  json j;

  if (results.count("input")) {
    const std::string input = results["input"].as<std::string>();
    if (ftk::ends_with(input, ".json")) {
      std::ifstream t(input);
      std::string str((std::istreambuf_iterator<char>(t)),
                       std::istreambuf_iterator<char>());
      t.close();
      return json::parse(str);
    }
    else 
      j["filenames"] = results["input"].as<std::string>();
  }

  if (results.count("synthetic")) {
    j["type"] = "synthetic";
    j["name"] = results["synthetic"].as<std::string>();
  } else 
    j["type"] = "file";

  if (results.count("input-format")) j["format"] = results["input-format"].as<std::string>();
  if (results.count("dim")) j["nd"] = results["dim"].as<std::string>();
  
  std::vector<size_t> dims;
  if (results.count("depth")) {
    dims.resize(3);
    dims[2] = results["depth"].as<size_t>();
    dims[1] = results["height"].as<size_t>();
    dims[0] = results["width"].as<size_t>();
  } else if (results.count("height")) {
    dims.resize(2);
    dims[1] = results["height"].as<size_t>();
    dims[0] = results["width"].as<size_t>();
  }
  if (dims.size())
    j["dimensions"] = dims;

  if (results.count("timesteps")) j["n_timesteps"] = results["timesteps"].as<size_t>();
  if (results.count("var")) {
    const auto var = results["var"].as<std::string>();
    const auto vars = ftk::split(var, ",");
    // if (vars.size() == 1) j["variable"] = var;
    // else if (var.size() > 1) j["variable"] = vars; 
    j["variables"] = vars;
  }

  if (results.count("temporal-smoothing-kernel")) j["temporal-smoothing-kernel"] = results["temporal-smoothing-kernel"].as<double>();
  if (results.count("temporal-smoothing-kernel-size")) j["temporal-smoothing-kernel-size"] = results["temporal-smoothing-kernel-size"].as<size_t>();
  if (results.count("spatial-smoothing-kernel")) j["spatial-smoothing-kernel"] = results["spatial-smoothing-kernel"].as<double>();
  if (results.count("spatial-smoothing-kernel-size")) j["spatial-smoothing-kernel-size"] = results["spatial-smoothing-kernel-size"].as<size_t>();
  if (results.count("perturbation")) j["perturbation"] = results["perturbation"].as<double>();

  return j;
}

#endif
