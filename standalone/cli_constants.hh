#ifndef _FTK_CLI_CONSTANTS_HH
#define _FTK_CLI_CONSTANTS_HH

#include <ftk/ndarray/stream.hh>
#include <ftk/external/cxxopts.hpp>

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

static inline bool ends_with(std::string const & value, std::string const & ending)
{
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

// https://stackoverflow.com/questions/9435385/split-a-string-using-c11
static inline std::vector<std::string> split(const std::string& input, const std::string& regex) {
    // passing -1 as the submatch index parameter performs splitting
    std::regex re(regex);
    std::sregex_token_iterator
        first{input.begin(), input.end(), re, -1},
        last;
    return {first, last};
}

static inline nlohmann::json parse_input_json(cxxopts::ParseResult& results)
{
  using nlohmann::json;
  json j;

  if (results.count("synthetic")) {
    j["type"] = "synthetic";
    j["name"] = results["synthetic"].as<std::string>();
  } else 
    j["type"] = "file";

  if (results.count("input")) j["filenames"] = results["input"].as<std::string>();
  if (results.count("input-format")) j["format"] = results["format"].as<std::string>();
  if (results.count("dim")) j["nd"] = results["dim"].as<std::string>();
  if (results.count("width")) j["width"] = results["width"].as<size_t>();
  if (results.count("height")) j["height"] = results["height"].as<size_t>();
  if (results.count("depth")) j["depth"] = results["depth"].as<size_t>();
  if (results.count("timesteps")) j["n_timesteps"] = results["timesteps"].as<size_t>();
  if (results.count("var")) {
    const auto var = results["var"].as<std::string>();
    const auto vars = split(var, ",");
    if (var.size() == 1) j["variable"] = var;
    else if (var.size() > 1) {
      json jv;
      for (int i = 0; i < vars.size(); i ++)
        jv[i] = vars[i];
      j["variable"] = jv;
    }
  }

  if (results.count("temporal-smoothing-kernel")) j["temporal-smoothing-kernel"] = results["temporal-smoothing-kernel"].as<double>();
  if (results.count("temporal-smoothing-kernel-size")) j["temporal-smoothing-kernel-size"] = results["temporal-smoothing-kernel-size"].as<size_t>();
  if (results.count("spatial-smoothing-kernel")) j["spatial-smoothing-kernel"] = results["spatial-smoothing-kernel"].as<double>();
  if (results.count("spatial-smoothing-kernel-size")) j["spatial-smoothing-kernel-size"] = results["spatial-smoothing-kernel-size"].as<size_t>();

  return j;
}

#endif
