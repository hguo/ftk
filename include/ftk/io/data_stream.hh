#ifndef _FTK_REGULAR_DATA_FEED_HH
#define _FTK_REGULAR_DATA_FEED_HH

#include <ftk/ndarray.hh>
#include <ftk/io/data_group.hh>
#include <ftk/external/json.hh>

namespace ftk {

using nlohmann::json;

enum {
  SYNTHETIC_NONE,
  SYNTHETIC_WOVEN,
  SYNTHETIC_DOUBLE_GYRE,
  SYNTHETIC_ABC_FLOW
};

struct data_stream {
  data_stream(const json& j_) : j(j_) {}

  virtual void initialize();
  virtual void finalize() {};

  void set_callback(std::function<void(int, data_group*)> f);

  virtual void advance_timestep() {};

  template <typename T> const ndarray<T>& get(const std::string& key, int offset=0) {
    const size_t i = staged_data.size() - offset - 1;
    return staged_data[i]->get<T>();
  }

  void push_timestep(data_group*);
  void pop_timestep();

protected:
  static void fatal(const std::string& str) {
    std::cerr << "FATAL: " << str << std::endl;
    exit(1);
  }

  static void warn(const std::string& str) {
    std::cerr << "WARN: " << str << std::endl;
  }

  static bool ends_with(std::string const & value, std::string const & ending)
  {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
  }

protected:
  json j; // configs, metadata, and everything

  int current_timestep = 0;

protected: 
  std::deque<data_group*> staged_data;
};

//////////////////////////////
void data_stream::initialize()
{
  std::cerr << j << std::endl;

  if (j.contains("type")) {
    if (j["type"] == "synthetic") {
      if (j.contains("name")) {
        if (j["name"] == "woven") {
          j["nd"] = 2;
          if (!j.contains("scaling_factor")) j["scalaring_factor"] = 15.0;
        } else if (j["name"] == "double_gyre") {
          j["nd"] = 2;
        } else if (j["name"] == "merger") {
          j["nd"] = 2;
        } else fatal("synthetic case not available.");
      } else fatal("synthetic case name not given.");
     
#if 0
      // default dimensions
      if (j.contains("width")) assert(j["width"] != 0);
      else j["width"] = 32;
      
      if (j.contains("height")) assert(j["height"] != 0);
      else j["height"] = 32;
        
      if (j["nd"] == 3) {
        if (j.contains("depth")) assert(j["depth"] != 0);
        else j["depth"] = 32;
      }
#endif

      if (j.contains("n_timesteps")) assert(j["n_timesteps"] != 0);
      else j["n_timesteps"] = 32;
    } else if (j["type"] == "file") {
      if (j.contains("filenames")) {
        if (!j["filenames"].is_array()) {
          auto filenames = ftk::ndarray<double>::glob(j["filenames"]);
          if (filenames.empty()) fatal("unable to find matching filename(s).");
          if (j.contains("n_timesteps")) filenames.resize(j["n_timesteps"]);
          else j["n_timesteps"] = filenames.size();
          j["filenames"] = filenames;
        }
        const std::string filename0 = j["filenames"][0];

        if (!j.contains("format")) { // probing file format
          if (ends_with(filename0, "vti")) j["format"] = "vti";
          else if (ends_with(filename0, "nc")) j["format"] = "nc";
          else if (ends_with(filename0, "h5")) j["format"] = "h5";
          else fatal("unabled to determine file format.");
        }
        
        if (j.contains("variables")) { // sanity check of variables
          if (j["variables"].is_array()) {
            for (const auto &v : j["variables"]) {
              if (!v.contains("name")) fatal("missing variable name");
              if (!v["name"].is_string()) fatal("invalid variable name");
              if (!v.contains("components")) fatal("missing variable component list");
              if (!v["components"].is_array()) fatal("variable components must be in an array");
              for (const auto &c : j["components"])
                if (!c.is_string()) fatal("invalid variable component");
            }
          } else fatal("variables must be an array");
        } else fatal("missing variable list");

        if (j["format"] == "float32" || j["format"] == "float64") {
          if (j.contains("nd")) {
            if (j["nd"] != 2 && j["nd"] != 3) fatal("unsupported spatial dimensionality");
          } else fatal("unable to determine spatial dimensionality");

          // if ((j["nd"] == 2 && ((!j.contains("width") || !j.contains("height")))) || 
          //     (j["nd"] == 3 && ((!j.contains("width") || !j.contains("height") || !j.contains("depth")))))
          //   fatal("width, height, and/or depth not specified.");
        } else if (j["format"] == "vti") {

        } else if (j["format"] == "nc") {

        } else if (j["format"] == "h5") {

        }

      } else fatal("missing filenames");
    } else fatal("invalid input type");
  } else fatal("missing `type'");

#if 0
  if (j.contains("synthetic")) {
    if (j.contains("nd")) warn("overriding nd.");

  } else {
    } else {
      if (!j.contains("input_source")) fatal("input_source empty.");
      
    }
      
    const std::string filename0 = j["filenames"][0];
  }
#endif
  std::cerr << j << std::endl;
}







#if 0
struct data_stream_synthetic : public data_stream {
  void set_input_source(const std::string&) {} // nothing todo

  void advance_timestep() {
    data_group *g = data_group::create();
    const ndarray<double>& array = synthesize_timestep(current_timestep);
    g->set(default_variable_name(), array);
    push_timestep(g);
  }

  void set_input_parameters(const std::map<std::string, std::string>& param) {
    data_stream::set_input_parameters(param);
    if (param.find(str_dw) != param.end()) DW = std::stoi(param.at(str_dw));
    if (param.find(str_dh) != param.end()) DH = std::stoi(param.at(str_dh));
    if (param.find(str_dd) != param.end()) DD = std::stoi(param.at(str_dd));
    if (param.find(str_dt) != param.end()) n_timesteps = std::stoi(param.at(str_dt));
    if (param.find(str_time_scale) != param.end()) time_scale = stod(param.at(str_time_scale));
  }

  virtual std::string default_variable_name() const {return std::string();};

  virtual ndarray<double> synthesize_timestep(int t) = 0;

protected:
  size_t DW = 32, DH = 32, DD = 32;
  double time_scale = 1.0;
};

struct data_stream_synthetic_woven : public data_stream_synthetic {
  ndarray<double> synthesize_timestep(int t) {
    return synthetic_woven_2D(DW, DH, t*time_scale);
  };
};

struct data_stream_synthetic_double_gyre : public data_stream_synthetic {
  ndarray<double> synthesize_timestep(int t) {
    return synthetic_double_gyre(DW, DH, t * time_scale);
  };
};

struct data_stream_files : public data_stream {
  void set_input_source(const std::string& pattern) {
    filenames = ndarray<float>::glob(pattern);
    if (n_timesteps == 0) n_timesteps = filenames.size();
    else filenames.resize(n_timesteps);
  }

  void set_input_filenames(const std::vector<std::string>& f) {filenames = f;}

  void load_timestep(int k);

  void advance_timestep();

protected:
  std::vector<std::string> filenames;
};

struct data_stream_raw : public data_stream_files {
  void initialize() {}
};

struct data_stream_vti : public data_stream_files {
  void initialize() {

  }
};

struct data_stream_nc : public data_stream_files {

};

struct data_stream_pnc : public data_stream_files {

};

struct data_stream_factory {
  static data_stream* new_data_stream(int type, const std::string& source);
  static data_stream* new_data_stream(const std::string& source);
  static data_stream* new_synthetic_data_stream(int c) {
    switch (c) {
    case SYNTHETIC_WOVEN: return new data_stream_synthetic_woven;
    case SYNTHETIC_DOUBLE_GYRE: return new data_stream_synthetic_double_gyre;
    default: return NULL;
    }
  }
};
#endif

}

#endif
