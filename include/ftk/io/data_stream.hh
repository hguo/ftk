#ifndef _FTK_REGULAR_DATA_FEED_HH
#define _FTK_REGULAR_DATA_FEED_HH

#include <ftk/ndarray.hh>
#include <ftk/io/data_group.hh>

namespace ftk {

struct data_stream {
  virtual void set_input_source(const std::string&) = 0;
  virtual void set_input_parameters(const std::map<std::string, std::string>& parameters); // key-value pairs for additional parameters
  virtual void initialize();
  virtual void finalize();

  virtual void advance_timestep() = 0;

  template <typename T> const ndarray<T>& get(const std::string& key, int offset=0) {
    const size_t i = staged_data.size() - offset - 1;
    return staged_data[i]->get<T>();
  }

  void push_timestep(data_group*);
  void pop_timestep();

protected:
  int n_timesteps = 0, current_timestep = 0;
  std::deque<data_group*> staged_data;
};

struct data_stream_synthetic : public data_stream {
  void advance_timestep() {
    data_group *g = data_group::create();
    const ndarray<double>& array = synthesize_timestep(current_timestep);
    g->set(default_variable_name(), array);
    push_timestep(g);
  }

  virtual std::string default_variable_name() const {return std::string();};

  virtual ndarray<double> synthesize_timestep(int t) = 0;

protected:
  size_t DW = 0, DH = 0, DD = 0, DT = 0;
};

struct data_stream_synthetic_double_gyre {

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
};

}

#endif
