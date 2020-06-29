#ifndef _FTK_REGULAR_DATA_FEED_HH
#define _FTK_REGULAR_DATA_FEED_HH

#include <ftk/ndarray.hh>

namespace ftk {

struct data_stream {
  virtual void set_input_source(const std::string&);
  virtual void initialize();
  virtual void finalize();

  virtual void advance_timestep() = 0;

  std::map<std::string, ndarray<double>> get_current_timestep(); // <key, array>
  std::map<std::string, std::vector<size_t>> shapes; // <key, shape>
  // std::map<std::string, 

protected:
  int n_timesteps = 0;
};

struct data_stream_synthetic : public data_stream {

};

struct data_stream_file : public data_stream {
  void set_input_source(const std::string& pattern) {
    filenames = ndarray<float>::glob(pattern);
    if (n_timesteps == 0) n_timesteps = filenames.size();
    else filenames.resize(n_timesteps);
  }

  void advance_timestep();

protected:
  std::vector<std::string> filenames;
};

struct data_stream_float : public data_stream_file {
  void initialize();
};

struct data_stream_vti : public data_stream_file {
  void initialize() {

  }
};

struct data_stream_nc : public data_stream_file {

};

struct data_stream_pnc : public data_stream_file {

};

struct data_stream_factory {
  static data_stream* new_data_stream(int type, const std::string& source);
  static data_stream* new_data_stream(const std::string& source);
};

}

#endif
