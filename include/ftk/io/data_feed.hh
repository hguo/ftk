#ifndef _FTK_REGULAR_DATA_FEED_HH
#define _FTK_REGULAR_DATA_FEED_HH

#include <ftk/ndarray.hh>

namespace ftk {

struct data_feed {
  virtual void advance_timestep() = 0;

  std::map<std::string, ndarray<double>> get_current_timestep(); // <key, array>

  std::map<std::string, std::vector<size_t>> shapes; // <key, shape>
  // std::map<std::string, 
};

struct data_feed_synthetic : public data_feed {

};

struct data_feed_io : public data_feed {
  void set_input_filename_pattern(const std::string& pattern);

  void advance_timestep();
};

struct data_feed_vti : public data_feed_io {

};

struct data_feed_nc : public data_feed_io {

};

struct data_feed_pnc : public data_feed_io {

};

}

#endif
