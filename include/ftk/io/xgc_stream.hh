#ifndef _FTK_XGC_STREAM_HH
#define _FTK_XGC_STREAM_HH

#include <ftk/object.hh>
#include <ftk/ndarray.hh>

namespace ftk {
using nlohmann::json;

struct xgc_stream : public object {
  enum { XGC_FORMAT_UNKNOWN, XGC_FORMAT_H5, XGC_FORMAT_BP3, XGC_FORMAT_BP4 };

  void set_path(const std::string& p) 
  void probe_format();
  void initialize();

  std::string postfix() const { return format == XGC_FORMAT_H5 ? ".h5" : ".bp"; }
  std::string mesh_filename() const { return path + "/xgc.mesh" + postfix(); }
  std::string oneddiag_filename() const { return path + "/xgc.oneddiag" + postfix(); }
  std::string bfield_filename() const { return path + "/xgc.bfield" + postfix(); }
  std::string units_filename() const { return path + "/units.m"; }
  std::string filename(int t) const { return series_filename(path + "/xgc.3d.%05d" + postfix, t); }
  
  const xgc_units_t& get_units() const { return units; }

protected:
  std::string path;
  int format = XGC_FORMAT_UNKNOWN;

  xgc_units_t units;
};

/////

inline void xgc_stream::set_path(const std::string& p)
{
  path = p;
  if (path.length() > 0 && path.back() == '/') path.pop_back();
  if (format == XGC_FORMAT_UNKNOWN) probe_format();
}

inline void xgc_stream::probe_format() 
{
  if (file_exists(path + "/xgc.mesh.h5")) return XGC_FORMAT_H5;
  else return XGC_FORMAT_BP4;
}

inline void xgc_stream::initialize()
{

}

}


#endif
