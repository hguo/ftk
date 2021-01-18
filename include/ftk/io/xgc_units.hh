#ifndef _XGC_UNITS_HH
#define _XGC_UNITS_HH

#include <ftk/ftk_config.hh>
#include <ftk/utils/string.hh>
#include <fstream>
#include <string>

namespace ftk {

struct xgc_units_t {
  double sml_totalpe = 4096;
  double sml_dt = 1.5656919591263202E-007;
  double sml_tran = 7.8284597956316011E-005;
  double vth = 268083.9777733014;
  double ptl_num = 90000;
  double eq_axis_r = 1.724851386561000;
  double eq_axis_z = 2.0561998799429999E-002;
  double eq_axis_b = 2.105926080000000;
  double eq_tempi_v1 = 1500.000000000000;
  double ptl_ion_mass_au = 2.000000000000000;
  double ptl_ion_charge_eu = 1.000000000000000;
  double diag_1d_period = 10;
  double neu_col_period = 0;
  double psi_x = 0.2661956235889000;

public:
  xgc_units_t() {}
  xgc_units_t(const std::string& filename) { read(filename); }
  bool read(const std::string& filename);
  friend std::ostream& operator<<(std::ostream& os, const xgc_units_t&);
};

/////
inline bool xgc_units_t::read(const std::string& filename)
{
  std::ifstream ifs(filename, std::ifstream::in);
  if (ifs.is_open()) {
    std::string str;
    std::string varname, dummy;
    double value;
    while (std::getline(ifs, str)) {
      str.erase(std::remove(str.begin(), str.end(), ' '), str.end()); // remove spaces
      str.pop_back(); // remove semicolon
      auto strs = split(str, "=");
      const auto var = strs[0];
      double val = std::stod(strs[1]);

      if (var == "sml_totalpe") sml_totalpe = val;
      else if (var == "sml_dt") sml_dt = val; 
      else if (var == "sml_tran") sml_tran = val;
      else if (var == "vth") vth = val;
      else if (var == "ptl_num") ptl_num = val;
      else if (var == "eq_axis_r") eq_axis_r = val;
      else if (var == "eq_axis_z") eq_axis_z = val;
      else if (var == "eq_axis_b") eq_axis_b = val;
      else if (var == "eq_tempi_v1") eq_tempi_v1 = val;
      else if (var == "ptl_ion_mass_au") ptl_ion_mass_au = val;
      else if (var == "ptl_ion_charge_eu") ptl_ion_charge_eu = val;
      else if (var == "diag_1d_period") diag_1d_period = val;
      else if (var == "neu_col_period") neu_col_period = val;
      else if (var == "psi_x") psi_x = val;
      else return false;
      // std::cerr << var << "=" << val << std::endl;
    }
    ifs.close();
    return true;
  } else return false;
}

}

#endif
