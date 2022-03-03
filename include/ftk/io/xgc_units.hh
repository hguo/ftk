#ifndef _XGC_UNITS_HH
#define _XGC_UNITS_HH

#include <ftk/config.hh>
#include <ftk/utils/string.hh>
#include <fstream>
#include <string>
#include <ftk/external/json.hh>

namespace ftk {

struct xgc_units_t {
  int sml_totalpe = 4096;
  double sml_dt = 1.5656919591263202E-007;
  double sml_tran = 7.8284597956316011E-005;
  double vth = 268083.9777733014;
  int ptl_num = 90000;
  double eq_axis_r = 1.724851386561000;
  double eq_axis_z = 2.0561998799429999E-002;
  double eq_axis_b = 2.105926080000000;
  double eq_tempi_v1 = 1500.000000000000;
  double ptl_ion_mass_au = 2.000000000000000;
  double ptl_ion_charge_eu = 1.000000000000000;
  int diag_1d_period = 10;
  int neu_col_period = 0;
  double psi_x = 0.2661956235889000;

  // double eq_x_psi = 0.0697345; // psi_x
  double eq_x_r = 2.8;
  double eq_x_z = -0.99988;
  double eq_den_v1 = 3.5e+19;
  double ptl_charge_eu = 1;
  int sml_wedge_n = 1;
  double ptl_mass_au = 1;

public:
  xgc_units_t() {}
  xgc_units_t(const std::string& filename) { read(filename); }
  bool read(const std::string& filename);
  bool read_bp(const std::string& filename, diy::mpi::communicator comm = MPI_COMM_WORLD);
};

}

namespace nlohmann
{
  using namespace ftk;
  template <>
  struct adl_serializer<xgc_units_t> {
    static void to_json(json& j, const xgc_units_t& s) {
      j = {
        {"sml_totalpe", s.sml_totalpe},
        {"sml_dt", s.sml_dt},
        {"sml_tran", s.sml_tran},
        {"vth", s.vth},
        {"ptl_num", s.ptl_num},
        {"eq_axis_r", s.eq_axis_r},
        {"eq_axis_z", s.eq_axis_z},
        {"eq_axis_b", s.eq_axis_b},
        {"eq_tempi_v1", s.eq_tempi_v1},
        {"ptl_ion_mass_au", s.ptl_ion_mass_au},
        {"ptl_ion_charge_eu", s.ptl_ion_charge_eu},
        {"diag_1d_period", s.diag_1d_period},
        {"neu_col_period", s.neu_col_period},
        {"psi_x", s.psi_x},
        // {"eq_x_psi", s.eq_x_psi},
        {"eq_x_r", s.eq_x_r},
        {"eq_x_z", s.eq_x_z},
        {"eq_den_v1", s.eq_den_v1},
        {"ptl_charge_eu", s.ptl_charge_eu},
        {"sml_wedge_n", s.sml_wedge_n},
        {"ptl_mass_au", s.ptl_mass_au}
      };
    }
   
    static void from_json(const json&j, xgc_units_t& s) {
      // TODO
    }
  };
}


/////
namespace ftk {
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
      else if (var == "eq_x_r") eq_x_r = val;
      else if (var == "eq_x_z") eq_x_z = val;
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
  } else 
    return false; 
}

bool xgc_units_t::read_bp(const std::string& filename, diy::mpi::communicator comm)
{
#if FTK_HAVE_ADIOS2
#if ADIOS2_USE_MPI
  adios2::ADIOS adios(comm);
#else
  adios2::ADIOS adios;
#endif
  adios2::IO io = adios.DeclareIO("BPReader");
  adios2::Engine reader = io.Open(filename, adios2::Mode::Read);
  if (!reader) return false;

  auto var_diag_1d_period = io.InquireVariable<int>("diag_1d_period");
  reader.Get<int>(var_diag_1d_period, &diag_1d_period);

  auto var_eq_axis_b = io.InquireVariable<double>("eq_axis_b");
  reader.Get<double>(var_eq_axis_b, &eq_axis_b);

  auto var_eq_axis_r = io.InquireVariable<double>("eq_axis_r");
  reader.Get<double>(var_eq_axis_r, &eq_axis_r);

  auto var_eq_axis_z = io.InquireVariable<double>("eq_axis_z");
  reader.Get<double>(var_eq_axis_z, &eq_axis_z);

  auto var_eq_den_v1 = io.InquireVariable<double>("eq_den_v1");
  reader.Get<double>(var_eq_den_v1, &eq_den_v1);

  auto var_eq_tempi_v1 = io.InquireVariable<double>("eq_tempi_v1");
  reader.Get<double>(var_eq_tempi_v1, &eq_tempi_v1);

  auto var_eq_x_psi = io.InquireVariable<double>("eq_x_psi");
  reader.Get<double>(var_eq_x_psi, &psi_x);
  fprintf(stderr, "psi_x=%f\n", psi_x);

  auto var_eq_x_r = io.InquireVariable<double>("eq_x_r");
  reader.Get<double>(var_eq_x_r, &eq_x_r);

  auto var_eq_x_z = io.InquireVariable<double>("eq_x_z");
  reader.Get<double>(var_eq_x_z, &eq_x_z);

  auto var_neu_col_period = io.InquireVariable<int>("neu_col_period");
  reader.Get<int>(var_neu_col_period, &neu_col_period);

  auto var_ptl_charge_eu = io.InquireVariable<double>("ptl_charge_eu");
  reader.Get<double>(var_ptl_charge_eu, &ptl_charge_eu);

  auto var_ptl_mass_au = io.InquireVariable<double>("ptl_mass_au");
  reader.Get<double>(var_ptl_mass_au, &ptl_mass_au);

  auto var_ptl_num = io.InquireVariable<int>("ptl_num");
  reader.Get<int>(var_ptl_num, &ptl_num);

  auto var_sml_dt = io.InquireVariable<double>("sml_dt");
  reader.Get<double>(var_sml_dt, &sml_dt);

  auto var_sml_totalpe = io.InquireVariable<int>("sml_totalpe");
  reader.Get<int>(var_sml_totalpe, &sml_totalpe);

  auto var_sml_tran = io.InquireVariable<double>("sml_tran");
  reader.Get<double>(var_sml_tran, &sml_tran);

  auto var_sml_wedge_n = io.InquireVariable<int>("sml_wedge_n");
  reader.Get<int>(var_sml_wedge_n, &sml_wedge_n);

  reader.Close();
  return true;
#else
  return false;
#endif
}

}

#endif
