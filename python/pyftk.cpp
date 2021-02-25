#include <ftk/config.hh>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <ftk/ndarray.hh>
#include <ftk/filters/json_interface.hh>

diy::mpi::environment env;
diy::mpi::communicator comm;

namespace py = pybind11;

PYBIND11_MODULE(pyftk, m) {
  m.doc() = R"pbdoc(FTK Python bindings)pbdoc";
 
  py::module extractors = m.def_submodule("extractors", "Feature extraction algorithms");
  extractors.def("extract_critical_points_2d_scalar", [](py::array_t<double, py::array::c_style | py::array::forcecast> array) {
    if (array.ndim() != 4)
      throw std::runtime_error("Number of dimensions must be 4 (1, width, height, 1)");

    ftk::ndarray<double> data(array);
    const size_t DW = data.dim(0), DH = data.dim(1);
    data.reshape(DW, DH);

    ftk::critical_point_tracker_2d_regular tracker(comm);
    tracker.set_scalar_field_source( ftk::SOURCE_GIVEN );
    tracker.set_vector_field_source( ftk::SOURCE_DERIVED );
    tracker.set_jacobian_field_source( ftk::SOURCE_DERIVED );
    tracker.set_jacobian_symmetric( true );
    tracker.set_domain(ftk::lattice({2, 2}, {DW-3, DH-3}));
    tracker.initialize();

    tracker.push_scalar_field_snapshot(data);
    tracker.update_timestep();
    // tracker.finalize();
    // tracker.write_critical_points_text(std::cerr);

    const auto critical_points = tracker.get_critical_points();
    py::list results;

    for (const auto &cp : critical_points) {
      py::dict pcp;
      pcp["x"] = cp.x[0];
      pcp["y"] = cp.x[1];
      pcp["t"] = cp.t;
      pcp["type"] = ftk::critical_point_type_to_string(2, cp.type, true);
      pcp["scalar"] = cp.scalar[0];
      results.append(pcp);
    }

    return results;
  }, R"pbdoc(Extract critical points in a 2D scalar field)pbdoc");
  
  extractors.def("extract_critical_points_2d_vector", [](py::array_t<double, py::array::c_style | py::array::forcecast> array) {
    if (array.ndim() != 4)
      throw std::runtime_error("Number of dimensions must be 4 (2, width, height, 1)");

    ftk::ndarray<double> data(array);
    if (data.dim(0) != 2) 
      throw std::runtime_error("The first dimension must be 2");
    const size_t DW = data.dim(1), DH = data.dim(2);
    data.reshape(2, DW, DH);

    ftk::critical_point_tracker_2d_regular tracker(comm);
    tracker.set_scalar_field_source( ftk::SOURCE_NONE );
    tracker.set_vector_field_source( ftk::SOURCE_DERIVED );
    tracker.set_jacobian_field_source( ftk::SOURCE_DERIVED );
    tracker.set_jacobian_symmetric( false );
    tracker.set_input_array_partial( false );
    tracker.set_domain(ftk::lattice({1, 1}, {DW-2, DH-2}));
    tracker.initialize();

    tracker.push_vector_field_snapshot(data);
    tracker.update_timestep();
    // tracker.finalize();
    // tracker.write_critical_points_text(std::cerr);

    const auto critical_points = tracker.get_critical_points();
    py::list results;

    for (const auto &cp : critical_points) {
      py::dict pcp;
      pcp["x"] = cp.x[0];
      pcp["y"] = cp.x[1];
      pcp["t"] = cp.t;
      pcp["type"] = ftk::critical_point_type_to_string(2, cp.type, false);
      results.append(pcp);
    }

    return results;
  }, R"pbdoc(Extract critical points in a vector field)pbdoc");

  py::module trackers = m.def_submodule("trackers", "Feature tracking algorithms");
  trackers.def("track_critical_points_2d_scalar", [](py::array_t<double, py::array::c_style | py::array::forcecast> array) {
    if (array.ndim() != 4)
      throw std::runtime_error("Number of dimensions must be 4: (1, width, height, time)");

    ftk::ndarray<double> data(array);
    const size_t DW = data.dim(1), DH = data.dim(2), DT = data.dim(3);
    data.reshape(DW, DH, DT);

    ftk::critical_point_tracker_2d_regular tracker(comm);
    tracker.set_scalar_field_source( ftk::SOURCE_GIVEN );
    tracker.set_vector_field_source( ftk::SOURCE_DERIVED );
    tracker.set_jacobian_field_source( ftk::SOURCE_DERIVED );
    tracker.set_jacobian_symmetric( true );
    tracker.set_domain(ftk::lattice({2, 2}, {DW-3, DH-3}));
    tracker.set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
    tracker.initialize();

    for (int k = 0; k < DT; k ++) {
      const auto slice = data.slice_time(k);
      tracker.push_scalar_field_snapshot(slice);
      if (k != 0) tracker.advance_timestep();
      if (k == DT-1) tracker.update_timestep();
    }

    tracker.finalize();
    // tracker.write_traced_critical_points_text(std::cerr);

    const auto traces = tracker.get_traced_critical_points();
    py::list result;

    for (const auto &kv : traces) {
      py::dict traj;
      py::list trace;
      for (auto j = 0; j < kv.second.size(); j ++) {
        const auto &cp = kv.second[j];
        py::dict pcp;
        pcp["x"] = cp.x[0];
        pcp["y"] = cp.x[1];
        pcp["t"] = cp.t;
        pcp["type"] = ftk::critical_point_type_to_string(2, cp.type, true);
        pcp["scalar"] = cp.scalar[0];
        trace.append(pcp);
      }
      traj["length"] = kv.second.size();
      traj["trace"] = trace;
      result.append(traj);
    }

    return result;
  }, R"pbdoc(Track 2D critical points)pbdoc");

  py::module synth = m.def_submodule("synthesizers", "Synthetic data generator");
  synth.def("spiral_woven", [](int DW, int DH, int DT) {
    auto array = ftk::synthetic_woven_2Dt<double>(DW, DH, DT);
    array.reshape(1, DW, DH, DT);
    return array.to_numpy();
  }, R"pbdoc(Generate 2D spiral woven data)pbdoc");

  synth.def("double_gyre_flow", [](int DW, int DH, int DT) {
    if (DT == 1) {
      auto array = ftk::synthetic_double_gyre(DW, DH, 0.0);
      array.reshape(2, DW, DH, 1);
      return array.to_numpy();
    } else if (DT > 1) {
      std::vector<ftk::ndarray<double>> arrays;
      for (int i = 0; i < DT; i ++)
        arrays.push_back( ftk::synthetic_double_gyre(DW, DH, (double)i*0.1) );
      auto array = ftk::ndarray<double>::stack(arrays);
      array.reshape(2, DW, DH, DT);
      return array.to_numpy();
    } else {
      throw std::runtime_error("DT must be an integer greater than 1");
    }
  }, R"pbdoc(Generate 2D double gyre flow)pbdoc");

  synth.def("moving_extremum", [](int DW, int DH, int DT, double x0, double y0, double dir_x, double dir_y) {
    // fprintf(stderr, "%f, %f, %f,  %f\n", x0, y0, dir_x, dir_y);
    const double xc[2] = {x0, y0}, dir[2] = {dir_x, dir_y};
    std::vector<ftk::ndarray<double>> arrays;
    std::vector<size_t> shape = {(size_t)DW, (size_t)DH};
    for (int k = 0; k < DT; k ++) {
      auto a = ftk::synthetic_moving_extremum<double, 2>(shape, xc, dir, (double)k);
      a.reshape(1, DW, DH);
      arrays.push_back( a );
    }
    return ftk::ndarray<double>::stack(arrays).to_numpy();
  }, R"pbdoc(Generate moving extremum data)pbdoc");

  py::module numeric = m.def_submodule("numeric", "Numeric functions");
  numeric.def("det4", [](py::array_t<double, py::array::c_style | py::array::forcecast> array) {
    if (array.ndim() != 2)
      throw std::runtime_error("Number of dimensions must be 2");

    ftk::ndarray<double> m(array);
    double M[4][4];
    for (int i = 0; i < 4; i ++)
      for (int j = 0; j < 4; j ++)
        M[i][j] = m(j, i);

    return ftk::det4(M);
  });
  
  numeric.def("det3", [](py::array_t<double, py::array::c_style | py::array::forcecast> array) {
    if (array.ndim() != 2)
      throw std::runtime_error("Number of dimensions must be 2");

    ftk::ndarray<double> m(array);
    double M[3][3];
    for (int i = 0; i < 3; i ++)
      for (int j = 0; j < 3; j ++)
        M[i][j] = m(j, i);

    return ftk::det3(M);
  });

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
