#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <ftk/ftk_config.hh>
#include <ftk/ndarray.hh>
#include <ftk/filters/critical_point_tracker_wrapper.hh>

namespace py = pybind11;

template <typename T, int I>
static ftk::ndarray<T> to_ndarray(const py::array_t<T, I> &array)
{
  py::buffer_info buf = array.request();
  std::vector<size_t> shape;
  for (auto i = 0; i < buf.ndim; i ++)
    shape.push_back(array.shape(i));
  std::reverse(std::begin(shape), std::end(shape));

  ftk::ndarray<T> array1;
  array1.from_array((double*)buf.ptr, shape);

  return array1;
}

PYBIND11_MODULE(pyftk, m) {
    m.doc() = R"pbdoc(FTK Python bindings)pbdoc";
    
    m.def("extract_critical_points_2d_scalar", [](py::array_t<double, py::array::c_style | py::array::forcecast> array) {
      if (array.ndim() != 2)
        throw std::runtime_error("Number of dimensions must be 2");

      ftk::ndarray<double> data = to_ndarray(array);
      const size_t DW = data.dim(0), DH = data.dim(1);

      ftk::critical_point_tracker_2d_regular tracker;
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

      auto result = py::array_t<double>(critical_points.size() * 2);
      result.resize({(int)critical_points.size(), 2});
      py::buffer_info buf = result.request();

      double *ptr = (double*)buf.ptr;
      for (auto i = 0; i < critical_points.size(); i ++) {
        ptr[i*2] = critical_points[i].x[0];
        ptr[i*2+1] = critical_points[i].x[1];
      }

      return result;
    }, R"pbdoc(Extract 2D critical points)pbdoc");

    m.def("track_critical_points_2d_scalar", [](py::array_t<double, py::array::c_style | py::array::forcecast> array) {
      if (array.ndim() != 3)
        throw std::runtime_error("Number of dimensions must be 3");

      ftk::ndarray<double> data = to_ndarray(array);
      const size_t DW = data.dim(0), DH = data.dim(1), DT = data.dim(2);

      ftk::critical_point_tracker_2d_regular tracker;
      tracker.set_scalar_field_source( ftk::SOURCE_GIVEN );
      tracker.set_vector_field_source( ftk::SOURCE_DERIVED );
      tracker.set_jacobian_field_source( ftk::SOURCE_DERIVED );
      tracker.set_jacobian_symmetric( true );
      tracker.set_domain(ftk::lattice({2, 2}, {DW-3, DH-3}));
      tracker.initialize();

      for (int k = 0; k < DT; k ++) {
        const auto slice = data.slice_time(k);
        tracker.push_scalar_field_snapshot(slice);
        if (k != 0) tracker.advance_timestep();
        if (k == DT-1) tracker.update_timestep();
      }

      tracker.finalize();
      // tracker.write_traced_critical_points_text(std::cerr);

    }, R"pbdoc(Track 2D critical points)pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
