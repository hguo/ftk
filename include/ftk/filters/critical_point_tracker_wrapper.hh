#ifndef _FTK_CRITICAL_POINT_TRACKER_WRAPPER_HH
#define _FTK_CRITICAL_POINT_TRACKER_WRAPPER_HH

#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/ndarray/stream.hh>

namespace ftk {

struct critical_point_tracker_wrapper {
  void configure(const json& j);

  void consume(ndarray_stream<> &stream);
  std::shared_ptr<critical_point_tracker_regular> get_tracker() {return tracker;};
  
private:
  std::shared_ptr<critical_point_tracker_regular> tracker;
  json j; // config
};

///////////////
void critical_point_tracker_wrapper::consume(ndarray_stream<> &stream)
{
  const auto j = stream.get_json();
  const size_t nd = j["nd"], DW = j["width"], DH = j["height"], DD = (nd == 2 ? 0 : j["depth"].get<size_t>()), DT = j["n_timesteps"];
  const size_t nv = stream.n_components();

  if (j["nd"] == 2) {
    tracker = std::shared_ptr<critical_point_tracker_regular>(new ftk::critical_point_tracker_2d_regular);
    tracker->set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
  } else {
    tracker = std::shared_ptr<critical_point_tracker_regular>(new ftk::critical_point_tracker_3d_regular);
    tracker->set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
  }
  
  // tracker->set_number_of_threads(nthreads);
  tracker->set_input_array_partial(false); // input data are not distributed

  // if (use_type_filter)
  //   tracker->set_type_filter(type_filter);
  
  if (nv == 1) { // scalar field
    tracker->set_scalar_field_source( ftk::SOURCE_GIVEN );
    tracker->set_vector_field_source( ftk::SOURCE_DERIVED );
    tracker->set_jacobian_field_source( ftk::SOURCE_DERIVED );
    if (nd == 2) { // 2D
      tracker->set_domain(ftk::lattice({0, 0}, {DW-1, DH-1})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
    } else { // 3D
      tracker->set_domain(ftk::lattice({0, 0, 0}, {DW-1, DH-1, DD-1})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
    }
  } else { // vector field
    tracker->set_scalar_field_source( ftk::SOURCE_NONE );
    tracker->set_vector_field_source( ftk::SOURCE_GIVEN );
    tracker->set_jacobian_field_source( ftk::SOURCE_DERIVED );
    if (nd == 2) { // 2D
      tracker->set_domain(ftk::lattice({1, 1}, {DW-2, DH-2})); // the indentation is needed becase the jacoobian field will be automatically derived
    } else {
      tracker->set_domain(ftk::lattice({1, 1, 1}, {DW-2, DH-2, DD-2})); // the indentation is needed becase the jacoobian field will be automatically derived
    }
  }
  tracker->initialize();
   

  auto push_timestep = [&](const ftk::ndarray<double>& field_data) {
    if (nv == 1) { // scalar field
#if 0
      if (spatial_smoothing) {
        ftk::ndarray<double> scalar = 
          ftk::conv2D_gaussian(field_data, spatial_smoothing, 
              spatial_smoothing_kernel_size, spatial_smoothing_kernel_size, 2);
        tracker->push_scalar_field_snapshot(scalar);
      } else 
#endif
        tracker->push_scalar_field_snapshot(field_data);
    }
    else // vector field
      tracker->push_vector_field_snapshot(field_data);
  };

  stream.set_callback([&](int k, ftk::ndarray<double> field_data) {
    push_timestep(field_data);
    if (k != 0) tracker->advance_timestep();
    if (k == DT-1) tracker->update_timestep();
  });

  stream.start();
  stream.finish();
  tracker->finalize();
  // delete tracker;
}

}

#endif
