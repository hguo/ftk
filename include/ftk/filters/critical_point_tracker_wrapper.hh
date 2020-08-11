#ifndef _FTK_CRITICAL_POINT_TRACKER_WRAPPER_HH
#define _FTK_CRITICAL_POINT_TRACKER_WRAPPER_HH

#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/ndarray/stream.hh>

namespace ftk {

struct critical_point_tracker_wrapper : public object {
  // json options:
  // - accelerator, string, 
  void configure(const json& j);

  void consume(ndarray_stream<> &stream, 
      diy::mpi::communicator comm = diy::mpi::communicator()/*MPI_COMM_WORLD*/);

  std::shared_ptr<critical_point_tracker_regular> get_tracker() {return tracker;};
  
private:
  std::shared_ptr<critical_point_tracker_regular> tracker;
  json j; // config
};

///////////////
void critical_point_tracker_wrapper::configure(const json& j0) 
{
  j = j0;
  if (j.contains("root_proc")) {
    if (j["root_proc"].is_number()) {
      // OK
    } else
      fatal("invalid root_proc");
  } else 
    j["root_proc"] = 0; // default root proc
}

void critical_point_tracker_wrapper::consume(ndarray_stream<> &stream, diy::mpi::communicator comm)
{
  if (j.is_null())
    configure(j); // make default options

  const auto js = stream.get_json();
  // std::cerr << j << std::endl;
  const size_t nd = stream.n_dimensions(),
               DW = js["dimensions"][0], 
               DH = js["dimensions"][1],
               DT = js["n_timesteps"];
  size_t DD;
  if (nd == 3) DD = js["dimensions"][2];
  else DD = 0;
  const size_t nv = stream.n_components();

  if (nd == 2) {
    tracker = std::shared_ptr<critical_point_tracker_regular>(new ftk::critical_point_tracker_2d_regular);
    tracker->set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
  } else {
    tracker = std::shared_ptr<critical_point_tracker_regular>(new ftk::critical_point_tracker_3d_regular);
    tracker->set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
  }

  tracker->set_root_proc(j["root_proc"]);
  
  // tracker->set_number_of_threads(nthreads);
  tracker->set_input_array_partial(false); // input data are not distributed
  tracker->set_communicator(comm);

  // if (use_type_filter)
  //   tracker->set_type_filter(type_filter);
  
  if (nv == 1) { // scalar field
    tracker->set_scalar_field_source( ftk::SOURCE_GIVEN );
    tracker->set_vector_field_source( ftk::SOURCE_DERIVED );
    tracker->set_jacobian_field_source( ftk::SOURCE_DERIVED );
    tracker->set_jacobian_symmetric( true );
    if (nd == 2) { // 2D
      tracker->set_domain(ftk::lattice({2, 2}, {DW-3, DH-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
    } else { // 3D
      tracker->set_domain(ftk::lattice({2, 2, 2}, {DW-3, DH-3, DD-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
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

  stream.set_callback([&](int k, const ftk::ndarray<double> &field_data) {
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
