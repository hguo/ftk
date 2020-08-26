#ifndef _FTK_CRITICAL_POINT_TRACKER_WRAPPER_HH
#define _FTK_CRITICAL_POINT_TRACKER_WRAPPER_HH

#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/ndarray/stream.hh>
#include <ftk/ndarray/writer.hh>

namespace ftk {

struct critical_point_tracker_wrapper : public object {
  // json options:
  // - accelerator, string, 
  void configure(const json& j);

  void consume(ndarray_stream<> &stream, 
      diy::mpi::communicator comm = diy::mpi::communicator()/*MPI_COMM_WORLD*/);

  std::shared_ptr<critical_point_tracker_regular> get_tracker() {return tracker;};

  json get_json() const {return j;}

private:
  std::shared_ptr<critical_point_tracker_regular> tracker;
  json j; // config

  static bool ends_with(std::string const & value, std::string const & ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
  }
};

///////////////
void critical_point_tracker_wrapper::configure(const json& j0) 
{
  j = j0;

  /////////
  auto add_boolean_option = [&](const std::string& key, bool default_value) {
    if (j.contains(key)) {
      if (j[key].is_boolean()) {
        // OK
      } else 
        fatal("invalid " + key);
    } else 
      j[key] = default_value;
  };

  add_boolean_option("enable_streaming_trajectories", false);
  add_boolean_option("enable_discarding_interval_points", false);
  add_boolean_option("enable_discarding_degenerate_points", false);
  add_boolean_option("enable_ignoring_degenerate_points", false);

  //// 
  if (j.contains("root_proc")) {
    if (j["root_proc"].is_number()) {
      // OK
    } else
      fatal("invalid root_proc");
  } else 
    j["root_proc"] = 0; // default root proc

  if (j.contains("type_filter")) {
    if (j["type_filter"].is_string()) {
      // OK
    } else 
      fatal("invalid type_filter");
  }

  if (j.contains("output")) {
    if (j["output"].is_string()) {
      // OK
    } else 
      fatal("invalid output");
  } // else warn("empty output");

  // output type
  static const std::set<std::string> valid_output_types = {
    "discrete", // discrete and un-traced critical points.  TODO: currently ignored
    "traced", // critical point trajectories
    "sliced" // trajectories sliced into individual timesteps
  };
  static const std::string default_output_type = "traced";
  if (j.contains("output_type")) {
    if (j["output_type"].is_string()) {
      if (valid_output_types.find(j["output_type"]) != valid_output_types.end()) {
        // OK
      } else 
        fatal("invalid output_type");
    } else
      fatal("invalid output_type");
  } else 
    j["output_type"] = "traced";

  if (j["output_type"] == "sliced")
    j["enable_streaming_trajectories"] = true;

  // output format
  static const std::set<std::string> valid_output_formats = {"text", "vtp"};
  bool output_format_determined = false;
  
  if (j.contains("output_format")) {
    if (j["output_format"].is_string()) {
      if (valid_output_formats.find(j["output_format"]) != valid_output_formats.end())
        output_format_determined = true;
      else
        output_format_determined = false;
    } else 
      fatal("invalid output_format");
  } else 
    output_format_determined = false;

  if (!output_format_determined) {
    if (j.contains("output")) {
      if (ends_with(j["output"], "vtp")) j["output_format"] = "vtp";
      else j["output_format"] = "text";
    }
  }
}

void critical_point_tracker_wrapper::consume(ndarray_stream<> &stream, diy::mpi::communicator comm)
{
  if (j.is_null())
    configure(j); // make default options
  // std::cerr << j << std::endl;

  const auto js = stream.get_json();
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
 
  if (j["enable_streaming_trajectories"] == true)
    tracker->set_enable_streaming_trajectories(true);

  if (j["enable_discarding_interval_points"] == true)
    tracker->set_enable_discarding_interval_points(true);

  if (j["enable_discarding_degenerate_points"] == true)
    tracker->set_enable_discarding_degenerate_points(true);

  if (j["enable_ignoring_degenerate_points"] == true)
    tracker->set_enable_ignoring_degenerate_points(true);

  if (j.contains("type_filter")) {
    const std::string str = j["type_filter"];
    unsigned int type_filter = 0;
    if (str.find("min") != std::string::npos)
      type_filter |= ftk::CRITICAL_POINT_2D_MINIMUM;
    if (str.find("max") != std::string::npos)
      type_filter |= ftk::CRITICAL_POINT_2D_MAXIMUM;
    if (str.find("saddle") != std::string::npos)
      type_filter |= ftk::CRITICAL_POINT_2D_SADDLE;

    if (type_filter)
      tracker->set_type_filter(type_filter);
  }

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

    if (j.contains("output")) {
      if (k > 0 && j["output_type"] == "sliced") {
        const std::string pattern = j["output"];
        const std::string filename = ndarray_writer<double>::filename(pattern, k-1);
        if (j["output_format"] == "vtp")
          tracker->write_sliced_critical_points_vtk(k-1, filename);
        else 
          tracker->write_sliced_critical_points_text(k-1, filename);
      }
    }
  });

  stream.start();
  stream.finish();
  tracker->finalize();
  // delete tracker;

  if (j.contains("output")) {
    if (j["output_type"] == "traced") {
      if (j["output_format"] == "vtp") tracker->write_traced_critical_points_vtk(j["output"]);
      else tracker->write_traced_critical_points_text(j["output"]);
    }
  }
}

}

#endif
