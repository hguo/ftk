#ifndef _FTK_CRITICAL_POINT_TRACKER_WRAPPER_HH
#define _FTK_CRITICAL_POINT_TRACKER_WRAPPER_HH

#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/filters/critical_point_tracker_2d_unstructured.hh>
#include <ftk/filters/critical_point_tracker_3d_unstructured.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_2d_mesh.hh>
#include <ftk/ndarray/stream.hh>
#include <ftk/ndarray/writer.hh>
#include <ftk/utils/io_util.hh>

namespace ftk {

typedef std::chrono::high_resolution_clock clock_type;

struct critical_point_tracker_wrapper : public object {
  // json options:
  // - accelerator, string, 
  void configure(const json& j);

  void consume(ndarray_stream<> &stream, 
      diy::mpi::communicator comm = diy::mpi::communicator()/*MPI_COMM_WORLD*/);

  void post_process();
  void xgc_post_process();

  void write();

  std::shared_ptr<critical_point_tracker> get_tracker() {return tracker;};

  json get_json() const {return j;}

private:
  void configure_tracker_general(diy::mpi::communicator comm);
  void consume_regular(ndarray_stream<> &stream, diy::mpi::communicator comm);
  void consume_xgc(ndarray_stream<> &stream, diy::mpi::communicator comm);

  void write_sliced_results(int k);
  void write_intercepted_results(int k, int nt);

private:
  std::shared_ptr<critical_point_tracker> tracker;
  json j, js; // config

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

  auto add_number_option = [&](const std::string& key, double default_value) {
    if (j.contains(key)) {
      if (j[key].is_number()) { // OK
      } else fatal("invalid +" + key);
    } else j[key] = default_value;
  };

  auto add_string_option = [&](json &j, const std::string& key, bool required=true) {
    if (j.contains(key)) {
      if (j[key].is_string()) { // OK
      } else 
        fatal("invalid " + key);
    } else if (required) 
      fatal("missing " + key);
  };

  add_boolean_option("enable_robust_detection", true);
  add_boolean_option("enable_post_processing", true);
  add_boolean_option("enable_streaming_trajectories", false);
  add_boolean_option("enable_discarding_interval_points", false);
  add_boolean_option("enable_discarding_degenerate_points", false);
  add_boolean_option("enable_ignoring_degenerate_points", false);
  add_boolean_option("enable_timing", false);

  add_number_option("duration_pruning_threshold", 0);
  add_number_option("nblocks", 1);
  
  /// application specific
  if (j.contains("xgc")) {
    // mesh_filename
    // smoothing_kernel_filename
    // smoothing_kernel_size
    auto jx = j["xgc"];
    if (jx.is_object()) {
      add_string_option(jx, "mesh_filename");
      add_string_option(jx, "smoothing_kernel_filename");
      if (jx.contains("smoothing_kernel_size")) {
        if (jx["smoothing_kernel_size"].is_number()) { // OK
        } else 
          fatal("missing xgc smoothing_kernel_size");
      }
    } else 
      fatal("invalid xgc configuration");
  }

  /// general options
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

  add_string_option(j, "archived_discrete_critical_points_filename",  false);
  add_string_option(j, "archived_traced_critical_points_filename", false);

  add_string_option(j, "output", false);

  add_string_option(j, "accelerator", false);

  // output type
  static const std::set<std::string> valid_output_types = {
    "discrete", // discrete and un-traced critical points.  TODO: currently ignored
    "traced", // critical point trajectories
    "sliced", // trajectories sliced into individual timesteps
    "intercepted"
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

  // if (j["output_type"] == "sliced")
  //   j["enable_streaming_trajectories"] = true;

  // output format
  static const std::set<std::string> valid_output_formats = {"text", "vtp", "json"}; // , "binary"};
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
      else if (ends_with(j["output"], "txt")) j["output_format"] = "text";
      else if (ends_with(j["output"], "json")) j["output_format"] = "json";
      else j["output_format"] = "binary";
    }
  }

  if (j.contains("mesh")) { // TODO
    static const std::set<std::string> valid_mesh_file_formats = {"vtu", "hdf5", "netcdf"};
    // filename (required)
    // format (optional if format can be derived)
    // connectivity (not required for vtu)
    // coordinates (not requried for vtu)
    auto jm = j["mesh"];  
    if (jm.is_object()) {
      if (jm.contains("filename")) {
        if (jm.is_string()) {
          // OK
        } else
          fatal("missing mesh filename");
      }

      if (jm.contains("format")) {
        if (jm["format"].is_string()) {
          if (valid_mesh_file_formats.find(jm["format"]) != valid_mesh_file_formats.end()) {
            // OK
          } else fatal("mesh format not supported");
        } else fatal("invalid mesh format");
      } else { // determine mesh format
        const std::string f = jm["filename"];
        if (ends_with(f, "vtu")) jm["format"] = "vtu";
        else if (ends_with(f, "h5")) jm["format"] = "hdf5";
        else if (ends_with(f, "nc")) jm["format"] = "netcdf";
        else fatal("unable to determin mesh format");
      }

      if (jm.contains("connectivity")) {
        if (jm["connectivity"].is_string()) {
          // OK
        } else fatal("invalid connectivity");
      } else {
        if (jm["format"] == "vtu") {
          // OK
        } else fatal("missing connectivity");
      }
      
      if (jm.contains("coordinates")) {
        if (jm["coordinates"].is_string()) {
          // OK
        } else fatal("invalid coordinates");
      } else {
        if (jm["format"] == "vtu") {
          // OK
        } else fatal("missing coordinates");
      }
    } else
      fatal("invalid mesh spec");
  }
}

void critical_point_tracker_wrapper::configure_tracker_general(diy::mpi::communicator comm)
{
  tracker->set_communicator(comm);
  tracker->set_root_proc(j["root_proc"]);
 
  if (j.contains("nthreads") && j["nthreads"].is_number())
    tracker->set_number_of_threads(j["nthreads"]);

  if (j.contains("accelerator")) {
    if (j["accelerator"] == "cuda")
      tracker->use_accelerator( FTK_XL_CUDA );
    else if (j["accelerator"] == "sycl")
      tracker->use_accelerator( FTK_XL_SYCL );
    else 
      fatal(FTK_ERR_ACCELERATOR_UNSUPPORTED);
  }

  if (j.contains("thread_backend")) {
    std::string backend = j["thread_backend"];
    tracker->use_thread_backend( backend );
  }

  if (j.contains("nblocks"))
    tracker->set_number_of_blocks(j["nblocks"]);

  tracker->set_input_array_partial(false); // input data are not distributed

  // if (use_type_filter)
  //   tracker->set_type_filter(type_filter);

  if (j.contains("enable_robust_detection"))
    tracker->set_enable_robust_detection( j["enable_robust_detection"].get<bool>() );

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
}

void critical_point_tracker_wrapper::consume(ndarray_stream<> &stream, diy::mpi::communicator comm)
{
  if (j.is_null())
    configure(j); // make default options
  // std::cerr << j << std::endl;

  if (j.contains("xgc"))
    consume_xgc(stream, comm);
  else 
    consume_regular(stream, comm);
}

void critical_point_tracker_wrapper::consume_xgc(ndarray_stream<> &stream, diy::mpi::communicator comm)
{
  // const auto js = stream.get_json();
  js = stream.get_json();
  const size_t DT = js["n_timesteps"];

  const std::string mesh_filename = j["xgc"]["mesh_filename"];
  const std::string smoothing_kernel_filename = j["xgc"]["smoothing_kernel_filename"];
  const double smoothing_kernel_size = j["xgc"]["smoothing_kernel_size"];

  // load mesh & data from hdf5
  ftk::ndarray<int> triangles;
  ftk::ndarray<double> coords, psi;

  triangles.read_h5(mesh_filename, "/cell_set[0]/node_connect_list");
  coords.read_h5(mesh_filename, "/coordinates/values");
  psi.read_h5(mesh_filename, "psi");
  
  // build mesh
  ftk::simplicial_unstructured_2d_mesh<> m(coords, triangles);
  m.build_edges();

  // tracker set up
  tracker = std::shared_ptr<critical_point_tracker_2d_unstructured>(
      new critical_point_tracker_2d_unstructured(comm, m));
  configure_tracker_general(comm);

  tracker->set_scalar_components({"dneOverne0", "psi"});

  if (j.contains("archived_traced_critical_points_filename")) {
    fprintf(stderr, "reading archived traced critical points...\n");
    const std::string filename = j["archived_traced_critical_points_filename"];
    if (ends_with(filename, "json")) tracker->read_traced_critical_points_json(filename);
    else tracker->read_traced_critical_points_binary(filename);
    // fprintf(stderr, "done reading.\n");
    return;
  } else if (j.contains("archived_discrete_critical_points_filename")) {
    const std::string filename = j["archived_discrete_critical_points_filename"];
    if (ends_with(filename, "json")) tracker->read_critical_points_json(filename);
    else tracker->read_critical_points_binary(filename);
    tracker->finalize();
    return;
  }

  // load smoothing kernel
  bool succ = m.read_smoothing_kernel(smoothing_kernel_filename);
  if (!succ) {
    m.build_smoothing_kernel(smoothing_kernel_size);
    m.write_smoothing_kernel(smoothing_kernel_filename);
  }
  // fprintf(stderr, "mesh loaded., %zu, %zu, %zu\n", m.n(0), m.n(1), m.n(2));

  auto push_timestep = [&](int k, const ftk::ndarray<double>& data) {
    auto dpot = data.get_transpose();
    dpot.reshape(dpot.dim(0));

    ftk::ndarray<double> scalar, grad, J;
    m.smooth_scalar_gradient_jacobian(dpot, /*smoothing_kernel_size,*/ scalar, grad, J);
  
    ftk::ndarray<double> scalars = ftk::ndarray<double>::concat({scalar, psi});

    tracker->push_field_data_snapshot(scalars, grad, J);

    if (j["xgc"].contains("write_back_filename")) { // write data back to vtu files
      const std::string pattern = j["xgc"]["write_back_filename"];
      const std::string filename = series_filename(pattern, k);
      // m.scalar_to_vtk_unstructured_grid_data_file(filename, "dneOverne0", dpot);
      m.scalar_to_vtk_unstructured_grid_data_file(filename, "dneOverne0", scalar);
    }
  };
  
  stream.set_callback([&](int k, const ftk::ndarray<double> &field_data) {
    if (j["xgc"].contains("torus") && j["xgc"]["torus"] == true) { // tracking over torus
      auto dpot = field_data.get_transpose();
      for (int k = 0; k < dpot.dim(1); k ++) {
        ftk::ndarray<double> dpot_slice = dpot.slice_time(k), scalar, grad, J;
        m.smooth_scalar_gradient_jacobian(dpot_slice, /*smoothing_kernel_size,*/ scalar, grad, J);

        ftk::ndarray<double> scalars = ftk::ndarray<double>::concat({scalar, psi});
        tracker->push_field_data_snapshot(scalars, grad, J);

        if (k != 0) tracker->advance_timestep();
        if (k == dpot.dim(1)-1) tracker->update_timestep();
      
        if (j["xgc"].contains("write_back_filename")) { // write data back to vtu files
          const std::string pattern = j["xgc"]["write_back_filename"];
          const std::string filename = series_filename(pattern, k);
          // m.scalar_to_vtk_unstructured_grid_data_file(filename, "dneOverne0", dpot);
          m.scalar_to_vtk_unstructured_grid_data_file(filename, "dneOverne0", scalar);
        }
      }
    } else { // tracking over time
      push_timestep(k, field_data);
      if (k != 0) tracker->advance_timestep();
      if (k == DT-1) tracker->update_timestep();

      if (k>0 && j.contains("output") && j["output_type"] == "sliced" && j["enable_streaming_trajectories"] == true)
        write_sliced_results(k-1);
    }
  });

  stream.start();
  stream.finish();
  tracker->finalize();
}

void critical_point_tracker_wrapper::write_sliced_results(int k)
{
  const std::string pattern = j["output"];
  const std::string filename = series_filename(pattern, k);
  if (j["output_format"] == "vtp")
    tracker->write_sliced_critical_points_vtk(k, filename);
  else 
    tracker->write_sliced_critical_points_text(k, filename);
}

void critical_point_tracker_wrapper::write_intercepted_results(int k, int nt)
{
  const std::string pattern = j["output"];
  const std::string filename = series_filename(pattern, k);
  if (j["output_format"] == "vtp") {
    int nt = 2;
    if (j.contains("intercept_length") && j["intercept_length"].is_number())
      nt = j["intercept_length"];
    tracker->write_intercepted_critical_points_vtk(k-nt, k, filename);
  }
}

void critical_point_tracker_wrapper::consume_regular(ndarray_stream<> &stream, diy::mpi::communicator comm)
{
  auto t0 = clock_type::now();

  const auto js = stream.get_json();
  const size_t nd = stream.n_dimensions(),
               DW = js["dimensions"][0], 
               DH = js["dimensions"].size() > 1 ? js["dimensions"][1].get<int>() : 0,
               DD = js["dimensions"].size() > 2 ? js["dimensions"][2].get<int>() : 0,
               DT = js["n_timesteps"];
  const size_t nv = stream.n_components();

  std::shared_ptr<critical_point_tracker_regular> rtracker;
  if (nd == 2) {
    rtracker.reset(new critical_point_tracker_2d_regular(comm));
    rtracker->set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
  } else {
    rtracker.reset(new critical_point_tracker_3d_regular(comm));
    rtracker->set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
  }

  if (nv == 1) { // scalar field
    rtracker->set_scalar_field_source( ftk::SOURCE_GIVEN );
    rtracker->set_vector_field_source( ftk::SOURCE_DERIVED );
    rtracker->set_jacobian_field_source( ftk::SOURCE_DERIVED );
    rtracker->set_jacobian_symmetric( true );
    if (nd == 2) { // 2D
      rtracker->set_domain(ftk::lattice({2, 2}, {DW-3, DH-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
    } else { // 3D
      rtracker->set_domain(ftk::lattice({2, 2, 2}, {DW-3, DH-3, DD-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
    }
    fprintf(stderr, "treating input data as scalar field.\n");
  } else { // vector field
    rtracker->set_scalar_field_source( ftk::SOURCE_NONE );
    rtracker->set_vector_field_source( ftk::SOURCE_GIVEN );
    rtracker->set_jacobian_field_source( ftk::SOURCE_DERIVED );
    rtracker->set_jacobian_symmetric( false );
    if (nd == 2) { // 2D
      rtracker->set_domain(ftk::lattice({1, 1}, {DW-2, DH-2})); // the indentation is needed becase the jacoobian field will be automatically derived
    } else {
      rtracker->set_domain(ftk::lattice({1, 1, 1}, {DW-2, DH-2, DD-2})); // the indentation is needed becase the jacoobian field will be automatically derived
    }
    fprintf(stderr, "treating input data as vector field.\n");
  }
  tracker = rtracker;
  
  configure_tracker_general(comm);
  tracker->initialize();
 
  if (j.contains("archived_traced_critical_points_filename")) {
    fprintf(stderr, "reading archived traced critical points...\n");
    const std::string filename = j["archived_traced_critical_points_filename"];
    if (ends_with(filename, "json")) tracker->read_traced_critical_points_json(filename);
    else tracker->read_traced_critical_points_binary(filename);
    // fprintf(stderr, "done reading.\n");
    return;
  } else if (j.contains("archived_discrete_critical_points_filename")) {
    const std::string filename = j["archived_discrete_critical_points_filename"];
    if (ends_with(filename, "json")) tracker->read_critical_points_json(filename);
    else tracker->read_critical_points_binary(filename);
    tracker->finalize();
    return;
  }

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
    
    if (k>0 && j.contains("output") && j["output_type"] == "sliced" && j["enable_streaming_trajectories"] == true)
      write_sliced_results(k-1);
  });

  auto t1 = clock_type::now();

  stream.start();
  stream.finish();

  auto t2 = clock_type::now();
  
  tracker->finalize();
  auto t3 = clock_type::now();

  if (comm.rank() == 0 && j["enable_timing"]) {
    float t_init = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9,
          t_compute = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() * 1e-9,
          t_finalize = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() * 1e-9;
    fprintf(stderr, "t_init=%f, t_compute=%f, t_finalize=%f\n", t_init, t_compute, t_finalize);
  }
  // delete tracker;
}

void critical_point_tracker_wrapper::xgc_post_process()
{
  fprintf(stderr, "post processing for xgc...\n");

  tracker->select_trajectories([](const ftk::feature_curve_t& traj) {
    // if (traj.tmax - traj.tmin /*duration*/< 2.0) return false;
    
    if (traj.consistent_type == ftk::CRITICAL_POINT_2D_MAXIMUM && traj.max[0] > 5.0) return true;
    else if (traj.consistent_type == ftk::CRITICAL_POINT_2D_MINIMUM && traj.min[0] < -5.0) return true;

    // if (traj.max[1] /*max of psi*/ < 0.2) return false;
    // if (traj.min[1] /*min of psi*/ > 0.3) return false;
    //// if (traj.max[0] /*max of dpot*/ < 0.0) return false;
    return false;
  });
  
  auto &trajs = tracker->get_traced_critical_points();
  trajs.foreach([](ftk::feature_curve_t& t) {
    t.discard_interval_points();
    t.derive_velocity();
  });
  
  tracker->slice_traced_critical_points();
  tracker->select_sliced_critical_points([](const ftk::feature_point_t& cp) {
    // if (cp.scalar[1] < 0.18 || cp.scalar[1] > 0.3) return false;
    if (cp.type == ftk::CRITICAL_POINT_2D_MAXIMUM && cp.scalar[0] > 5.0) return true;
    else if (cp.type == ftk::CRITICAL_POINT_2D_MINIMUM && cp.scalar[0] < -5.0) return true;
    else return false;
  });
}

void critical_point_tracker_wrapper::post_process()
{
  auto &trajs = tracker->get_traced_critical_points();

  trajs.foreach([](ftk::feature_curve_t& t) {
    t.discard_high_cond();
    t.smooth_ordinal_types();
    t.smooth_interval_types();
    t.rotate();
    // t.discard_interval_points();
    t.update_statistics();
  });
  if (j["duration_pruning_threshold"] > 0) {
    const double threshold = j["duration_pruning_threshold"];
    trajs.filter([&](const feature_curve_t& traj) {
      if (traj.tmax - traj.tmin < threshold) return false;
      else return true;
    });
  }
  trajs.split_all();
  if (j["enable_discarding_interval_points"] == true) {
    trajs.foreach([](ftk::feature_curve_t& t) {
      t.discard_interval_points();
    });
  }
  trajs.foreach([](ftk::feature_curve_t& t) {
    t.reorder();
    t.adjust_time();
    // t.relabel(k);
    t.update_statistics();
  });
  
  if (j.contains("xgc") && j["xgc"].contains("post_process") && j["xgc"]["post_process"] == true)
    xgc_post_process();
  
  if (j.contains("enable_deriving_velocities")) {
    trajs.foreach([](ftk::feature_curve_t& t) {
      t.discard_interval_points();
      t.derive_velocity();
      t.update_statistics();
    });
  }
}

void critical_point_tracker_wrapper::write()
{
  if (j.contains("output")) {
    if (j["output_type"] == "sliced") {
      fprintf(stderr, "slicing and writing..\n");
      if (tracker->get_sliced_critical_points().empty())
        tracker->slice_traced_critical_points();
      for (const auto &kv : tracker->get_sliced_critical_points()) 
        write_sliced_results(kv.first);
    } else if (j["output_type"] == "intercepted") {
      fprintf(stderr, "writing intercepted critical points..\n");
      // for (int t = 0; t < tracker->get_current_timestep(); t ++)
      for (int t = 0; t < js["n_timesteps"]; t ++)
        write_intercepted_results(t, 2);
    } else if (j["output_type"] == "traced") {
      fprintf(stderr, "writing traced critical points..\n");
      if (j["output_format"] == "vtp") tracker->write_traced_critical_points_vtk(j["output"]);
      else if (j["output_format"] == "text") tracker->write_traced_critical_points_text(j["output"]);
      else if (j["output_format"] == "json") tracker->write_traced_critical_points_json(j["output"]);
      else tracker->write_traced_critical_points_binary(j["output"]);
    } else if (j["output_type"] == "discrete") {
      fprintf(stderr, "writing discrete critical points..\n");
      if (j["output_format"] == "vtp") tracker->write_critical_points_vtk(j["output"]);
      else if (j["output_format"] == "text") tracker->write_critical_points_text(j["output"]);
      else if (j["output_format"] == "json") tracker->write_critical_points_json(j["output"]);
      else tracker->write_critical_points_binary(j["output"]);
    }
  }
}

}

#endif
