#include <ftk/filters/critical_point_tracker_2d_unstructured.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_2d_mesh.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/external/cxxopts.hpp>
#include <vtkXMLPolyDataWriter.h>
#include "../../cli/constants.hh"

std::string input_filename_pattern,
  mesh_filename, 
  kernel_filename = "xgc.kernel",
  output_filename;
std::vector<std::string> input_filenames;
bool enable_streaming_trajectories = false;
double sigma(0.02);

void parse_arguments(int argc, char **argv)
{
  const int argc0 = argc; 
  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "Input file name pattern: a single file or a series of file, e.g. 'scalar.raw', 'cm1out_000*.nc'",
     cxxopts::value<std::string>(input_filename_pattern))
    ("o,output", "Output file", cxxopts::value<std::string>(output_filename))
    ("m,mesh", "Input mesh file", cxxopts::value<std::string>(mesh_filename))
    ("k,kernel", "Input/output smoothing kernel file", cxxopts::value<std::string>(kernel_filename))
    ("stream", "Streaming trajectories", cxxopts::value<bool>(enable_streaming_trajectories))
    ("s,sigma", "Kernel bandwidth", cxxopts::value<double>(sigma));
  auto results = options.parse(argc, argv);

  if (argc0 < 2) {
    std::cerr << options.help() << std::endl;
    exit(1);
  }

  if (!output_filename.length()) {
    fprintf(stderr, "missing output filename.\n");
    exit(1);
  }

  if (!mesh_filename.length()) {
    fprintf(stderr, "missing mesh filename.\n");
    exit(1);
  }
}

int main(int argc, char **argv)
{
  parse_arguments(argc, argv);

  // load mesh & data from hdf5
  ftk::ndarray<int> triangles;
  ftk::ndarray<double> coords, psi;
  
  triangles.from_h5(mesh_filename, "/cell_set[0]/node_connect_list");
  coords.from_h5(mesh_filename, "/coordinates/values");
  psi.from_h5(mesh_filename, "psi");

  // build mesh
  ftk::simplicial_unstructured_2d_mesh<> m(coords, triangles);
  m.build_edges();

  if (kernel_filename.length()) {
    bool succ = m.read_smoothing_kernel(kernel_filename);
    if (!succ) {
      m.build_smoothing_kernel(sigma);
      m.write_smoothing_kernel(kernel_filename);
    }
  }

  fprintf(stderr, "mesh loaded., %zu, %zu, %zu\n", m.n(0), m.n(1), m.n(2));

  input_filenames = ftk::ndarray<double>::glob(input_filename_pattern);
 
  ftk::critical_point_tracker_2d_unstructured tracker(m);
  // tracker.set_type_filter( ftk::CRITICAL_POINT_2D_MAXIMUM );

  if (enable_streaming_trajectories)
    tracker.set_enable_streaming_trajectories(true);

  tracker.set_scalar_components({"dneOverne0", "psi"}); // dpot and psi

  if (input_filenames.size() > 1) { // track over time
    for (int t = 0; t < input_filenames.size(); t ++) {
      ftk::ndarray<double> dpot;
      dpot.from_h5(input_filenames[t], "/dneOverne0");
      dpot = dpot.transpose();
      dpot.reshape(dpot.dim(0)); // only use the first slice

      ftk::ndarray<double> scalar, grad, J;
      m.smooth_scalar_gradient_jacobian(dpot, sigma, scalar, grad, J);
  
      ftk::ndarray<double> scalars = ftk::ndarray<double>::concat({scalar, psi});
      tracker.push_field_data_snapshot(scalars, grad, J);

      if (t != 0) tracker.advance_timestep();
      if (t == input_filenames.size()-1) tracker.update_timestep();
    }
  } else { // track over poloidal planes
    ftk::ndarray<double> dpot;
    dpot.from_h5(input_filenames[0], "/dneOverne0");
    dpot = dpot.transpose();

    for (int k = 0; k < dpot.dim(1); k ++) {
      ftk::ndarray<double> dpot_slice = dpot.slice_time(k), scalar, grad, J;
      m.smooth_scalar_gradient_jacobian(dpot_slice, sigma, scalar, grad, J);
   
      ftk::ndarray<double> scalars = ftk::ndarray<double>::concat({scalar, psi});
      tracker.push_field_data_snapshot(scalars, grad, J);

      if (k != 0) tracker.advance_timestep();
      if (k == dpot.dim(1)-1) tracker.update_timestep();
    }
  }
  tracker.finalize();

  tracker.select_traj([](const ftk::critical_point_traj_t& traj) {
    if (traj.max[1] /*max of psi*/ < 0.2) return false;
    // if (traj.min[1] /*min of psi*/ > 0.27) return false;
    // if (traj.max[0] /*max of dpot*/ < 0.0) return false;
    // if (traj.bbmax[2] - traj.bbmin[2] /*duration*/< 2) return false;
    return true;
  });

  if (ends_with(output_filename, str_vtp)) {
    // auto poly = tracker.get_discrete_critical_points_vtk();
    auto poly = tracker.get_traced_critical_points_vtk();
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
    writer->SetFileName(output_filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
  } else {
    tracker.write_traced_critical_points_text(output_filename);
  }

  return 0;
}
