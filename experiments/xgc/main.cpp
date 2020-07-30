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
#include <gmpxx.h>
#include <vtkXMLPolyDataWriter.h>

std::string input_filename_pattern,
  mesh_filename, 
  kernel_filename = "xgc.kernel";
std::vector<std::string> input_filenames;
double sigma(0.02);

void parse_arguments(int argc, char **argv)
{
  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "Input file name pattern: a single file or a series of file, e.g. 'scalar.raw', 'cm1out_000*.nc'",
     cxxopts::value<std::string>(input_filename_pattern))
    ("m,mesh", "Input mesh file", cxxopts::value<std::string>(mesh_filename))
    ("k,kernel", "Input/output smoothing kernel file", cxxopts::value<std::string>(kernel_filename))
    ("s,sigma", "Kernel bandwidth", cxxopts::value<double>(sigma));
  auto results = options.parse(argc, argv);

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
  ftk::ndarray<double> coords;
  
  triangles.from_h5(mesh_filename, "/cell_set[0]/node_connect_list");
  coords.from_h5(mesh_filename, "/coordinates/values");
 
  // build mesh
  ftk::simplicial_unstructured_2d_mesh<> m(coords, triangles);
  m.build_edges();
  m.build_vertex_links();

  if (kernel_filename.length()) {
    bool succ = m.read_smoothing_kernel(kernel_filename);
    if (!succ) {
      m.build_smoothing_kernel(sigma);
      m.write_smoothing_kernel(kernel_filename);
    }
  }

  fprintf(stderr, "mesh loaded.\n");

  input_filenames = ftk::ndarray<double>::glob(input_filename_pattern);
  
  ftk::critical_point_tracker_2d_unstructured tracker(m);
  for (int t = 0; t < std::min(3, (int)input_filenames.size()); t ++) {
    ftk::ndarray<double> data;
    data.from_h5(input_filenames[t], "/dpot");
    data = data.transpose();
    data.reshape(data.dim(0)); // only use the first slice

    ftk::ndarray<double> scalar, grad, J;
    m.smooth_scalar_gradient_jacobian(data, sigma, scalar, grad, J);
 
    tracker.push_field_data_snapshot(scalar, grad, J);

    if (t != 0) tracker.advance_timestep();
    if (t == input_filenames.size()-1) tracker.update_timestep();
  }
  tracker.finalize();
   
  // auto poly = tracker.get_discrete_critical_points_vtk();
  auto poly = tracker.get_traced_critical_points_vtk();
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName("out.vtp");
  writer->SetInputData(poly);
  writer->Write();

  return 0;
}
