// #include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/mesh/simplex_2d_mesh.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>
#include <ftk/external/cxxopts.hpp>

std::string input_filename, 
  mesh_filename, 
  kernel_filename = "xgc.kernel";
double sigma(0.01);

void parse_arguments(int argc, char **argv)
{
  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "Input file name pattern: a single file or a series of file, e.g. 'scalar.raw', 'cm1out_000*.nc'",
     cxxopts::value<std::string>(input_filename))
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
  ftk::simplex_2d_mesh<> m(coords, triangles);
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

  if (input_filename.length()) {
    ftk::ndarray<double> dpot;
    dpot.from_h5(input_filename, "/dpot");
    dpot = dpot.transpose();
    dpot.reshape(dpot.dim(0));

    auto dpot1 = m.smooth_scalar_field(dpot);
    m.scalar_to_vtk_unstructured_grid_data_file("out.vtu", "dpot", dpot1);
  }

  return 0;
}
