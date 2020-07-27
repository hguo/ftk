// #include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/mesh/simplex_2d_mesh.hh>
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

std::string input_filename, 
  mesh_filename, 
  kernel_filename = "xgc.kernel";
double sigma(0.02);

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

    auto grad = m.smooth_gradient_field(dpot);
    grad *= (1.0/(sigma*sigma));
    m.vector_to_vtk_unstructured_grid_data_file("out-grad.vtu", "ddpot", grad);
    // std::cerr << grad << std::endl;

    std::vector<double> cps;

    m.element_for(2, [&](int i) {
      int tri[3];
      m.get_triangle(i, tri);

      // typedef ftk::fixed_point<> fp_t;
      typedef mpf_class fp_t;
      double V[3][2], X[3][2];
      fp_t Vf[3][2];
      for (int k = 0; k < 3; k ++)
        for (int j = 0; j < 2; j ++) {
          V[k][j] = grad(j, tri[k]);
          Vf[k][j] = V[k][j];
          X[k][j] = coords(j, tri[k]);
        }

      bool succ = ftk::robust_critical_point_in_simplex2(Vf, tri);
      if (!succ) return;

      // ftk::print3x2("V", V);
      double mu[3], x[2];
      bool succ2 = ftk::inverse_lerp_s2v2(V, mu);
      // if (!succ2) return;
  
      ftk::lerp_s2v2(X, mu, x);
      fprintf(stderr, "mu=%f, %f, %f, x=%f, %f\n", 
          mu[0], mu[1], mu[2], x[0], x[1]);
      cps.push_back(x[0]);
      cps.push_back(x[1]);
    });

    auto poly = ftk::points2vtk(cps, 2);
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
    writer->SetFileName("out.vtp");
    writer->SetInputData(poly);
    writer->Write();
  }

  return 0;
}
