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

  // build extruded mesh
  // ftk::simplicial_unstructured_extruded_2d_mesh<> m1(m);
  // fprintf(stderr, "extrusion mesh built.\n");

  // smoothing data; extract cps
  if (input_filename.length()) {
    ftk::ndarray<double> data;
    data.from_h5(input_filename, "/dpot");
    data = data.transpose();
    // dpot.reshape(dpot.dim(0));

    ftk::critical_point_tracker_2d_unstructured tracker(m);

    const int nt = data.dim(1);
    for (int t = 0; t < nt; t ++) {
      fprintf(stderr, "current timestep %d\n", t);
      ftk::ndarray<double> slice;
      slice.reshape(data.dim(0));
      for (auto j = 0; j < data.dim(0); j ++)
        slice[j] = data(j, t);

      ftk::ndarray<double> scalar, grad, J;
      m.smooth_scalar_gradient_jacobian(slice, sigma, scalar, grad, J);
   
      tracker.push_field_data_snapshot(scalar, grad, J);

      if (t != 0) tracker.advance_timestep();
      else tracker.update_timestep();
    }

#if 0
    std::vector<double> cps;
    std::vector<unsigned int> types;

    m1.element_for_ordinal(2, 0, [&](int i) {
      int tri[3];
      m1.get_simplex(2, i, tri); 

      // typedef ftk::fixed_point<> fp_t;
      typedef mpf_class fp_t;
      double V[3][2], X[3][3];
      fp_t Vf[3][2];
      for (int k = 0; k < 3; k ++) {
        int t = tri[k] >= m.n(0) ? 1 : 0;
        int v = tri[k] % m.n(0);
        for (int j = 0; j < 2; j ++) {
          V[k][j] = grad[t](j, v);
          Vf[k][j] = V[k][j];
          X[k][j] = coords(j, v);
        }
        X[k][2] = t;
      }

      bool succ = ftk::robust_critical_point_in_simplex2(Vf, tri);
      if (!succ) return;

      // ftk::print3x2("V", V);
      double mu[3], x[3];
      bool succ2 = ftk::inverse_lerp_s2v2(V, mu);
      // if (!succ2) return;
      ftk::lerp_s2v3(X, mu, x);

      double Js[3][2][2], H[2][2];
      for (int k = 0; k < 3; k ++) {
        int t = tri[k] >= m.n(0) ? 1 : 0;
        int v = tri[k] % m.n(0);
        for (int j = 0; j < 2; j ++)
          for (int i = 0; i < 2; i ++)
            Js[k][j][i] = J[t](i, j, v); 
      }
      ftk::lerp_s2m2x2(Js, mu, H);
      // ftk::print2x2("H", H);
      const int type = ftk::critical_point_type_2d(H, true);

      fprintf(stderr, "mu=%f, %f, %f, x=%f, %f, %f, type=%d\n", 
          mu[0], mu[1], mu[2], x[0], x[1], x[2], type);
      cps.push_back(x[0]);
      cps.push_back(x[1]);
      cps.push_back(x[2]);
      types.push_back(type);
    });

    {
      auto poly = ftk::points2vtk(cps, 3);
      vtkSmartPointer<vtkUnsignedIntArray> vtypes = vtkUnsignedIntArray::New();
      vtypes->SetNumberOfValues(types.size());
      for (auto i = 0; i < types.size(); i ++)
        vtypes->SetValue(i, types[i]);
      vtypes->SetName("type");
      poly->GetPointData()->AddArray(vtypes);

      vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
      writer->SetFileName("out.vtp");
      writer->SetInputData(poly);
      writer->Write();
    }
#endif
  }

  return 0;
}
