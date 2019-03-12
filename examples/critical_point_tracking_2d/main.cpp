#include <cxxopts.hpp>
#include <mutex>

#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <hypermesh/ndarray.hh>
#include <hypermesh/regular_simplex_mesh.hh>

#if HAVE_QT
#include "widget.h"
#include <QApplication>
#endif

#if HAVE_VTK
#include <ftk/geometry/curve2vtk.hh>
#include <vtkPolyDataMapper.h>
#include <vtkTubeFilter.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#endif

int DW, DH; // the dimensionality of the data is DW*DH
int DT; // number of timesteps

hypermesh::ndarray<float> scalar, grad, hess;
hypermesh::regular_simplex_mesh m(3); // the 3D space-time mesh

std::mutex mutex;

struct intersection_t {
  float x[3]; // the spacetime coordinates of the trajectory
};
 
std::map<hypermesh::regular_simplex_mesh_element, intersection_t> intersections;

// the output trajectories
std::vector<std::vector<std::vector<float>>> trajectories;


template <typename T> // the synthetic function
T f(T x, T y, T t) 
{
  return cos(x*cos(t)-y*sin(t))*sin(x*sin(t)+y*cos(t));
}

template <typename T>
hypermesh::ndarray<T> generate_synthetic_data(int DW, int DH, int DT)
{
  hypermesh::ndarray<T> scalar;
  scalar.reshape(DW, DH, DT);

  const T scaling_factor = 15; // the factor that controls the shape of the synthesize data
  for (int k = 0; k < DT; k ++) {
    for (int j = 0; j < DH; j ++) {
      for (int i = 0; i < DW; i ++) {
        const T x = ((T(i) / (DW-1)) - 0.5) * scaling_factor,
                y = ((T(j) / (DH-1)) - 0.5) * scaling_factor, 
                t = (T(k) / (DT-1)) + 1e-4;
        scalar(i, j, k) = f(x, y, t);
      }
    }
  }

  return scalar;
}

template <typename T>
hypermesh::ndarray<T> derive_gradients2(const hypermesh::ndarray<T>& scalar)
{
  hypermesh::ndarray<T> grad;
  grad.reshape(2, scalar.dim(0), scalar.dim(1), scalar.dim(2));
  
  for (int k = 0; k < DT; k ++) {
    for (int j = 1; j < DH-1; j ++) {
      for (int i = 1; i < DW-1; i ++) {
        grad(0, i, j, k) = 0.5 * (scalar(i+1, j, k) - scalar(i-1, j, k)) * (DW-1);
        grad(1, i, j, k) = 0.5 * (scalar(i, j+1, k) - scalar(i, j-1, k)) * (DH-1);
      }
    }
  }
  return grad;
}

template <typename T>
hypermesh::ndarray<T> derive_hessians2(const hypermesh::ndarray<T>& grad)
{
  hypermesh::ndarray<T> hess;
  hess.reshape(2, grad.dim(0), grad.dim(1), grad.dim(2), grad.dim(3));

  for (int k = 0; k < DT; k ++) {
    for (int j = 2; j < DH-2; j ++) {
      for (int i = 2; i < DW-2; i ++) {
        const T H00 = hess(0, 0, i, j, k) = // ddf/dx2
          0.5 * (grad(0, i+1, j, k) - grad(0, i-1, j, k)) * (DW-1) / 15;
        const T H01 = hess(0, 1, i, j, k) = // ddf/dxdy
          0.5 * (grad(0, i, j+1, k) - grad(0, i, j-1, k)) * (DH-1) / 15;
        const T H10 = hess(1, 0, i, j, k) = // ddf/dydx
          0.5 * (grad(1, i+1, j, k) - grad(1, i-1, j, k)) * (DW-1) / 15;
        const T H11 = hess(1, 1, i, j, k) = // ddf/dy2
          0.5 * (grad(1, i, j+1, k) - grad(1, i, j-1, k)) * (DH-1) / 15;
      }
    }
  }
  return hess;
}

void check_simplex(const hypermesh::regular_simplex_mesh_element& f)
{
  if (!f.valid()) return; // check if the 2-simplex is valid
  const auto &vertices = f.vertices(); // obtain the vertices of the simplex
  float g[3][2], value[3];

  for (int i = 0; i < 3; i ++) {
    g[i][0] = grad(0, vertices[i][0], vertices[i][1], vertices[i][2]);
    g[i][1] = grad(1, vertices[i][0], vertices[i][1], vertices[i][2]);
    value[i] = scalar(vertices[i][0], vertices[i][1], vertices[i][2]);
  }
 
  float mu[3];
  bool succ = ftk::inverse_linear_interpolation_2simplex_vector2(g, mu);
  float val = ftk::linear_interpolation_2simplex(value, mu);
  
  if (!succ) return;

  float hessxx[3], hessxy[3], hessyy[3];
  for (int i = 0; i < vertices.size(); i ++) {
    hessxx[i] = hess(0, 0, vertices[i][0], vertices[i][1], vertices[i][2]);
    hessxy[i] = hess(0, 1, vertices[i][0], vertices[i][1], vertices[i][2]);
    hessyy[i] = hess(1, 1, vertices[i][0], vertices[i][1], vertices[i][2]);
  }
  float hxx = ftk::linear_interpolation_2simplex(hessxx, mu),
        hxy = ftk::linear_interpolation_2simplex(hessxy, mu), 
        hyy = ftk::linear_interpolation_2simplex(hessyy, mu);
  float eig[2];
  ftk::solve_eigenvalues_symmetric2x2(hxx, hxy, hyy, eig);

  if (eig[0] < 0 && eig[1] < 0) { 
    float X[3][3];
    for (int i = 0; i < vertices.size(); i ++)
      for (int j = 0; j < 3; j ++)
        X[i][j] = vertices[i][j];

    intersection_t intersection;
    ftk::linear_interpolation_2simplex_vector3(X, mu, intersection.x);

    {
      std::lock_guard<std::mutex> guard(mutex);
      intersections[f] = intersection;
    }
  }
}

void trace_intersections()
{
  typedef hypermesh::regular_simplex_mesh_element element_t;

  std::set<element_t> qualified_elements;
  for (const auto &f : intersections)
    qualified_elements.insert(f.first);

  auto neighbors = [](element_t f) {
    std::set<element_t> neighbors;
    const auto cells = f.side_of();
    for (const auto c : cells) {
      const auto elements = c.sides();
      for (const auto f1 : elements)
        neighbors.insert(f1);
    }
    return neighbors;
  };

  // connected components
  auto cc = ftk::extract_connected_components<element_t, std::set<element_t>>(neighbors, qualified_elements);
  // fprintf(stderr, "#cc=%lu\n", cc.size());

  for (int i = 0; i < cc.size(); i ++) {

    std::vector<std::vector<float>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(cc[i], neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      std::vector<float> mycurve, mycolors;
      for (int k = 0; k < linear_graphs[j].size(); k ++) {
        auto p = intersections[linear_graphs[j][k]];
        mycurve.push_back(p.x[0] / (DW-1));
        mycurve.push_back(p.x[1] / (DH-1));
        mycurve.push_back(p.x[2] / (DT-1));
      }
      mycurves.emplace_back(mycurve);
    }
    trajectories.emplace_back(mycurves);
  }
}

void scan_intersections() 
{
  m.element_for(2, check_simplex); // iterate over all 2-simplices
}

void print_trajectories()
{
  printf("We found %lu trajectories:\n", trajectories.size());
  for (int i = 0; i < trajectories.size(); i ++) {
    const auto &curves = trajectories[i];
    printf("-Trajectory %d has %lu components:\n", i, curves.size());
    for (int j = 0; j < curves.size(); j ++) {
      const auto &curve = curves[j];
      printf("--Curve %d:\n", j);
      for (int k = 0; k < curve.size()/3; k ++) {
        printf("---x=(%f, %f), t=%f\n", curve[k*3], curve[k*3+1], curve[k*3+2]);
      }
    }
  }
}

void read_dump_file(const std::string& f)
{
  FILE *fp = fopen(f.c_str(), "rb");
  if (!fp) return;

  fseek(fp, 0L, SEEK_END);
  const auto sz = ftell(fp);
  const auto n = sz / sizeof(std::pair<size_t, intersection_t>);
  fseek(fp, 0L, SEEK_SET);
  
  std::vector<std::pair<size_t, intersection_t>> vector;
  vector.resize(n);

  fread((void*)(&vector[0]), sizeof(std::pair<size_t, intersection_t>), vector.size(), fp);
  fclose(fp);

  for (const auto &i : vector) {
    hypermesh::regular_simplex_mesh_element e(m, 2, i.first);
    intersections[e] = i.second;
  }
}

void write_dump_file(const std::string& f)
{
  std::vector<std::pair<size_t, intersection_t>> vector;
  for (const auto &i : intersections)
    vector.push_back(std::make_pair(i.first.to_integer<size_t>(), i.second));

  FILE *fp = fopen(f.c_str(), "wb");
  fwrite((void*)(&vector[0]), sizeof(std::pair<size_t, intersection_t>), vector.size(), fp);
  fclose(fp);
}

#if HAVE_VTK
void start_vtk_window()
{
  auto vtkcurves = ftk::curves2vtk(trajectories);
  // vtkcurves->Print(std::cerr);

  vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
  tubeFilter->SetInputData(vtkcurves);
  tubeFilter->SetRadius(.01);
  tubeFilter->SetNumberOfSides(50);
  tubeFilter->Update();
  
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  // mapper->SetInputData(vtkcurves);
  mapper->SetInputConnection(tubeFilter->GetOutputPort());

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // a renderer and render window
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // add the actors to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(1, 1, 1); // Background color white

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = 
      vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle( style );

  renderWindowInteractor->Start();
}
#endif

int main(int argc, char **argv)
{
  std::string filename, format;
  std::string filename_dump_r, filename_dump_w;
  bool show_qt = false, show_vtk = false;

  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "input file name", cxxopts::value<std::string>(filename))
    ("read-dump", "read dump file", cxxopts::value<std::string>(filename_dump_r))
    ("write-dump", "write dump file", cxxopts::value<std::string>(filename_dump_w))
    ("f,format", "input file format", cxxopts::value<std::string>(format))
    ("w,width", "width", cxxopts::value<int>(DW)->default_value("128"))
    ("h,height", "height", cxxopts::value<int>(DH)->default_value("128"))
    ("t,timesteps", "timesteps", cxxopts::value<int>(DT)->default_value("10"))
    ("vtk", "visualization with vtk", cxxopts::value<bool>(show_vtk))
    ("qt", "visualization with qt", cxxopts::value<bool>(show_qt))
    ("d,debug", "enable debugging");

  auto results = options.parse(argc, argv);
  // fprintf(stderr, "dims=%d, %d, %d\n", DW, DH, DT);

  scalar = generate_synthetic_data<float>(DW, DH, DT);
  m.set_lb_ub({1, 1, 0}, {DW-2, DH-2, DT-1}); // set the lower and upper bounds of the mesh

  if (!filename_dump_r.empty()) {
    read_dump_file(filename_dump_r);
  } else {
    grad = derive_gradients2(scalar);
    hess = derive_hessians2(grad);
    scan_intersections();
  }

  if (!filename_dump_w.empty())
    write_dump_file(filename_dump_w);

  trace_intersections();

  if (show_qt) {
#if HAVE_QT
    QApplication app(argc, argv);
    QGLFormat fmt = QGLFormat::defaultFormat();
    fmt.setSampleBuffers(true);
    fmt.setSamples(16);
    QGLFormat::setDefaultFormat(fmt);

    CGLWidget *widget = new CGLWidget(scalar);
    widget->show();
    return app.exec();
#else
    fprintf(stderr, "Error: the executable is not compiled with Qt\n");
#endif
  } else if (show_vtk) {
#if HAVE_VTK
    start_vtk_window();
    // ftk::write_curves_vtk(trajectories, "trajectories.vtp");
#else
    fprintf(stderr, "Error: the executable is not compiled with VTK\n");
#endif
  } else {
    print_trajectories();
  }

  return 0;
}