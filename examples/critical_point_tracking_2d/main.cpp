#include <fstream>
#include <mutex>
#include <cassert>
#include <cxxopts.hpp>

#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
// #include <ftk/algorithms/cca.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <hypermesh/ndarray.hh>
#include <hypermesh/regular_simplex_mesh.hh>

#if FTK_HAVE_QT5
#include "widget.h"
#include <QApplication>
#endif

#if FTK_HAVE_VTK
#include <ftk/geometry/curve2vtk.hh>
#include <vtkPolyDataMapper.h>
#include <vtkTubeFilter.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#endif

// for serialization
#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>


int DW, DH; // the dimensionality of the data is DW*DH
int DT; // number of timesteps

hypermesh::ndarray<float> scalar, grad, hess;
hypermesh::regular_simplex_mesh m(3); // the 3D space-time mesh

std::mutex mutex;

struct intersection_t {
  size_t eid;
  float x[3]; // the spacetime coordinates of the trajectory
  float val; // scalar value at the intersection

  template <class Archive> void serialize(Archive & ar) {
    ar(eid, x[0], x[1], x[2], val);
  }
};
 
std::map<hypermesh::regular_simplex_mesh_element, intersection_t> intersections;

// the output trajectories
std::vector<std::vector<float>> trajectories;


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
  bool succ = ftk::inverse_lerp_s2v2(g, mu);
  float val = ftk::lerp_s2(value, mu);
  
  if (!succ) return;

  float hessxx[3], hessxy[3], hessyy[3];
  for (int i = 0; i < vertices.size(); i ++) {
    hessxx[i] = hess(0, 0, vertices[i][0], vertices[i][1], vertices[i][2]);
    hessxy[i] = hess(0, 1, vertices[i][0], vertices[i][1], vertices[i][2]);
    hessyy[i] = hess(1, 1, vertices[i][0], vertices[i][1], vertices[i][2]);
  }
  float hxx = ftk::lerp_s2(hessxx, mu),
        hxy = ftk::lerp_s2(hessxy, mu), 
        hyy = ftk::lerp_s2(hessyy, mu);
  float eig[2];
  ftk::solve_eigenvalues_symmetric2x2(hxx, hxy, hyy, eig);

  if (eig[0] < 0 && eig[1] < 0) { 
    float X[3][3];
    for (int i = 0; i < vertices.size(); i ++)
      for (int j = 0; j < 3; j ++)
        X[i][j] = vertices[i][j];

    intersection_t I;
    I.eid = f.to_integer();
    ftk::lerp_s2v3(X, mu, I.x);
    I.val = ftk::lerp_s2(value, mu);

    {
      std::lock_guard<std::mutex> guard(mutex);
      intersections[f] = I;
      // fprintf(stderr, "x={%f, %f}, t=%f, val=%f\n", I.x[0], I.x[1], I.x[2], I.val);
    }
  }
}


void extract_connected_components(std::vector<std::set<hypermesh::regular_simplex_mesh_element>>& components)
{
  typedef hypermesh::regular_simplex_mesh_element element_t;

  // Initialization
  ftk::union_find<std::string> uf; 
  std::map<std::string, element_t> id2ele; 
  for (const auto &f : intersections) {
    std::string eid = f.first.to_string(); 
    uf.add(eid); 
    id2ele.insert(std::make_pair(eid, f.first));
  }

  // Connected Component Labeling by using union-find. 
  m.element_for(3, [&](const hypermesh::regular_simplex_mesh_element& f) {
    const auto elements = f.sides();
    std::set<std::string> features; 

    for (const auto& ele : elements) {
      std::string eid = ele.to_string(); 

      if(uf.has(eid)) {
        features.insert(eid); 
      }
    }

    if(features.size()  > 1) {
      for(std::set<std::string>::iterator ite_i = std::next(features.begin(), 1); ite_i != features.end(); ++ite_i) {
        std::lock_guard<std::mutex> guard(mutex); // Use a lock for thread-save. 
        uf.unite(*(features.begin()), *ite_i); 
      }
    }
  }); 


  // Get disjoint sets of element IDs
  std::vector<std::set<std::string>> components_str;
  uf.get_sets(components_str);

  // Convert element IDs to elements
  for(auto& comp_str : components_str) {
    std::set<element_t> comp; 
    for(auto& ele_id : comp_str) {
      comp.insert(id2ele.find(ele_id)->second); 
    }

    components.push_back(comp); 
  }
}

void trace_intersections()
{
  typedef hypermesh::regular_simplex_mesh_element element_t; 

  std::vector<std::set<element_t>> cc; // connected components 
  extract_connected_components(cc);

  // Convert connected components to geometries

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

  for (int i = 0; i < cc.size(); i ++) {
    std::vector<std::vector<float>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(cc[i], neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      std::vector<float> mycurve, mycolors;
      for (int k = 0; k < linear_graphs[j].size(); k ++) {
        auto p = intersections[linear_graphs[j][k]];
        mycurve.push_back(p.x[0]); //  / (DW-1));
        mycurve.push_back(p.x[1]); //  / (DH-1));
        mycurve.push_back(p.x[2]); //  / (DT-1));
        mycurve.push_back(p.val);
      }
      trajectories.emplace_back(mycurve);
    }
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
    printf("--Curve %d:\n", i);
    const auto &curve = trajectories[i];
    for (int k = 0; k < curve.size()/4; k ++) {
      printf("---x=(%f, %f), t=%f, val=%f\n", curve[k*4], curve[k*4+1], curve[k*4+2], curve[k*4+3]);
    }
  }
}

void read_traj_file(const std::string& f)
{
  std::ifstream ifs(f);
  cereal::BinaryInputArchive ar(ifs);
  ar(trajectories);
  ifs.close();
}

void write_traj_file(const std::string& f)
{
  std::ofstream ofs(f);
  cereal::BinaryOutputArchive ar(ofs);
  ar(trajectories);
  ofs.close();
}

void read_dump_file(const std::string& f)
{
  std::vector<intersection_t> vector;

  std::ifstream ifs(f);
  cereal::BinaryInputArchive ar(ifs);
  ar(vector);
  ifs.close();

  for (const auto &i : vector) {
    hypermesh::regular_simplex_mesh_element e(m, 2, i.eid);
    intersections[e] = i;
  }
}

void write_dump_file(const std::string& f)
{
  std::vector<intersection_t> vector;
  for (const auto &i : intersections)
    vector.push_back(i.second);
  
  std::ofstream ofs(f);
  cereal::BinaryOutputArchive ar(ofs);
  ar(vector);
  ofs.close();
}

#if FTK_HAVE_VTK
void start_vtk_window()
{
  auto vtkcurves = ftk::curves2vtk(trajectories, 4);
  vtkcurves->Print(std::cerr);

  vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
  tubeFilter->SetInputData(vtkcurves);
  tubeFilter->SetRadius(1);
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
  std::string pattern, format;
  std::string filename_dump_r, filename_dump_w;
  std::string filename_traj_r, filename_traj_w;
  float threshold;
  bool show_qt = false, show_vtk = false;

  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "input file name pattern", cxxopts::value<std::string>(pattern))
    ("f,format", "input file format", cxxopts::value<std::string>(format)->default_value("float32"))
    ("read-dump", "read dump file", cxxopts::value<std::string>(filename_dump_r))
    ("write-dump", "write dump file", cxxopts::value<std::string>(filename_dump_w))
    ("read-traj", "read traj file", cxxopts::value<std::string>(filename_traj_r))
    ("write-traj", "write traj file", cxxopts::value<std::string>(filename_traj_w))
    ("w,width", "width", cxxopts::value<int>(DW)->default_value("128"))
    ("h,height", "height", cxxopts::value<int>(DH)->default_value("128"))
    ("t,timesteps", "timesteps", cxxopts::value<int>(DT)->default_value("10"))
    ("threshold", "threshold", cxxopts::value<float>(threshold)->default_value("0"))
    ("vtk", "visualization with vtk", cxxopts::value<bool>(show_vtk))
    ("qt", "visualization with qt", cxxopts::value<bool>(show_qt))
    ("d,debug", "enable debugging");
  auto results = options.parse(argc, argv);

  if (pattern.empty()) { // if the input data is not given, generate a synthetic data for the demo
    scalar = generate_synthetic_data<float>(DW, DH, DT);
  } else { // load the binary data
    scalar.reshape(DW, DH, 0);
    scalar.from_binary_file_sequence(pattern);
    DT = scalar.dim(2);
  }
 
  m.set_lb_ub({2, 2, 0}, {DW-3, DH-3, DT-1}); // update the mesh; set the lower and upper bounds of the mesh

  if (!filename_traj_r.empty()) { // if the trajectory file is given, skip all the analysis and visualize/print the trajectories
    read_traj_file(filename_traj_r);
  } else { // otherwise do the analysis
    if (!filename_dump_r.empty()) { // if the dump file is given, skill the sweep step; otherwise do sweep-and-trace
      read_dump_file(filename_dump_r);
    } else { // derive gradients and do the sweep
      grad = derive_gradients2(scalar);
      hess = derive_hessians2(grad);
      scan_intersections();
    }

    if (!filename_dump_w.empty())
      write_dump_file(filename_dump_w);

    trace_intersections();

    if (!filename_traj_w.empty())
      write_traj_file(filename_traj_w);
  }

  if (show_qt) {
#if FTK_HAVE_QT5
    QApplication app(argc, argv);
    QGLFormat fmt = QGLFormat::defaultFormat();
    fmt.setSampleBuffers(true);
    fmt.setSamples(16);
    QGLFormat::setDefaultFormat(fmt);

    CGLWidget *widget = new CGLWidget(scalar);
    widget->set_trajectories(trajectories, threshold);
    widget->show();
    return app.exec();
#else
    fprintf(stderr, "Error: the executable is not compiled with Qt\n");
#endif
  } else if (show_vtk) {
#if FTK_HAVE_VTK
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
