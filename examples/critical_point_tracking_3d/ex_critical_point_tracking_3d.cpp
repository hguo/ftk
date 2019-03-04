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

const int DW = 128, DH = 128, DD = 128;// the dimensionality of the data is DW*DH
const int DT = 10; // number of timesteps

hypermesh::ndarray<float> scalar, grad, hess;
hypermesh::regular_simplex_mesh m(3); // the 3D space-time mesh

std::mutex mutex;

struct intersection_t {
  float x[3]; // the spacetime coordinates of the trajectory
};
 
std::map<hypermesh::regular_simplex_mesh_element, intersection_t> intersections;

// the output trajectories
std::vector<std::vector<std::vector<float>>> trajectories;

void derive_gradients()
{
  grad.reshape({2, (size_t)DW, (size_t)DH, (size_t)DT});
  for (int k = 0; k < DT; k ++) {
    for (int j = 1; j < DH-1; j ++) {
      for (int i = 1; i < DW-1; i ++) {
        grad(0, i, j, k) = 0.5 * (scalar(i+1, j, k) - scalar(i-1, j, k)) * (DW-1);
        grad(1, i, j, k) = 0.5 * (scalar(i, j+1, k) - scalar(i, j-1, k)) * (DH-1);
      }
    }
  }
}

void derive_hessians()
{
  hess.reshape({2, 2, (size_t)DW, (size_t)DH, (size_t)DT});

  for (int k = 0; k < DT; k ++) {
    for (int j = 2; j < DH-2; j ++) {
      for (int i = 2; i < DW-2; i ++) {
        const float H00 = hess(0, 0, i, j, k) = // ddf/dx2
          0.5 * (grad(0, i+1, j, k) - grad(0, i-1, j, k)) * (DW-1) / 15;
        const float H01 = hess(0, 1, i, j, k) = // ddf/dxdy
          0.5 * (grad(0, i, j+1, k) - grad(0, i, j-1, k)) * (DH-1) / 15;
        const float H10 = hess(1, 0, i, j, k) = // ddf/dydx
          0.5 * (grad(1, i+1, j, k) - grad(1, i-1, j, k)) * (DW-1) / 15;
        const float H11 = hess(1, 1, i, j, k) = // ddf/dy2
          0.5 * (grad(1, i, j+1, k) - grad(1, i, j-1, k)) * (DH-1) / 15;
      }
    }
  }
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

void track_critical_points()
{
  m.set_lb_ub({1, 1, 0}, {DW-2, DH-2, DT-1}); // set the lower and upper bounds of the mesh
  m.element_for(2, check_simplex); // iterate over all 2-simplices

  trace_intersections();
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
  size_t starts[4] = {0, 0, 0, 0}, 
         sizes[4]  = {size_t(DT), size_t(DD), size_t(DH), size_t(DW)};
  scalar.from_netcdf(argv[1], "vorts", starts, sizes);

#if 0
  derive_gradients();
  derive_hessians();
  track_critical_points();
#endif

#if HAVE_VTK
  start_vtk_window();
#else
  print_trajectories();
#endif

  return 0;
}
