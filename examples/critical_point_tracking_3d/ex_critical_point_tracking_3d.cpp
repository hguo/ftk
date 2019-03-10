#include <mutex>

#if HAVE_NETCDF
#include <netcdf.h>
#endif

#if HAVE_VTK
#include <ftk/geometry/curve2vtk.hh>
#include <vtkPolyDataMapper.h>
#include <vtkImageData.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkTubeFilter.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#endif

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

const int DW = 128, DH = 128, DD = 128;// the dimensionality of the data is DW*DH
const int DT = 2; // number of timesteps

hypermesh::ndarray<float> scalar, grad, hess;
hypermesh::regular_simplex_mesh m(4); // the 4D space-time mesh

std::mutex mutex;

struct intersection_t {
  float x[4]; // the spacetime coordinates of the trajectory
};
 
std::map<hypermesh::regular_simplex_mesh_element, intersection_t> intersections;

// the output trajectories
std::vector<std::vector<std::vector<float>>> trajectories;

void derive_gradients()
{
  grad.reshape({3, DW, DH, DD, DT});
  for (int t = 0; t < DT; t ++) {
    for (int k = 1; k < DD-1; k ++) {
      for (int j = 1; j < DH-1; j ++) {
        for (int i = 1; i < DW-1; i ++) {
          grad(0, i, j, k, t) = 0.5 * (scalar(i+1, j, k, t) - scalar(i-1, j, k, t));
          grad(1, i, j, k, t) = 0.5 * (scalar(i, j+1, k, t) - scalar(i, j-1, k, t));
          grad(2, i, j, k, t) = 0.5 * (scalar(i, j, k+1, t) - scalar(i, j, k-1, t));
        }
      }
    }
  }
}

void derive_hessians()
{
  hess.reshape({3, 3, DW, DH, DD, DT});

  for (int t = 0; t < DT; t ++) {
    for (int k = 0; k < DD; k ++) {
      for (int j = 2; j < DH-2; j ++) {
        for (int i = 2; i < DW-2; i ++) {
          const float H00 = hess(0, 0, i, j, k, t) = // ddf/dx2
            0.5 * (grad(0, i+1, j, k, t) - grad(0, i-1, j, k, t));
          const float H01 = hess(0, 1, i, j, k, t) = // ddf/dxdy
            0.5 * (grad(0, i, j+1, k, t) - grad(0, i, j-1, k, t));
          const float H02 = hess(0, 2, i, j, k, t) = // ddf/dxdz
            0.5 * (grad(0, i, j, k+1, t) - grad(0, i, j, k-1, t));

          const float H10 = hess(1, 0, i, j, k, t) = // ddf/dydx
            0.5 * (grad(1, i+1, j, k, t) - grad(1, i-1, j, k, t));
          const float H11 = hess(1, 1, i, j, k, t) = // ddf/dy2
            0.5 * (grad(1, i, j+1, k, t) - grad(1, i, j-1, k, t));
          const float H12 = hess(1, 2, i, j, k, t) = // ddf/dydz
            0.5 * (grad(1, i, j, k+1, t) - grad(1, i, j, k-1, t));

          const float H20 = hess(2, 0, i, j, k, t) = // ddf/dydx
            0.5 * (grad(2, i+1, j, k, t) - grad(2, i-1, j, k, t));
          const float H21 = hess(2, 1, i, j, k, t) = // ddf/dy2
            0.5 * (grad(2, i, j+1, k, t) - grad(2, i, j-1, k, t));
          const float H22 = hess(2, 2, i, j, k, t) = // ddf/dydz
            0.5 * (grad(2, i, j, k+1, t) - grad(2, i, j, k-1, t));
        }
      }
    }
  }
}

void check_simplex(const hypermesh::regular_simplex_mesh_element& s)
{
  if (!s.valid()) return; // check if the 3-simplex is valid
  
  const auto &vertices = s.vertices();
  float X[4][4], g[4][3], value[4];

  for (int i = 0; i < 4; i ++) {
    for (int j = 0; j < 3; j ++)
      g[i][j] = grad(j, vertices[i][0], vertices[i][1], vertices[i][2], vertices[i][3]);
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
    value[i] = scalar(vertices[i][0], vertices[i][1], vertices[i][2], vertices[i][3]);
  }

  // check intersection
  float mu[4], x[4];
  bool succ = ftk::inverse_linear_interpolation_3simplex_vector3(g, mu);
  if (!succ) return;

  // check hessian
  float H[4][3][3], h[3][3];
  for (int i = 0; i < 4; i ++)
    for (int j = 0; j < 3; j ++)
      for (int k = 0; k < 3; k ++)
        H[i][j][k] = hess(j, k, vertices[i][0], vertices[i][1], vertices[i][2], vertices[i][3]);
  ftk::linear_interpolation_3simplex_matrix3x3(H, mu, h);

  float eig[3];
  ftk::solve_eigenvalues_symmetric3x3(h, eig);
  // fprintf(stderr, "eig=%f, %f, %f\n", eig[0], eig[1], eig[2]);

  if (eig[0] < 0 && eig[1] < 0 && eig[2] < 0) { // local maxima
    // dump results
    float val = ftk::linear_interpolation_3simplex(value, mu);
    ftk::linear_interpolation_3simplex_vector4(X, mu, x);
   
    intersection_t p;
    p.x[0] = x[0]; p.x[1] = x[1]; p.x[2] = x[2]; p.x[3] = x[3];
    {
      std::lock_guard<std::mutex> guard(mutex);
      intersections[s] = p;
    
      std::cerr << s << std::endl;
      fprintf(stderr, "x={%f, %f, %f, %f}\n", x[0], x[1], x[2], x[3]);
    }
  }
}

void trace_intersections()
{
  typedef hypermesh::regular_simplex_mesh_element element_t;

  std::set<element_t> qualified_elements;
  for (const auto &f : intersections)
    qualified_elements.insert(f.first);

  // std::function<std::set<face_t>(face_t)> neighbors = 
  auto neighbors = [](element_t e) {
    std::set<element_t> neighbors;
    const auto hypercells = e.side_of();
    for (const auto c : hypercells) {
      const auto sides = c.sides();
      for (const auto s : sides)
        neighbors.insert(s);
    }
    return neighbors;
  };

  auto cc = ftk::extract_connected_components<element_t, std::set<element_t>>(neighbors, qualified_elements);
  fprintf(stderr, "#cc=%lu\n", cc.size());

  for (int i = 0; i < cc.size(); i ++) {
    // auto mycolor = QColor::fromHslF((float)rand()/RAND_MAX, 0.5, 0.5);
    // cc_colors.push_back(mycolor);

    std::vector<std::vector<float>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(cc[i], neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      std::vector<float> mycurve, mycolors;
      // fprintf(stderr, "curve:\n");
      for (int k = 0; k < linear_graphs[j].size(); k ++) {
        auto p = intersections[linear_graphs[j][k]];
        mycurve.push_back(p.x[0] / (DW-1));
        mycurve.push_back(p.x[1] / (DH-1));
        mycurve.push_back(p.x[2] / (DD-1));
        // fprintf(stderr, "p={%f, %f, %f}\n", p.x[0], p.x[1], p.x[2]);
      }
      mycurves.emplace_back(mycurve);
      // tubes.push_back(ftk::curve2tube<float>(mycurve, mycolors, 12, 0.005));
    }
    trajectories.emplace_back(mycurves);
  }

#if 0
  int i = 0;
  for (const auto c : cc) {
    for (const auto f : c)
      intersections[f].label = i;
    i ++;
  }
#endif
}

void track_critical_points()
{
  m.set_lb_ub({2, 2, 2, 0}, {DW-3, DH-3, DD-3, DT-1}); // set the lower and upper bounds of the mesh
  m.element_for(3, check_simplex); // iterate over all 3-simplices

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
        printf("---x=(%f, %f, %f)\n", curve[k*3], curve[k*3+1], curve[k*3+2]);
      }
    }
  }
}

#if HAVE_VTK
void start_vtk_window()
{
  // initialize volume data
  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
  imageData->SetDimensions(DW, DH, DD);
  imageData->AllocateScalars(VTK_FLOAT, 1);

  float *ptr = static_cast<float*>(imageData->GetScalarPointer(0, 0, 0));
  for (int i = 0; i < DW*DH*DD; i ++)
    ptr[i] = scalar[i];

  vtkSmartPointer<vtkVolume> volume = vtkVolume::New();
  vtkSmartPointer<vtkSmartVolumeMapper> volumeMapper = vtkSmartVolumeMapper::New();
  vtkSmartPointer<vtkColorTransferFunction> colorFun = vtkColorTransferFunction::New();
  vtkSmartPointer<vtkPiecewiseFunction> opacityFun = vtkPiecewiseFunction::New();

  vtkSmartPointer<vtkVolumeProperty> property = vtkVolumeProperty::New();
  property->SetColor(colorFun);
  property->SetScalarOpacity(opacityFun);
  property->SetInterpolationTypeToLinear();

  volumeMapper->SetInputData(imageData);
  volume->SetProperty(property);
  volume->SetMapper(volumeMapper);

  // render curves
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

  renderer->SetUseDepthPeelingForVolumes(true);
  renderer->AddVolume(volume);

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
  scalar.reshape(DW, DH, DD, DT);
  scalar.from_netcdf(argv[1], "vort", starts, sizes);

#if 1
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
