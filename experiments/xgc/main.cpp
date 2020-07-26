// #include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/mesh/simplex_2d_mesh.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

#if FTK_HAVE_QT5
#include <QApplication>
#include "widget.h"
#endif

int main(int argc, char **argv)
{
  if (argc < 3) return 1; // mesh.h5, dpot.h5

  // load mesh & data from hdf5
  ftk::ndarray<int> triangles;
  triangles.from_h5(argv[1], "/cell_set[0]/node_connect_list");
  // std::cerr << triangles << std::endl; 

  ftk::ndarray<double> coords;
  coords.from_h5(argv[1], "/coordinates/values");
  // std::cerr << coords << std::endl;

  ftk::ndarray<double> dpot;
  dpot.from_h5(argv[2], "/dpot");
  dpot = dpot.transpose();
  // std::cerr << dpot << std::endl;

  ftk::simplex_2d_mesh<> m(coords, triangles);
#if 0
  m.build_edges();
  m.build_smoothing_kernel(0.04);
  m.smooth_scalar_field(dpot);
#endif

#if FTK_HAVE_QT5
  QApplication app(argc, argv);
  QGLFormat fmt = QGLFormat::defaultFormat();
  fmt.setSampleBuffers(true);
  fmt.setSamples(16); 
  QGLFormat::setDefaultFormat(fmt); 
  
  CGLWidget *widget = new CGLWidget;
  widget->set_mesh(&m);
#if 0
  widget->loadMeshFromJsonFile("xgc.mesh.json");
  widget->loadBranchesFromJsonFile("xgc.branches.json");
  widget->loadLabels("xgc.labels.bin"); // TODO: load labels from ADIOS
#endif
  widget->show();
  return app.exec();
#else
  return 0;
#endif
}
