#include <ftk/ftk_config.hh>
#include <ftk/basic/union_find.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/fmod.hh>
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
#include <cxxopts.hpp>
#include "io/GLHeader.h"
#include "io/GLGPU_IO_Helper.h"

#include "puncture.h"

#if FTK_HAVE_QT5
#include "widget.h"
#include <QApplication>
#endif

GLHeader hdr;
hypermesh::ndarray<float> Re, Im, Rho, Phi;
hypermesh::regular_simplex_mesh m(3);

std::mutex mutex;

// std::map<hypermesh::regular_simplex_mesh_element, punctured_face_t> punctures;
std::map<size_t, punctured_face_t> punctures;
ftk::union_find<size_t> uf;
std::map<size_t, std::set<size_t>> links;
std::vector<std::vector<float>> vortices;

bool load_data(const std::string& filename)
{
  float *rho = NULL, *phi = NULL, *re = NULL, *im = NULL;

  if (GLGPU_IO_Helper_ReadBDAT(filename, hdr, &rho, &phi, &re, &im, NULL, NULL, NULL, false, false)) {
  } else if (GLGPU_IO_Helper_ReadLegacy(filename, hdr, &rho, &phi, &re, &im, NULL, NULL, NULL, false, false)) {
  } else return false;

  const size_t n = hdr.dims[0] * hdr.dims[1] * hdr.dims[2];
  const std::vector<size_t> shape = {(size_t)hdr.dims[0], (size_t)hdr.dims[1], (size_t)hdr.dims[2]};

  Re.copy(re, re+n);
  Im.copy(im, im+n);
  Rho.copy(rho, rho+n);
  Phi.copy(phi, phi+n);

  Re.reshape(shape);
  Im.reshape(shape);
  Rho.reshape(shape);
  Phi.reshape(shape);
  
  return true;
}

void magnetic_potential(const GLHeader &h, const float X[3], float A[3]) 
{
  if (h.B[1] > 0) {
    A[0] = -h.Kex;
    A[1] = X[0] * h.B[2];
    A[2] = -X[0] * h.B[1];
  } else {
    A[0] = -X[1] * h.B[2] - h.Kex;
    A[1] = 0;
    A[2] = X[1] * h.B[0];
  }
}

float line_integral(const float X0[], const float X1[], const float A0[], const float A1[])
{
  float dX[3] = {X1[0] - X0[0], X1[1] - X0[1], X1[2] - X0[2]};
  float A[3] = {A0[0]+A1[0], A0[1]+A1[1], A0[2]+A1[2]};

  return 0.5 * ftk::inner_product3(A, dX);
}

void extract_vortices()
{
  m.set_lb_ub({0, 0, 0}, {hdr.dims[0]-2, hdr.dims[1]-2, hdr.dims[2]-2});

  fprintf(stderr, "sweeping 2-simplices...\n");
  m.element_for(2, [&](const hypermesh::regular_simplex_mesh_element &f) {
      const auto &vertices = f.vertices();
      float X[3][3], A[3][3];
      float rho[3], phi[3], re[3], im[3];

      for (int i = 0; i < 3; i ++) {
        const size_t index = Rho.index(vertices[i]);
        rho[i] = Rho[index];
        phi[i] = Phi[index];
        re[i] = Re[index];
        im[i] = Im[index];
        for (int j = 0; j < 3; j ++) 
          X[i][j] = vertices[i][j] * hdr.cell_lengths[j] + hdr.origins[j];
        magnetic_potential(hdr, X[i], A[i]);
      }

      float delta[3], phase_shift = 0;
      for (int i = 0; i < 3; i ++) {
        const int j = (i+1) % 3;
        delta[i] = phi[j] - phi[i];
        float li = line_integral(X[i], X[j], A[i], A[j]);
        delta[i] = ftk::mod2pi1(delta[i] - li);
        phase_shift -= delta[i];
      }

      float critera = phase_shift / (2 * M_PI);
      if (fabs(critera) < 0.5) return;

      for (int i = 0; i < 3; i ++) {
        if (i != 0) phi[i] = phi[i-1] + delta[i-1];
        re[i] = rho[i] * cos(phi[i]);
        im[i] = rho[i] * sin(phi[i]);
      }

      float mu[3], pos[3];
      float re_im[3][2] = {{re[0], im[0]}, {re[1], im[1]}, {re[2], im[2]}};
      ftk::inverse_lerp_s2v2(re_im, mu);
      ftk::lerp_s2v3(X, mu, pos);
      // fprintf(stderr, "mu={%f, %f, %f}, pos={%f, %f, %f}\n", mu[0], mu[1], mu[2], pos[0], pos[1], pos[2]);

      punctured_face_t puncture;
      puncture.x[0] = pos[0];
      puncture.x[1] = pos[1];
      puncture.x[2] = pos[2];

      {
        std::lock_guard<std::mutex> guard(mutex);
        punctures[f.to_integer()] = puncture;
      }
  });

  fprintf(stderr, "sweeping 3-simplices...\n");
  m.element_for(3, [&](const hypermesh::regular_simplex_mesh_element &c) {
      const auto &sides = c.sides();
      std::vector<size_t> vector;
      for (const auto& f : sides)
        if (punctures.find(f.to_integer()) != punctures.end())
          vector.push_back(f.to_integer());
      if (vector.size() == 2) {
        std::lock_guard<std::mutex> guard(mutex);
        uf.add(vector[0]);
        uf.add(vector[1]);
        uf.unite(vector[0], vector[1]);

        links[vector[0]].insert(vector[1]);
        links[vector[1]].insert(vector[0]);
      }
  });

  const auto cc = uf.get_sets();
  fprintf(stderr, "number of connected components: %zu\n", cc.size());

  auto neighbors = [&](size_t i) {
    std::set<size_t> neighbors;
    auto iterator = links.find(i);
    if (iterator == links.end()) return neighbors;
    else return iterator->second;
  };

  for (int i = 0; i < cc.size(); i ++) {
    std::vector<std::vector<float>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<size_t>(cc[i], neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      std::vector<float> mycurve; // , mycolors;
      for (int k = 0; k < linear_graphs[j].size(); k ++) {
        auto p = punctures[linear_graphs[j][k]];
        mycurve.push_back(p.x[0]);
        mycurve.push_back(p.x[1]);
        mycurve.push_back(p.x[2]);
      }
      vortices.emplace_back(mycurve);
    }
  }
}

void print_vortices()
{
  printf("We found %lu vortices:\n", vortices.size());
  for (int i = 0; i < vortices.size(); i ++) {
    printf("--Curve %d:\n", i);
    const auto &curve = vortices[i];
    for (int k = 0; k < curve.size()/4; k ++) {
      printf("---x=(%f, %f, %f)\n", curve[k*3], curve[k*3+1], curve[k*3+2]);
    }
  }
}

int main(int argc, char **argv)
{
  std::string input_filename;
  bool vis_qt = false;

  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "input file", cxxopts::value<std::string>(input_filename))
    ("qt", "visualization with qt", cxxopts::value<bool>(vis_qt));
  options.parse_positional("input");
  auto results = options.parse(argc, argv);

  bool succ = load_data(input_filename);
  if (succ) {
    extract_vortices();
    print_vortices();

    if (vis_qt) {
#if FTK_HAVE_QT5
      QApplication app(argc, argv);
      QGLFormat fmt = QGLFormat::defaultFormat();
      fmt.setSampleBuffers(true);
      fmt.setSamples(16);
      QGLFormat::setDefaultFormat(fmt);

      CGLWidget *widget = new CGLWidget;
      widget->show();
      return app.exec();
#else
      fprintf(stderr, "FATAL: FTK not compiled with Qt5.\n");
#endif
    }
  } else {
    fprintf(stderr, "failed to open file %s\n", input_filename.c_str());
  }

  return 0;
}
