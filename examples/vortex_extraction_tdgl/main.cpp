#include <ftk/ftk_config.hh>
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

std::map<hypermesh::regular_simplex_mesh_element, punctured_face_t> punctures;

template <typename T>
inline static T fmod1(T x, T y)
{
  T z = fmod(x, y);
  if (z<0) z += y;
  return z;
}

template <typename T>
inline static T mod2pi(T x)
{
  T y = fmod(x, 2*M_PI); 
  if (y<0) y+= 2*M_PI;
  return y; 
}

template <typename T>
inline static T mod2pi1(T x)
{
  return mod2pi(x + M_PI) - M_PI;
}

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

void evalA(const GLHeader &h, const float X[3], float A[3]) 
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
        evalA(hdr, X[i], A[i]);
      }

      float delta[3], phase_shift = 0;
      for (int i = 0; i < 3; i ++) {
        const int j = (i+1) % 3;
        delta[i] = phi[j] - phi[i];
        float li = line_integral(X[i], X[j], A[i], A[j]);
        delta[i] = mod2pi1(delta[i] - li);
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
      fprintf(stderr, "mu={%f, %f, %f}, pos={%f, %f, %f}\n", mu[0], mu[1], mu[2], pos[0], pos[1], pos[2]);

      punctured_face_t puncture;
      puncture.x[0] = pos[0];
      puncture.x[1] = pos[1];
      puncture.x[2] = pos[2];

      {
        std::lock_guard<std::mutex> guard(mutex);
        punctures[f] = puncture;
      }
  });
}

int main(int argc, char **argv)
{
  bool succ = load_data(argv[1]);
  if (succ)
    extract_vortices();

#if FTK_HAVE_QT5
    QApplication app(argc, argv);
    QGLFormat fmt = QGLFormat::defaultFormat();
    fmt.setSampleBuffers(true);
    fmt.setSamples(16);
    QGLFormat::setDefaultFormat(fmt);

    CGLWidget *widget = new CGLWidget;
    widget->show();
    return app.exec();
#endif

  return 0;
}
