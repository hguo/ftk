#ifndef _FTK_POINT_LOCATOR_2D_HH
#define _FTK_POINT_LOCATOR_2D_HH

#include <ftk/config.hh>
#include <ftk/mesh/aabb.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct point_locator_2d {
  point_locator_2d(const simplicial_unstructured_2d_mesh<I, F>& m) : m2(m) {}
  virtual ~point_locator_2d() {}

  virtual void initialize() = 0;
  virtual I locate(const F x[], F mu[]) const = 0;
  
  I locate(const F x[]) const { F mu[3]; return locate(x, mu);  }

protected:
  const simplicial_unstructured_2d_mesh<I, F>& m2;
};

}

#endif
