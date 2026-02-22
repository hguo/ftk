#ifndef _FTK_POINT_LOCATOR_3D_HH
#define _FTK_POINT_LOCATOR_3D_HH

#include <ftk/config.hh>
#include <ftk/mesh/aabb.hh>
#include <ftk/mesh/simplicial_unstructured_3d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct point_locator_3d {
  point_locator_3d(std::shared_ptr<simplicial_unstructured_3d_mesh<I, F>> m) : m3(m) {}
  virtual ~point_locator_3d() {}

  virtual void initialize() = 0;
  virtual I locate(const F x[], F mu[]) const = 0;
  
  I locate(const F x[]) const { F mu[4]; return locate(x, mu);  }

protected:
  std::shared_ptr<simplicial_unstructured_3d_mesh<I, F>> m3;
};

}

#endif
