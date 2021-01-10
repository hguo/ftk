#ifndef _FTK_POINT_LOCATOR_2D_HH
#define _FTK_POINT_LOCATOR_2D_HH

#include <ftk/ftk_config.hh>
#include <ftk/mesh/aabb.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct point_locator_2d {
  point_locator_2d(std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> m) : m2(m) {}

  virtual void initialize() = 0;
  virtual I locate(const F x[]) const = 0;

protected:
  std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> m2;
};

}

#endif
