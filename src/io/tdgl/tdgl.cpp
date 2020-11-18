#include <ftk/io/tdgl.hh>
#include "GLGPU_IO_Helper.h"

namespace ftk {

bool tdgl_reader::read()
{
  float *rho = NULL, *phi = NULL, *re = NULL, *im = NULL, *Jx = NULL, *Jy = NULL, *Jz = NULL;
  bool succ = GLGPU_IO_Helper_ReadBDAT(filename, meta, &rho, &phi, &re, &im, &Jx, &Jy, &Jz, !read_field_data, derive_J);
  if (!succ)
    succ = GLGPU_IO_Helper_ReadLegacy(filename, meta, &rho, &phi, &re, &im, &Jx, &Jy, &Jz, !read_field_data, derive_J);

  if (!succ) return false;

  if (read_field_data) {
    std::vector<size_t> shape;
    for (int i = 0; i < 3; i ++)
      shape.push_back(meta.dims[i]);
    
    this->rho.from_array(rho, shape);
    this->phi.from_array(phi, shape);
    this->re.from_array(re, shape);
    this->im.from_array(im, shape);

    // TODO: J
  }

  if (rho) free(rho);
  if (phi) free(phi);
  if (re) free(re);
  if (im) free(im);
  if (Jx) free(Jx);
  if (Jy) free(Jy);
  if (Jz) free(Jz);

  return true;
}

}
