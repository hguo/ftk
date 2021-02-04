#ifndef _FTK_XGC_3D_FF_MESH_HH
#define _FTK_XGC_3D_FF_MESH_HH

#include <ftk/config.hh>
#include <ftk/mesh/simplicial_xgc_3d_mesh.hh>
#include <ftk/io/util.hh>

namespace ftk {

/// legacy code: a magnetic field folowing mesh with simplicial elements
template <typename I=int, typename F=double>
struct simplicial_xgc_3dff_mesh : public simplicial_xgc_3d_mesh<I, F> {
  simplicial_xgc_3dff_mesh(std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> m2_, int nphi, int iphi=1, int vphi=1) 
    : simplicial_xgc_3d_mesh<I, F>(m2_, nphi, iphi, vphi) {}

  void initialize_ff_mesh() {}; // TODO
  void initialize_ff_mesh(const std::string& filename);
  void read_ff_mesh(const std::string& filename) { um = simplicial_unstructured_3d_mesh<I, F>::from_file(filename); }
  void write_ff_mesh(const std::string& filename) const { um->to_file(filename); }

public: // reimpl
  size_t n(int d) const;

  void get_simplex(int d, I i, I verts[]) const;
  bool find_simplex(int d, const I v[], I &i) const { return false; } // TODO

protected:
  std::shared_ptr<simplicial_unstructured_3d_mesh<I, F>> um; // the "unit" 3D mesh between poloidal planes
};

///////////////////

template <typename I, typename F>
void simplicial_xgc_3dff_mesh<I, F>::initialize_ff_mesh( const std::string& filename )
{
  fprintf(stderr, "initialize ff mesh from file %s\n", filename.c_str());
  if (file_exists(filename))
    read_ff_mesh(filename);
  else {
    initialize_ff_mesh();
    write_ff_mesh(filename);
  }
}

template <typename I, typename F>
size_t simplicial_xgc_3dff_mesh<I, F>::n(int d) const {
  if (d == 0) return this->m2->n(0) * this->np();
  else if (d == 3) return um->n(3) * this->np();
  else return 0; // TODO: 1- and 2-simplices
}

template <typename I, typename F>
void simplicial_xgc_3dff_mesh<I, F>::get_simplex(int d, I i, I verts[]) const {
  if (d == 3) {
    I j = i % um->n(3);
    I t = i / um->n(3);
    um->get_simplex(3, j, verts);

    for (int i = 0; i < d+1; i ++) 
      verts[i] = (t * this->m2->n(0) + verts[i]) % n(0);

    // um->get_simplex(3, this->transform(d, i), verts);
    // for (int i = 0; i < d+1; i ++)
    //   verts[i] = this->transform(0, verts[i]);
    std::sort(verts, verts+d+1);
    // fprintf(stderr, "j=%d, t=%d, verts=%d, %d, %d, %d\n", j, t, verts[0], verts[1], verts[2], verts[3]);
  } 
  // TODO: 1- and 2-simplicies
}

}

#endif
