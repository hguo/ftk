#ifndef FTK_XGC_3D_MESH_HH
#define FTK_XGC_3D_MESH_HH

#include <ftk/ftk_config.hh>
#include <ftk/mesh/simplicial_unstructured_3d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_xgc_3d_mesh : public simplicial_unstructured_3d_mesh<I, F> {
  simplicial_xgc_3d_mesh(std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> m2_, int nphi, int iphi=1, int vphi=1);
  // simplicial_xgc_3d_mesh(const std::string& mesh_filename);

  std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> get_m2() const {return m2;}

  size_t n(int d) const;
  size_t np() const {return nphi * iphi * vphi;} // number of poloidal planes, incl. virtual planes defined by vphi

  void set_nphi_iphi(int n, int i) {nphi = n; iphi = i;}
  void set_vphi(int v) { vphi = v; }

  int get_nphi() const {return nphi;}
  int get_iphi() const {return iphi;}
  int get_vphi() const {return vphi;}

  bool is_poloidal(int p) const { return p % vphi == 0; }
  bool is_poloidal(int d, I i) const { return m3->is_ordinal(d, i); }

public: 
  void element_for(int d, std::function<void(I)> f) {} // TODO
  
public:
  void get_simplex(int d, I i, I verts[]) const;
  bool find_simplex(int d, const I v[], I& i) const;
  
  void get_coords_rzp(I i, F coords[]) const { return m3->get_coords(i, coords); }
  void get_coords_xyz(I i, F coords[]) const;
 
public:
  std::set<I> sides(int d, I i) const;
  std::set<I> side_of(int d, I i) const;

  I transform(int d, I i) const;

private: // backend meshes
  int nphi, iphi, vphi;
  std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> m2;
  std::shared_ptr<simplicial_unstructured_extruded_2d_mesh<I, F>> m3;
};
///////
//

template <typename I, typename F>
simplicial_xgc_3d_mesh<I, F>::simplicial_xgc_3d_mesh(
    std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> m2_, 
    int nphi_, int iphi_, int vphi_) :
  m2(m2_),
  m3(new simplicial_unstructured_extruded_2d_mesh<I, F>(*m2_)),
  nphi(nphi_), iphi(iphi_), vphi(vphi_)
{

}

template <typename I, typename F>
size_t simplicial_xgc_3d_mesh<I, F>::n(int d) const
{
  return m3->n(d) * np();
}

template <typename I, typename F>
void simplicial_xgc_3d_mesh<I, F>::get_simplex(int d, I i, I verts[]) const
{
  // fprintf(stderr, "trying to get %d-simplex %d(%d)\n", d, i, transform(d, i));

  m3->get_simplex(d, transform(d, i), verts);
  for (int i = 0; i < d+1; i ++)
    verts[i] = transform(0, verts[i]);
  std::sort(verts, verts+d+1);
}

template <typename I, typename F>
bool simplicial_xgc_3d_mesh<I, F>::find_simplex(int d, const I v[], I& i) const
{
  I verts[d+1];
  int t[d+1];

  bool has0 = false, has1 = false;
  for (int i = 0; i < d+1; i ++) {
    verts[i] = transform(0, v[i]);
    t[i] = v[i] / m2->n(0); // wrong: m3->flat_vertex_time(v[i]);
    if (t[i] == 0) has0 = true;
    if (t[i] == np()-1) has1 = true;
  }
  bool crossed = has0 && has1;

  if (crossed)
  for (int i = 0; i < d+1; i ++)
    if (t[i] == 0) verts[i] += n(0);

  // fprintf(stderr, "xgc_mesh: trying to find %d-d vert %d, %d, %d, %d, t=%d, %d, %d, %d, crossed=%d\n", 
  //     d, verts[0], verts[1], verts[2], verts[3],
  //     t[0], t[1], t[2], t[3], crossed);

  std::sort(verts, verts+d+1);
  bool succ = m3->find_simplex(d, verts, i);
  assert(succ);
  // fprintf(stderr, "succ=%d, i=%d, %d\n", succ, i, transform(d, i));

  for (int i = 0; i < d+1; i ++)
    verts[i] = transform(0, verts[i]);

  i = transform(d, i);
  return succ;
}

template <typename I, typename F>
I simplicial_xgc_3d_mesh<I, F>::transform(int d, I i) const
{
  while (i < 0) i += n(d);
  return i % n(d);
}

template <typename I, typename F>
std::set<I> simplicial_xgc_3d_mesh<I, F>::sides(int d, I i) const 
{
  std::set<I> results;
  for (auto j : m3->sides(d, transform(d, i))) 
    results.insert(transform(d-1, j));
    // results.insert(j); 
  return results;
}

template <typename I, typename F>
std::set<I> simplicial_xgc_3d_mesh<I, F>::side_of(int d, I i) const 
{
  std::set<I> results;
  for (auto j : m3->side_of(d, transform(d, i))) 
    results.insert(transform(d+1, j));
    // results.insert(j % n(d+1));
    // results.insert(j); 
  return results;
}

}

#endif
