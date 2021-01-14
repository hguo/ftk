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

  virtual size_t n(int d) const;
  size_t np() const {return nphi * iphi * vphi;} // number of poloidal planes, incl. virtual planes defined by vphi

  void set_nphi_iphi(int n, int i) {nphi = n; iphi = i;}
  void set_vphi(int v) { vphi = v; }

  int get_nphi() const {return nphi;}
  int get_iphi() const {return iphi;}
  int get_vphi() const {return vphi;}

  bool is_poloidal(int p) const { return p % vphi == 0; }
  bool is_poloidal(int d, I i) const { return m3->is_ordinal(d, i); }

  std::set<I> get_vertex_edge_vertex(I i) const;
  std::set<I> get_vertex_edge_vertex_nextnodes(I i) const;

public: 
  void element_for(int d, std::function<void(I)> f) {} // TODO
  
public:
  virtual void get_simplex(int d, I i, I verts[]) const;
  virtual bool find_simplex(int d, const I v[], I& i) const;
  
  void get_coords_rzp(I i, F coords[]) const { return m3->get_coords(i, coords); }
  void get_coords_xyz(I i, F coords[]) const;
  void get_coords(I i, F coords[]) const { get_coords_xyz(i, coords); }
 
public:
  virtual std::set<I> sides(int d, I i) const;
  virtual std::set<I> side_of(int d, I i) const;

  I transform(int d, I i) const;

public: // vtk
  void to_vtu_slices_file(const std::string& filename) const;
  void scalar_to_vtu_slices_file(const std::string& filename, const std::string& varname, const ndarray<F>& data) const;
  void to_vtu_solid_file(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu_slices() const;
  vtkSmartPointer<vtkUnstructuredGrid> scalar_to_vtu_slices(const std::string& varname, const ndarray<F>& data) const;
  virtual vtkSmartPointer<vtkUnstructuredGrid> to_vtu_solid() const { return this->to_vtu(); }
#endif

public: // smoothing
  bool has_smoothing_kernel() const { return m2->has_smoothing_kernel(); }
  ndarray<F> smooth_scalar(const ndarray<F>& f) const;

public:
  struct interpolant_t { // poloidal interpolants
    I tri0[3], tri1[3];
    F mu0[3], mu1[3];
  };
 
  void initialize_interpolants();
  ndarray<F> interpolate(const ndarray<F>& scalar) const; // interpolate virtual planes
  F interpolate(const ndarray<F>& scalar, I i); // interpolate scalar value
  void interpolate(const ndarray<F>& scalar, const ndarray<F>& grad, const ndarray<F>& jacobian, 
      I i, F f[], F g[2], F j[2][2]) const;

  const interpolant_t& get_interpolant(int v /*virtual plane id*/, I i/*vertex id*/) const { return interpolants[v][i]; }
  
protected:
  std::vector<std::vector<interpolant_t>> interpolants;

protected: // backend meshes
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

template <typename I, typename F>
std::set<I> simplicial_xgc_3d_mesh<I, F>::get_vertex_edge_vertex(I i) const
{
  std::set<I> verts;
  for (const auto v : m3->get_vertex_edge_vertex(i))
    verts.insert( transform(0, v) );
  return verts;
}

template <typename I, typename F>
std::set<I> simplicial_xgc_3d_mesh<I, F>::get_vertex_edge_vertex_nextnodes(I i) const
{
  const int m2n0 = m2->n(0);
  const int t = i / m2n0;
  const I i0 = i % m2n0;

  std::set<I> set2 = m2->get_vertex_edge_vertex(i0);
  std::set<I> set3;
  for (const auto j : set2)
    set3.insert( j + t * m2n0 );
  set3.insert( (m2->nextnode(i0) + (t+1) * m2n0) % n(0) ); // assuming vphi=1

  return set3;
}

template <typename I, typename F>
void simplicial_xgc_3d_mesh<I, F>::get_coords_xyz(I i, F coords[]) const
{
  F rzp[3];
  get_coords_rzp(i, rzp);

  const F phi = rzp[2] * 2 * M_PI / np();
  coords[0] = rzp[0] * cos(phi);
  coords[1] = rzp[0] * sin(phi);
  coords[2] = rzp[1];
}

#if FTK_HAVE_VTK
template <typename I, typename F>
void simplicial_xgc_3d_mesh<I, F>::to_vtu_slices_file(const std::string& filename) const
{
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( to_vtu_slices() );
  writer->Write();
}

template <typename I, typename F>
void simplicial_xgc_3d_mesh<I, F>::
scalar_to_vtu_slices_file(const std::string& filename, const std::string& varname, const ndarray<F>& scalar) const
{
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( scalar_to_vtu_slices(varname, scalar) );
  writer->Write();
}

template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_xgc_3d_mesh<I, F>::
to_vtu_slices() const
{
  const int m2n0 = m2->n(0);
  const int m3n0 = n(0);

  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(m3n0);

  const auto &vertex_coords = m2->get_coords();
  const auto &triangles = m2->get_triangles();

  for (int p = 0; p < np(); p ++) {
    const F phi = p * 2 * M_PI / np();
    const vtkIdType offset = p * m2n0;
    for (int i=0; i < m2n0; i ++) 
      pts->SetPoint(offset + i, 
          vertex_coords[i*2] * cos(phi),
          vertex_coords[i*2] * sin(phi), 
          vertex_coords[i*2+1]);
  }
 
  for (int p = 0; p < np(); p ++) {
    const vtkIdType offset = p * m2n0;
    for (int i=0; i<m2->n(2); i ++) {
      vtkIdType ids[3] = {
        offset + triangles[i*3], 
        offset + triangles[i*3+1], 
        offset + triangles[i*3+2]
      };
      grid->InsertNextCell(VTK_TRIANGLE, 3, ids);
    }
  }

  grid->SetPoints(pts);
  return grid;
}

template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_xgc_3d_mesh<I, F>::
scalar_to_vtu_slices(const std::string& varname, const ndarray<F>& scalar) const 
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = to_vtu_slices();
  vtkSmartPointer<vtkDataArray> array = vphi == 1 ? 
    scalar.to_vtk_data_array(varname) : 
    interpolate(scalar).to_vtk_data_array(varname);
  
  grid->GetPointData()->AddArray(array);
  grid->GetPointData()->SetActiveScalars(varname.c_str());

  // grid->PrintSelf(std::cerr, vtkIndent(2));
  return grid;
}
#endif // HAVE_FTK

template <typename I, typename F>
void simplicial_xgc_3d_mesh<I, F>::initialize_interpolants()
{
  if (vphi == 1) return;

  const I m2n0 = m2->n(0);
  const F dphi = 2 * M_PI / (nphi * iphi);
 
  fprintf(stderr, "initializing interpolants...\n");
  interpolants.resize(vphi);
  for (int v = 1; v < vphi; v ++) {
    interpolants[v].resize(m2n0);
    const F beta = F(v) / vphi, alpha = F(1) - beta;
    this->parallel_for(m2n0, [&](int i) {
    // for (I i = 0; i < m2n0; i ++) {
      interpolant_t &l = interpolants[v][i];
      F rzp0[3], rzp1[3];
        
      // backward integration
      m2->get_coords(i, rzp0);
      rzp0[2] = beta * dphi;
      m2->magnetic_map(rzp0, 0);
      I tid0 = m2->locate(rzp0, l.mu0);
      if (tid0 < 0) { // invalid triangle
        l.tri0[0] = l.tri0[1] = l.tri0[2] = m2->nearest(rzp0);
        l.mu0[0] = l.mu0[1] = l.mu0[2] = F(1) / 3;
      } else {
        m2->get_simplex(2, tid0, l.tri0);
      }
     
      // forward integration
      m2->get_coords(i, rzp1);
      rzp1[2] = beta * dphi;
      m2->magnetic_map(rzp1, dphi);
      I tid1 = m2->locate(rzp1, l.mu1);
      if (tid1 < 0) { // invalid triangle
        l.tri1[0] = l.tri1[1] = l.tri1[2] = m2->nearest(rzp1);
        l.mu1[0] = l.mu1[1] = l.mu1[2] = F(1) / 3;
      } else {
        m2->get_simplex(2, tid1, l.tri1);
      }
    });
  }
  fprintf(stderr, "interpolants initialized.\n");
}
  
template <typename I, typename F>
void simplicial_xgc_3d_mesh<I, F>::interpolate(
    const ndarray<F>& scalar, const ndarray<F>& grad, const ndarray<F>& jacobian, 
    I i, F f[], F g[2], F j[2][2]) const
{
  const int m2n0 = m2->n(0);
  const int p = i / m2n0; // poloidal plane id

  if (p % vphi == 0) { // non-virtual plane
    const int p0 = p / vphi;
    const int k = i % m2n0;

    const auto idx = m2n0 * p0 + k;
    f[0] = scalar[idx];
    g[0] = grad[idx*2];
    g[1] = grad[idx*2+1];
    j[0][0] = jacobian[idx*4];
    j[0][1] = jacobian[idx*4+1];
    j[1][0] = jacobian[idx*4+2];
    j[1][1] = jacobian[idx*4+3];
  } else { // virtual plane
    const int p0 = p / vphi, p1 = (p0 + 1) % nphi;
    const F beta = F(p) / vphi - p0, alpha = F(1) - beta;
    const interpolant_t& l = interpolants[p % vphi][i % m2n0];

    // init
    f[0] = 0;
    for (int k0 = 0; k0 < 2; k0 ++) {
      g[k0] = 0;
      for (int k1 = 0; k1 < 2; k1 ++) 
        j[k0][k1] = 0;
    }

    // add values
    for (int k = 0; k < 3; k ++) {
      const auto idx0 = m2n0 * p0 + l.tri0[k], 
                 idx1 = m2n0 * p1 + l.tri1[k];

      f[0] += alpha * l.mu0[k] * scalar[idx0] + beta * l.mu1[k] * scalar[idx1];
      for (int k0 = 0; k0 < 2; k0 ++) {
        g[k0] += alpha * l.mu0[k] * grad[idx0*2+k0] + beta * l.mu1[k] * grad[idx1*2+k0];
        for (int k1 = 0; k1 < 2; k1 ++) 
          j[k0][k1] += alpha * l.mu0[k] * jacobian[idx0*4+k0*2+k1] + beta * l.mu1[k] * jacobian[idx1*4+k0*2+k1];
      }
    }
  }
}

template <typename I, typename F>
F simplicial_xgc_3d_mesh<I, F>::interpolate(const ndarray<F>& scalar, I i) // interpolate scalar value
{
  const int m2n0 = m2->n(0);
  const int p = i / m2n0; // poloidal plane id

  if (p % vphi == 0) { // non-virtual plane
    const int p0 = p / vphi;
    const int j = i % m2n0;
    return scalar[m2n0 * p0 + j];
  } else { // virtual plane
    const int p0 = p / vphi, p1 = (p0 + 1) % nphi;
    const F beta = F(p) / vphi - p0, alpha = F(1) - beta;
    const interpolant_t& l = interpolants[p % vphi][i % m2n0];

    F f0 = 0, f1 = 0;
    for (int k = 0; k < 3; k ++) {
      f0 += l.mu0[k] * scalar[m2n0 * p0 + l.tri0[k]];
      f1 += l.mu1[k] * scalar[m2n0 * p1 + l.tri1[k]];
    }

    return alpha * f0 + beta * f1;
  }
}

template <typename I, typename F>
ndarray<F> simplicial_xgc_3d_mesh<I, F>::smooth_scalar(const ndarray<F>& scalar) const
{
  ndarray<double> result;
  result.reshape(scalar);
  
  for (size_t i = 0; i < nphi; i ++) {
    auto slice = scalar.slice_time(i);
    auto f = m2->smooth_scalar(slice);
    for (size_t k = 0; k < m2->n(0); k ++)
      result(k, i) = f(k);
  }

  return result;
}

template <typename I, typename F>
ndarray<F> simplicial_xgc_3d_mesh<I, F>::interpolate(const ndarray<F>& scalar) const
{
  const int m2n0 = m2->n(0), 
            nvphi = nphi * vphi;
  const F dphi = 2 * M_PI / (nphi * iphi);

  ndarray<F> f;
  f.reshape(m2n0, nvphi);

  for (int i = 0; i < nvphi; i ++) {
    if (i % vphi == 0) { // non-virtual plane, just copy the values
      const int p = i / vphi;
      for (int j = 0; j < m2n0; j ++)
        f[m2n0 * i + j] = scalar[m2n0 * p + j];
    } else {
      const int p0 = i / vphi, p1 = (p0 + 1) % nphi;
      const F beta = F(i) / vphi - p0, alpha = F(1) - beta;
      // fprintf(stderr, "p0=%d, p1=%d, alpha=%f, beta=%f\n", p0, p1, alpha, beta);
      const std::vector<interpolant_t>& ls = interpolants[i % vphi];

      for (int j = 0; j < m2n0; j ++) {
        const interpolant_t& l = ls[j];
        // fprintf(stderr, "node %d, tri0=%d, %d, %d, mu0=%f, %f, %f, tri1=%d, %d, %d, mu1=%f, %f, %f\n", 
        //     j, l.tri0[0], l.tri0[1], l.tri0[2], l.mu0[0], l.mu0[1], l.mu0[2],
        //     l.tri1[0], l.tri1[1], l.tri1[2], l.mu1[0], l.mu1[1], l.mu1[2]);
        
        F f0 = 0, f1 = 0;
        for (int k = 0; k < 3; k ++) {
          f0 += l.mu0[k] * scalar[m2n0 * p0 + l.tri0[k]];
          f1 += l.mu1[k] * scalar[m2n0 * p1 + l.tri1[k]];
        }
        
        f[m2n0 * i + j] = alpha * f0 + beta * f1;
#if 0
        F rzp0[3], rzp1[3];
        F mu0[3], mu1[3];
        I tri0[3], tri1[3];

        // backward integration
        m2->get_coords(j, rzp0);
        rzp0[2] = alpha * dphi;
        m2->magnetic_map(rzp0, 0);
        I tid0 = m2->locate(rzp0, mu0);
        if (tid0 < 0) continue; // invalid triangle
       
        // forward integration
        m2->get_coords(j, rzp1);
        rzp1[2] = alpha * dphi;
        m2->magnetic_map(rzp1, dphi);
        I tid1 = m2->locate(rzp1, mu1);
        if (tid1 < 0) continue; // invalid triangle

        F f0 = 0, f1 = 0;
        m2->get_simplex(2, tid0, tri0);
        for (int k = 0; k < 3; k ++)
          f0 += mu0[k] * scalar[m2n0 * p0 + tri0[k]];
        
        m2->get_simplex(2, tid1, tri1);
        for (int k = 0; k < 3; k ++)
          f1 += mu1[k] * scalar[m2n0 * p1 + tri1[k]];
        
        f[m2n0 * i + j] = alpha * f0 + beta * f1;
#endif

#if 0
        const int next = m2->nextnode(j);
        F xcur[2], xnext[2], x[2];
        m2->get_coords(j, xcur);
        m2->get_coords(next, xnext);

        x[0] = alpha * xcur[0] + beta * xnext[0];
        x[1] = alpha * xcur[1] + beta * xnext[1];
        const int k = m2->nearest(x);

        F f0 = scalar[m2n0 * p0 + j], 
          f1 = scalar[m2n0 * p1 + next];

        f[m2n0 * i + k] = alpha * f0 + beta * f1;
#endif
      }
    }
  }
  return f;
}

}

#endif
