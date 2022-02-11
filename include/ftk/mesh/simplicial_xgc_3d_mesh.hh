#ifndef FTK_XGC_3D_MESH_HH
#define FTK_XGC_3D_MESH_HH

#include <ftk/config.hh>
#include <ftk/mesh/simplicial_xgc_2d_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_3d_mesh.hh>
#include <ftk/numeric/xgc_interpolant.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_xgc_3d_mesh : public simplicial_unstructured_3d_mesh<I, F> {
  simplicial_xgc_3d_mesh(std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> m2_); // nphi and iphi must be specified before using
  simplicial_xgc_3d_mesh(std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> m2_, int nphi, int iphi=1, int vphi=1);
  // simplicial_xgc_3d_mesh(const std::string& mesh_filename);

  std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> get_m2() const {return m2;}

  virtual size_t n(int d) const;
  size_t np() const {return nphi * iphi * vphi;} // number of poloidal planes, incl. virtual planes defined by vphi

  bool probe_nphi_iphi(const std::string& filename);
  void set_nphi_iphi(int n, int i) {nphi = n; iphi = i;}
  void set_vphi(int v) { vphi = v; }

  int get_nphi() const {return nphi;}
  int get_iphi() const {return iphi;}
  int get_vphi() const {return vphi;}

  bool is_poloidal(int p) const { return p % vphi == 0; }
  bool is_poloidal(int d, I i) const { return m3->is_ordinal(d, i); }
  bool is_actual_poloidal(int d, I i) const { 
    if (m3->is_ordinal(d, i)) {
      int v[4];
      get_simplex(d, i, v);
      const int t = m3->flat_vertex_time(v[0]);
      // fprintf(stderr, "t=%d, vphi=%d, t_mod_vphi=%d\n", t, vphi, t%vphi);
      return t % vphi == 0;
    } else return false;
  }

  int get_poloidal(int d, I i) const { return i / m3->n(d); }

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

public: // io
  ndarray<F> derive_turbulence(
      const ndarray<F>& dpot, 
      const ndarray<F>& pot0,
      const ndarray<F>& potm0,
      const ndarray<F>& eden, 
      const ndarray<F>& psi_mks,  // 80
      const ndarray<F>& e_gc_density_avg,
      const ndarray<F>& e_perp_temperature_avg,
      const ndarray<F>& e_parallel_mean_en_avg) const;

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
  void smooth_scalar_gradient_jacobian(
      const ndarray<F>& f, 
      ndarray<F>& fs, // smoothed scalar field
      ndarray<F>& g,  // smoothed gradient field
      ndarray<F>& j   // smoothed jacobian field
  ) const; 

public:
  void initialize_interpolants_cached(); // if the default interpolant file exists, load interpolants; otherwise initialize and write interpolants to the default file
  void initialize_interpolants();
  ndarray<F> interpolate(const ndarray<F>& scalar) const; // interpolate virtual planes
 
  std::string default_interpolant_filename() const;

  F interpolate(const ndarray<F>& scalar, I i); // interpolate scalar value at vertex i
  F interpolate(const ndarray<F>& scalar, F r, F z, I pi); // interpolate scalar value at rz and pi

  void interpolate(const ndarray<F>& scalar, const ndarray<F>& grad, const ndarray<F>& jacobian, 
      I i, F f[], F g[2], F j[2][2]) const;
  
  void interpolate_central_difference(const ndarray<F>& scalar, const ndarray<F>& grad, const ndarray<F>& jacobian, 
      I i, F f[], F g[2], F j[2][2], I delta) const;

  const xgc_interpolant_t<I, F>& get_interpolant(int v /*virtual plane id*/, I i/*vertex id*/) const { return interpolants[v][i]; }
 
  void write_interpolants(const std::string& filename) const { diy::serializeToFile(interpolants, filename); }
  void read_interpolants(const std::string& filename) { diy::unserializeFromFile(filename, interpolants); }
  const std::vector<std::vector<xgc_interpolant_t<I, F>>>& get_interpolants() const { return interpolants; }

protected:
  std::vector<std::vector<xgc_interpolant_t<I, F>>> interpolants;

public: // rotation
  void initialize_rotational_interpolants();
  ndarray<F> rotate(const ndarray<F>& scalar) const; // rotate poloidal planes to match the first plane

protected:
  std::vector<std::vector<std::tuple<I, F, F/* triangleId, mu0, mu1 */>>> rotational_interpolants;

protected: // backend meshes
  int nphi, iphi, vphi;
  std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> m2;
  std::shared_ptr<simplicial_unstructured_extruded_2d_mesh<I, F>> m3;
};
///////
//

template <typename I, typename F>
simplicial_xgc_3d_mesh<I, F>::simplicial_xgc_3d_mesh(
    std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> m2_) :
  m2(m2_),
  m3(new simplicial_unstructured_extruded_2d_mesh<I, F>(*m2_))
{
}

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
void simplicial_xgc_3d_mesh<I, F>::initialize_interpolants_cached()
{
  if (vphi == 1) return; // no interpolants needed

  const auto f = default_interpolant_filename();
  if (file_exists(f))
    read_interpolants(f);
  else {
    initialize_interpolants();
    write_interpolants(f);
  }
}

template <typename I, typename F>
std::string simplicial_xgc_3d_mesh<I, F>::default_interpolant_filename() const
{
  // hash of conn (2d), nphi, iphi, vphi
  const int arr[3] = {nphi, iphi, vphi};

  unsigned int h0 = m2->get_triangles().hash();
  unsigned int h1 = murmurhash2(arr, 3*sizeof(int), h0);

  std::stringstream ss;
  ss << "xgc.interpolants." 
     << std::hex << h1;
  return ss.str();
}

template <typename I, typename F>
bool simplicial_xgc_3d_mesh<I, F>::probe_nphi_iphi(const std::string& filename)
{
  // determine nphi and iphi
  const auto ext = file_extension(filename);
  const auto array_nphi = ndarray<int>::from_file(filename, "nphi");
  const auto array_iphi = ndarray<int>::from_file(filename, "iphi");

  nphi = array_nphi[0];
  iphi = std::max(1, array_iphi[0]);
  return true;
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
      xgc_interpolant_t<I, F> &l = interpolants[v][i];
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
    const xgc_interpolant_t<I, F>& l = interpolants[p % vphi][i % m2n0];

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
F simplicial_xgc_3d_mesh<I, F>::interpolate(
    const ndarray<F>& scalar, I i) // interpolate scalar value
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
    const xgc_interpolant_t<I, F>& l = interpolants[p % vphi][i % m2n0];

    F f0 = 0, f1 = 0;
    for (int k = 0; k < 3; k ++) {
      f0 += l.mu0[k] * scalar[m2n0 * p0 + l.tri0[k]];
      f1 += l.mu1[k] * scalar[m2n0 * p1 + l.tri1[k]];
    }

    return alpha * f0 + beta * f1;
  }
}

#if 0
template <typename I, typename F>
F simplicial_xgc_3d_mesh<I, F>::interpolate(
    const ndarray<F>& scalar, F r, F z, I pi) // interpolate scalar value for arbitrary rzp
{
  const int m2n0 = m2->n(0);
  const int p = pi / m2n0; // poloidal plane id

  const int p0 = p / vphi, p1 = (p0 + 1) % nphi;
  const F beta = F(p) / vphi - p0, alpha = F(1) - beta;
  const xgc_interpolant_t<I, F>& l = interpolants[pi % vphi][i % m2n0];

  if (p % vphi == 0) { // non-virtual plane
    const int p0 = p / vphi;
    const int j = i % m2n0;
    return scalar[m2n0 * p0 + j];
  } else { // virtual plane
    const int p0 = p / vphi, p1 = (p0 + 1) % nphi;
    const F beta = F(p) / vphi - p0, alpha = F(1) - beta;
    const xgc_interpolant_t<I, F>& l = interpolants[p % vphi][i % m2n0];

    F f0 = 0, f1 = 0;
    for (int k = 0; k < 3; k ++) {
      f0 += l.mu0[k] * scalar[m2n0 * p0 + l.tri0[k]];
      f1 += l.mu1[k] * scalar[m2n0 * p1 + l.tri1[k]];
    }

    return alpha * f0 + beta * f1;
  }
}

template <typename I, typename F>
void simplicial_xgc_3d_mesh<I, F>::interpolate_central_difference(
    const ndarray<F>& scalar, 
    const ndarray<F>& grad, 
    const ndarray<F>& jacobian, 
    I i, F val[], F g[2], F j[2][2], I delta) const
{
  F coords[3];
  get_coords_rzp(i, coords);
  const F r = coords[0], z = coords[1], p = coords[2];

  const F f[5][5] = {
    { 0, 0, interpolate(scalar, r, z+2*delta, p), 0, 0 }, 
    { 0, interpolate(scalar, r-delta, z+delta, p), interpolate(scalar, r, z+delta, p), interpolate(scalar, r+delta, z+delta, p), 0 },
    { interpolate(scalar, r-2*delta, z, p), interpolate(scalar, r-delta, z, p), interpolate(scalar, i), interpolate(scalar, r+delta, z, p), interpolate(scalar, r+2*delta, z, p) }, 
    { 0, interpolate(scalar, r-delta, z-delta, p), interpolate(scalar, r, z-delta, p), interpolate(scalar, r+delta, z-delta, p), 0 },
    { 0, 0, interpolate(scalar, r, z-2*delta, p), 0, 0 }
  };

  const F df[3][3][2] = {
    { { 0, 0 },                                 { f[1][3] - f[1][2], f[0][2] - f[2][2] }, { 0, 0 } }, 
    { { f[2][2] - f[2][0], f[1][1] - f[3][1] }, { f[2][3] - f[2][1], f[1][2] - f[3][2] }, { f[2][4] - f[2][2], f[1][3] - f[3][3] } }, 
    { { 0, 0 },                                 { f[3][3] - f[3][1], f[2][2] - f[4][2] }, { 0, 0 } }
  };

  j[0][0] = df[1][2][0] - df[1][0][0];
  j[0][1] = df[1][2][1] - df[1][0][1];
  j[1][0] = df[0][1][0] - df[2][1][0];
  j[1][1] = df[0][1][1] - df[2][1][1];

  g[0] = df[2][2][0];
  g[1] = df[2][2][1];

  val[0] = f[2][2];
}
#endif

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
void simplicial_xgc_3d_mesh<I, F>::smooth_scalar_gradient_jacobian(
      const ndarray<F>& scalar, 
      ndarray<F>& S, // smoothed scalar field
      ndarray<F>& G,  // smoothed gradient field
      ndarray<F>& J   // smoothed jacobian field
  ) const
{
  S.reshape(scalar);
  G.reshape(2, scalar.dim(0), scalar.dim(1));
  G.set_multicomponents(1);
  J.reshape(2, 2, scalar.dim(0), scalar.dim(1));
  J.set_multicomponents(2);
  
  for (size_t i = 0; i < nphi; i ++) {
    ndarray<double> f, grad, j;
    auto slice = scalar.slice_time(i);
    m2->smooth_scalar_gradient_jacobian(slice, f, grad, j);
    for (size_t k = 0; k < m2->n(0); k ++) {
      S(k, i) = f(k);
      G(0, k, i) = grad(0, k);
      G(1, k, i) = grad(1, k);
      J(0, 0, k, i) = j(0, 0, k);
      J(1, 0, k, i) = j(1, 0, k);
      J(1, 1, k, i) = j(1, 1, k);
      J(0, 1, k, i) = j(0, 1, k);
    }
  }
}

template <typename I, typename F>
ndarray<F> simplicial_xgc_3d_mesh<I, F>::derive_turbulence(
      const ndarray<F>& dpot, 
      const ndarray<F>& pot0,
      const ndarray<F>& potm0,
      const ndarray<F>& eden, 
      const ndarray<F>& psid, // psi_mks, 80
      const ndarray<F>& dens, // e_gc_density_avg, 80
      const ndarray<F>& temp1, // e_perp_temperature_avg,
      const ndarray<F>& temp2) const // e_parallel_mean_en_avg
{
  // see http://visit.ilight.com/svn/visit/trunk/src/databases/ADIOS/avtXGCFileFormat.C
  auto interpolate = [&](const ndarray<F>& x, const ndarray<F>& y, const ndarray<F>& xi) {
    const int n = x.size(), ni = xi.size();
    ndarray<F> yi(xi.shape());
    for (int i = 0; i < ni; i ++) {
      F val = xi[i];
      if (val < x[0])
        yi[i] = y[0];
      else if (val >= x[n-1])
        yi[i] = y[n-1];
      else {
        for (int j = 0; j < n-1; j ++) {
          if (val >= x[j] && val <= x[j+1]) {
            F dy = y[j+1] - y[j];
            F t = (xi[i] - x[j]) / (x[j+1] - x[j]);
            yi[i] = y[j] + t * dy;
            break;
          }
        }
      }
    }
    return yi;
  };

  const ndarray<F>& psi = m2->get_psifield();

  ndarray<F> temp(temp1.shape());
  for (int i = 0 ; i < temp.size(); i ++)
    temp[i] = 2.0 * (temp1[i] + temp2[i]) / 3.0;

  const int n = psid.size();
  const int ni = psi.size();

  ndarray<F> te = interpolate(psid, temp, psi);
  ndarray<F> de = interpolate(psid, dens, psi);

  const int m2n0 = m2->n(0);
  ndarray<F> mean_eden; 
  mean_eden.reshape(m2n0);
  for (int i = 0; i < nphi; i ++) 
    for (int j = 0; j < m2n0; j ++) 
      mean_eden[j] += eden[i*m2n0+j]; // eden(j, i);
  for (int j = 0; j < m2n0; j ++)
    mean_eden[j] /= nphi;

  ndarray<F> mean_pot;
  mean_pot.reshape(m2n0);
  for (int i = 0; i < nphi; i ++) 
    for (int j = 0; j < m2n0; j ++) 
      mean_pot[j] += dpot[i*m2n0+j]; // eden(j, i);
  for (int j = 0; j < m2n0; j ++)
    mean_pot[j] /= nphi;

  // ndarray<F> arr({size_t(m2n0), size_t(nphi)});
  ndarray<F> arr(eden.shape());
  for (int i = 0; i < nphi; i ++) {
    for (int j = 0; j < m2n0; j ++) {
      const int idx = i*m2n0 + j;
      // F v1 = dpot[idx] - (potm0[j] - pot0[j]);
      F v1 = dpot[idx] - mean_pot[j]; 
      v1 = v1 / te[j];

      F v2 = eden[idx] - mean_eden[j];
      v2 = v2 / de[j];

      // arr[idx] = v1 + v2;
      arr[idx] = v1;
    }
  }

  // check mean and std
#if 0
  F mean(0), var(0);
  for (int i = 0; i < arr.size(); i ++)
    mean += arr[i];
  mean /= arr.size();

  for (int i = 0; i < arr.size(); i ++) {
    F d = arr[i] - mean;
    var += d * d;
  }
  F s = std::sqrt(var / arr.size());
  fprintf(stderr, "mean=%f, std=%f\n", mean, s);
#endif

  return arr;
}

template <typename I, typename F>
ndarray<F> simplicial_xgc_3d_mesh<I, F>::interpolate(
    const ndarray<F>& scalar) const
{
  const size_t m2n0 = m2->n(0), 
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
      const std::vector<xgc_interpolant_t<I, F>>& ls = interpolants[i % vphi];

      for (int j = 0; j < m2n0; j ++) {
        const xgc_interpolant_t<I, F>& l = ls[j];
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

template <typename I, typename F>
void simplicial_xgc_3d_mesh<I, F>::initialize_rotational_interpolants()
{
  const I m2n0 = m2->n(0);
  const F dphi = 2 * M_PI / (nphi * iphi);
  const int nsteps_per_dphi = 128;

  fprintf(stderr, "initialize rotational interpolants\n");
  rotational_interpolants.resize(nphi);
  for (int i = 1; i < nphi; i ++) {
    rotational_interpolants[i].resize(m2n0);
    this->parallel_for(m2n0, [&](int k) {
      auto &l = rotational_interpolants[i][k];

      F rzp[3] = {0};
      m2->get_coords(k, rzp);
      m2->magnetic_map(rzp, dphi * i, nsteps_per_dphi * i);

      F mu[3] = {0};
      I tid = m2->locate(rzp, mu);

      rotational_interpolants[i][k] = std::make_tuple(tid, mu[0], mu[1]);
    });
  }
}

template <typename I, typename F>
ndarray<F> simplicial_xgc_3d_mesh<I, F>::rotate(
    const ndarray<F>& scalar) const
{
  const size_t m2n0 = m2->n(0);

  ndarray<F> f;
  f.reshape(scalar);

  for (int i = 0; i < nphi; i ++) {
    if (i == 0) {
      for (int j = 0; j < m2n0; j ++)
        f[j] = scalar[j];
    } else {
      for (int j = 0; j < m2n0; j ++) {
        const auto &l = rotational_interpolants[i][j];
        if (std::get<0>(l) < 0)
          f(j, i) = std::numeric_limits<F>::quiet_NaN();
        else {
          I tri[3];
          m2->get_triange(std::get<0>(l), tri);

          F mu[3] = {std::get<1>(l), std::get<2>(l), 0};
          mu[2] = F(1) - mu[0] - mu[1];

          f(j, i) =  mu[0] * scalar(tri[0], i) 
                   + mu[1] * scalar(tri[1], i)
                   + mu[2] * scalar(tri[2], i);
        }
      }
    }
  }

  return f;
}

}

#endif
