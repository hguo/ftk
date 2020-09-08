#ifndef _FTK_PARALLEL_VECTOR_TRACKER_HH
#define _FTK_PARALLEL_VECTOR_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/parallel_vector_curve_set.hh>
#include <ftk/numeric/parallel_vector_solver3.hh>
#include <ftk/ndarray/grad.hh>
#include <deque>

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

namespace ftk {

struct parallel_vector_tracker_3d_regular : public filter 
{
  parallel_vector_tracker_3d_regular() : m(4) {};
  
  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void push_field_data_snapshot(
      const ndarray<double> &v, 
      const ndarray<double> &w, 
      const ndarray<double> &Jv);

  void update_timestep();
  void advance_timestep();

  void initialize();
  void finalize();
  void update() {}

public:
  void write_discrete_pvs_vtk(const std::string& filename) const;
  void write_traced_pvs_vtk(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_discrete_pvs_vtk() const;
#endif

protected:
  simplicial_regular_mesh m;
  typedef simplicial_regular_mesh_element element_t;
  typedef parallel_vector_point_t pv_t;
  // std::multimap<element_t, parallel_vector_point_t> discrete_pvs; // discrete parallel vector points
  std::map<element_t, parallel_vector_point_t> discrete_pvs; // discrete parallel vector points

  parallel_vector_curve_set_t traced_pvs;
  
  struct field_data_snapshot_t {
    ndarray<double> v,  w;
    ndarray<double> Jv; // jacobian of v
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;

  template <int n, typename T>
  void simplex_values(const element_t &e, 
      T X[n][4], T V[n][3], T W[n][3], T J[n][3][3]) const;

private:
  void check_simplex(const element_t& s); // , pv_t& cp);

  void trace_curves();
  void trace_surfaces();

protected:
  lattice domain, array_domain;
  int current_timestep = 0;
  int start_timestep = 0, end_timestep = std::numeric_limits<int>::max();
};

/////
inline void parallel_vector_tracker_3d_regular::initialize()
{
  // initializing bounds
  m.set_lb_ub({
      static_cast<int>(domain.start(0)),
      static_cast<int>(domain.start(1)),
      static_cast<int>(domain.start(2)),
      start_timestep
    }, {
      static_cast<int>(domain.size(0)),
      static_cast<int>(domain.size(1)),
      static_cast<int>(domain.size(2)),
      end_timestep
    });
}

inline void parallel_vector_tracker_3d_regular::finalize()
{
  fprintf(stderr, "#dpvs=%zu\n", discrete_pvs.size());
  // trace_surfaces();
  trace_curves();
}

inline void parallel_vector_tracker_3d_regular::update_timestep()
{
  fprintf(stderr, "current_timestep = %d, #snapshots=%zu\n", current_timestep, field_data_snapshots.size());
  
  discrete_pvs.clear();
 
  using namespace std::placeholders;
  // m.element_for_ordinal(2, current_timestep, 
  //     std::bind(&parallel_vector_tracker_3d_regular::check_simplex, this, _1));
  
  m.element_for(2, lattice({ // ordinal
        domain.start(0), 
        domain.start(1), 
        domain.start(2), 
        static_cast<size_t>(current_timestep), 
      }, {
        domain.size(0), 
        domain.size(1), 
        domain.size(2), 
        1
      }), 
      ftk::ELEMENT_SCOPE_ORDINAL, 
      std::bind(&parallel_vector_tracker_3d_regular::check_simplex, this, _1), 
      nthreads);

  fprintf(stderr, "#dpvs=%zu\n", discrete_pvs.size());

#if 0
  if (field_data_snapshots.size() >= 2)
    m.element_for_interval(2, current_timestep, current_timestep+1,
        std::bind(&parallel_vector_tracker_3d_regular::check_simplex, this, _1));
#endif
}

inline void parallel_vector_tracker_3d_regular::advance_timestep()
{
  update_timestep();
  if (field_data_snapshots.size() >= 2)
    field_data_snapshots.pop_front();

  current_timestep ++;
}

inline void parallel_vector_tracker_3d_regular::push_field_data_snapshot(
    const ndarray<double>& v,
    const ndarray<double>& w, 
    const ndarray<double>& Jv)
{
  field_data_snapshot_t snapshot;
  snapshot.v = v;
  snapshot.w = w;
  snapshot.Jv = Jv;

  field_data_snapshots.emplace_back(snapshot);
}

template <int n/*number of vertices*/, typename T>
void parallel_vector_tracker_3d_regular::simplex_values(
    const element_t &e, 
    T X[n][4], T V[n][3], T W[n][3], T Jv[n][3][3]) const
{
  const auto &vertices = e.vertices(m);
  assert(n == vertices.size());
  
  for (int i = 0; i < n; i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    const auto &data = field_data_snapshots[iv];
   
    // v and w
    for (int j = 0; j < 3; j ++) {
      V[i][j] = data.v(j, vertices[i][0], vertices[i][1], vertices[i][2]);
      W[i][j] = data.w(j, vertices[i][0], vertices[i][1], vertices[i][2]);
    }

    // Jv
    if (!data.Jv.empty()) {
      for (int j = 0; j < 3; j ++)
        for (int k = 0; k < 3; k ++)
          Jv[i][j][k] = data.Jv(j, k, vertices[i][0], 
              vertices[i][1], vertices[i][2]);
    }
    
    // coordinates
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
  }
}

inline void parallel_vector_tracker_3d_regular::check_simplex(
    const simplicial_regular_mesh_element& e) // 2-simplex
{
  if (!e.valid(m)) return; // check if the 2-simplex is valid
  // std::cerr << e << std::endl;
  const auto &vertices = e.vertices(m);

  double X[3][4], V[3][3], W[3][3], JV[3][3][3];
  simplex_values<3>(e, X, V, W, JV);

  double lambdas[3], mus[3][3], cond;
  int ns = solve_pv_s2v3<double>(V, W, lambdas, mus, cond);

  // if (ns > 1) return;
  for (int i = 0; i < ns; i ++) {
    parallel_vector_point_t pv;

#if 0 // check eigsystem
    double Jv[3][3];
    lerp_s2v3x3(JV, mus[i], Jv);

    double P[4], roots[3];
    characteristic_polynomial_3x3(Jv, P);
    int nr = solve_cubic_real(P[2], P[1], P[0], roots);
    if (nr > 1) continue;
#endif

    double x[4], v[3], w[3];
    lerp_s2v4(X, mus[i], x);
    lerp_s2v3(V, mus[i], v);
    lerp_s2v3(W, mus[i], w);

    for (int i = 0; i < 3; i ++) {
      pv.x[i] = x[i];
      pv.v[i] = v[i];
      pv.w[i] = w[i];
    }
    
    pv.t = x[3];
    pv.lambda = lambdas[i];
    pv.cond = cond;
    pv.ordinal = e.is_ordinal(m);
    pv.tag = e.to_integer(m);
    if (!pv.ordinal) fprintf(stderr, "fuck\n");
     
    if (1) // pv.lambda > 0)
    {
      std::lock_guard<std::mutex> guard(mutex);
      discrete_pvs.insert(std::pair<element_t, pv_t>(e, pv));
    }
  }


#if 0
  double lambda, mu[3], cond;
  bool succ = solve_sujudi_haimes(V, W, lambda, mu, cond);
  // bool succ = solve_ridge(V, W, lambda, mu, cond);

  if (succ) {
    parallel_vector_point_t pv;
    double x[4], v[3], w[3];
    lerp_s2v4(X, mu, x);
    lerp_s2v3(V, mu, v);
    lerp_s2v3(W, mu, w);

    for (int i = 0; i < 3; i ++) {
      pv.x[i] = x[i];
      pv.v[i] = v[i];
      pv.w[i] = w[i];
    }
    
    pv.t = x[3];
    pv.lambda = lambda;
    pv.cond = cond;
    pv.ordinal = e.is_ordinal(m);
     
    if (1) // pv.lambda > 0)
    {
      std::lock_guard<std::mutex> guard(mutex);
      discrete_pvs.insert(std::pair<element_t, pv_t>(e, pv));
    }
  }
#endif
}

inline void parallel_vector_tracker_3d_regular::trace_surfaces()
{
  // Convert connected components to geometries
  auto neighbors = [&](element_t f) {
    std::set<element_t> neighbors;
    const auto cells = f.side_of(m);
    for (const auto c : cells) {
      const auto elements = c.sides(m);
      for (const auto f1 : elements)
        neighbors.insert(f1);
    }
    return neighbors;
  };

  std::set<element_t> elements;
  for (const auto &kv : discrete_pvs)
    elements.insert(kv.first);
  auto connected_components = extract_connected_components<element_t, std::set<element_t>>(
      neighbors, elements);

  fprintf(stderr, "#cc=%zu\n", connected_components.size());

  int k = 0;
  for (const auto &cc : connected_components) {
    if (cc.size() > 200) {
      for (const auto &e : cc) {
        if (discrete_pvs[e].ordinal)
          discrete_pvs[e].id = k;
        else 
          discrete_pvs.erase(e);
      }
      k ++;
    } else {
      for (const auto &e : cc)
        discrete_pvs.erase(e);
    }
  }
  fprintf(stderr, "#cc2=%d\n", k);

  trace_curves();
}

inline void parallel_vector_tracker_3d_regular::trace_curves()
{
  // Convert connected components to geometries
  auto neighbors = [&](element_t f) {
    std::set<element_t> neighbors;
    const auto cells = f.side_of(m);
    for (const auto c : cells) {
      const auto elements = c.sides(m);
      for (const auto f1 : elements)
        neighbors.insert(f1);
    }
    return neighbors;
  };

  std::set<element_t> elements;
  for (const auto &kv : discrete_pvs)
    elements.insert(kv.first);
  auto connected_components = extract_connected_components<element_t, std::set<element_t>>(
      neighbors, elements);

  fprintf(stderr, "#cc=%zu\n", connected_components.size());

  for (const auto &component : connected_components) {
    std::vector<std::vector<double>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(component, neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      parallel_vector_curve_t curve; 
      for (int k = 0; k < linear_graphs[j].size(); k ++)
        curve.push_back(discrete_pvs[linear_graphs[j][k]]);
      traced_pvs.add(curve);
    }
  }
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> parallel_vector_tracker_3d_regular::get_discrete_pvs_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  
  vtkIdType pid[1];
  
  // const auto critical_points = get_critical_points();
  for (const auto &kv : discrete_pvs) {
    const auto &pv = kv.second;
    double p[3] = {pv.x[0], pv.x[1], pv.x[2]}; // TODO: time
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
 
  const size_t nv = discrete_pvs.size();

  if (1) { // ids
    vtkSmartPointer<vtkUnsignedIntArray> ids = vtkSmartPointer<vtkUnsignedIntArray>::New();
    ids->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &kv : discrete_pvs)
      ids->SetValue(i ++, kv.second.id);
    ids->SetName("id");
    polyData->GetPointData()->AddArray(ids);
  }
  
  if (1) { // ordinal
    vtkSmartPointer<vtkUnsignedIntArray> ordinal = vtkSmartPointer<vtkUnsignedIntArray>::New();
    ordinal->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &kv : discrete_pvs)
      ordinal->SetValue(i ++, kv.second.ordinal);
    ordinal->SetName("ordinal");
    polyData->GetPointData()->AddArray(ordinal);
  }
  
  if (1) { // lambda
    vtkSmartPointer<vtkDoubleArray> lambdas = vtkSmartPointer<vtkDoubleArray>::New();
    lambdas->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &kv : discrete_pvs)
      lambdas->SetValue(i ++, kv.second.lambda);
    lambdas->SetName("lambda");
    polyData->GetPointData()->AddArray(lambdas);
  }
  
  if (1) { // v
    vtkSmartPointer<vtkDoubleArray> v = vtkSmartPointer<vtkDoubleArray>::New();
    v->SetNumberOfComponents(3);
    v->SetNumberOfTuples(nv);
    size_t i = 0;
    for (const auto &kv : discrete_pvs)
      v->SetTuple3(i ++, kv.second.v[0], kv.second.v[1], kv.second.v[2]);
    v->SetName("v");
    polyData->GetPointData()->AddArray(v);
  }
  
  if (1) { // w
    vtkSmartPointer<vtkDoubleArray> w = vtkSmartPointer<vtkDoubleArray>::New();
    w->SetNumberOfComponents(3);
    w->SetNumberOfTuples(nv);
    size_t i = 0;
    for (const auto &kv : discrete_pvs)
      w->SetTuple3(i ++, kv.second.w[0], kv.second.w[1], kv.second.w[2]);
    w->SetName("w");
    polyData->GetPointData()->AddArray(w);
  }

  if (1) { // condition numbers
    vtkSmartPointer<vtkDoubleArray> conds = vtkSmartPointer<vtkDoubleArray>::New();
    conds->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &kv : discrete_pvs)
      conds->SetValue(i ++, kv.second.cond);
    conds->SetName("cond");
    polyData->GetPointData()->AddArray(conds);
  }
  
  if (1) { // time
    vtkSmartPointer<vtkDoubleArray> time = vtkSmartPointer<vtkDoubleArray>::New();
    time->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &kv : discrete_pvs)
      time->SetValue(i ++, kv.second.t);
    time->SetName("time");
    polyData->GetPointData()->AddArray(time);
  }

  return polyData;
}

inline void parallel_vector_tracker_3d_regular::write_discrete_pvs_vtk(const std::string& filename) const
{
  if (is_root_proc()) {
    auto poly = get_discrete_pvs_vtk();
    write_vtp(filename, poly);
  }
}

inline void parallel_vector_tracker_3d_regular::write_traced_pvs_vtk(const std::string& filename) const
{
  if (is_root_proc()) {
    auto poly = traced_pvs.to_vtp();
    write_vtp(filename, poly);
  }
}
#else
inline void parallel_vector_tracker_3d_regular::write_discrete_pvs_vtk(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}

inline void parallel_vector_tracker_3d_regular::write_traced_pvs_vtk(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

}

#endif
