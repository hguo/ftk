#ifndef _FTK_CRITICAL_POINT_TRACKER_2D_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_2D_REGULAR_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/filters/critical_point_tracker_regular.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/hypermesh/regular_simplex_mesh.hh>
#include <ftk/external/diy/serialization.hpp>

#if FTK_HAVE_VTK
#include <vtkUnsignedIntArray.h>
#endif

#if FTK_HAVE_CUDA
extern std::vector<ftk::critical_point_t<3, double>> 
extract_cp2dt_cuda(
    int scope, int current_timestep, 
    const ftk::lattice& domain, // 3D
    const ftk::lattice& core, // 3D
    const ftk::lattice& ext, // 2D, array dimension
    const double *Vc, // current timestep
    const double *Vl, // last timestep
    const double *Jc, // jacobian of current timestep
    const double *Jl, // jacobian of last timestep
    const double *Sc, // scalar of current timestep
    const double *Sl, // scalar of last timestep
    bool use_explicit_coords,
    const double *coords // coords of vertices
  ); 
#endif

namespace ftk {

typedef critical_point_t<3, double> critical_point_2dt_t;

struct critical_point_tracker_2d_regular : public critical_point_tracker_regular {
  critical_point_tracker_2d_regular() : m(3) {}
  critical_point_tracker_2d_regular(int argc, char **argv) 
    : critical_point_tracker_regular(argc, argv), m(3) {}
  virtual ~critical_point_tracker_2d_regular() {}

  void initialize();
  void finalize();

  void update_timestep();
 
  void push_snapshot_scalar_field(const ndarray<double>&);
  void push_snapshot_vector_field(const ndarray<double>&);
  // void push_snapshot_jacobian_field(const ndarray<double>&); // TODO

#if FTK_HAVE_VTK
  virtual vtkSmartPointer<vtkPolyData> get_traced_critical_points_vtk() const;
  virtual vtkSmartPointer<vtkPolyData> get_discrete_critical_points_vtk() const;
#endif

  void write_discrete_critical_points(const std::string& filename) const;
  void write_traced_critical_points(const std::string& filename) const;

protected:
  regular_simplex_mesh m;
  
  unsigned int type_filter = 0xffffffff;

  typedef regular_simplex_mesh_element element_t;
  
  std::map<element_t, critical_point_2dt_t> discrete_critical_points;
  std::vector<std::set<element_t>> connected_components;
  std::vector<std::vector<critical_point_2dt_t>> traced_critical_points;

protected:
  bool check_simplex(const element_t& s, critical_point_2dt_t& cp);
  void trace_intersections();
  void trace_connected_components();

  template <typename I=int> void simplex_indices(const std::vector<std::vector<int>>& vertices, I indices[]) const;
  virtual void simplex_coordinates(const std::vector<std::vector<int>>& vertices, double X[][3]) const;
  template <typename T=double> void simplex_vectors(const std::vector<std::vector<int>>& vertices, T v[][2]) const;
  virtual void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[]) const;
  virtual void simplex_jacobians(const std::vector<std::vector<int>>& vertices, 
      double Js[][2][2]) const;

protected: // working in progress
  bool robust_check_simplex0(const element_t& s, critical_point_2dt_t &cp);
  bool robust_check_simplex1(const element_t& s, critical_point_2dt_t &cp);
  bool robust_check_simplex2(const element_t& s, critical_point_2dt_t &cp);
};


////////////////////
inline void critical_point_tracker_2d_regular::initialize()
{
  // initializing bounds
  m.set_lb_ub({
      static_cast<int>(domain.start(0)),
      static_cast<int>(domain.start(1)),
      start_timestep
    }, {
      static_cast<int>(domain.size(0)),
      static_cast<int>(domain.size(1)),
      end_timestep
    });

  if (use_default_domain_partition) {
    lattice_partitioner partitioner(domain);
    
    // a ghost size of 2 is necessary for jacobian derivaition; 
    // even if jacobian is not necessary, a ghost size of 1 is 
    // necessary for accessing values on boundaries
    partitioner.partition(comm.size(), {}, {2, 2});

    local_domain = partitioner.get_core(comm.rank());
    local_array_domain = partitioner.get_ext(comm.rank());
  }

  if (!is_input_array_partial)
    local_array_domain = array_domain;
}

inline void critical_point_tracker_2d_regular::finalize()
{
  diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, 0);

  if (comm.rank() == 0) {
    fprintf(stderr, "finalizing...\n");
    // trace_intersections();
    trace_connected_components();
  }
}

inline void critical_point_tracker_2d_regular::push_snapshot_scalar_field(const ndarray<double>& s)
{
  scalar.push_back(s);
  if (scalar_field_source == SOURCE_GIVEN) {
    if (vector_field_source == SOURCE_DERIVED) push_snapshot_vector_field(gradient2D(scalar[scalar.size()-1])); // 0 is the current timestep; 1 is the last timestep
    if (jacobian_field_source == SOURCE_DERIVED) push_snapshot_jacobian_field(jacobian2D(V[V.size()-1]));
  }
}

inline void critical_point_tracker_2d_regular::push_snapshot_vector_field(const ndarray<double>& v)
{
  V.push_back(v);
  if (vector_field_source == SOURCE_GIVEN) {
    if (jacobian_field_source == SOURCE_DERIVED) push_snapshot_jacobian_field(jacobian2D(V[0]));
  }
}

inline void critical_point_tracker_2d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  auto func0 = [=](element_t e) {
      critical_point_2dt_t cp;
      if (robust_check_simplex0(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        if (filter_critical_point_type(cp))
          discrete_critical_points[e] = cp;
      }
    };
  
  auto func1 = [=](element_t e) {
      critical_point_2dt_t cp;
      if (robust_check_simplex1(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        if (filter_critical_point_type(cp))
          discrete_critical_points[e] = cp;
      }
    };

  // scan 2-simplices
  // fprintf(stderr, "tracking 2D critical points...\n");
  auto func2 = [=](element_t e) {
      critical_point_2dt_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        if (filter_critical_point_type(cp))
          discrete_critical_points[e] = cp;
      }
    };

  if (xl == FTK_XL_NONE) {
    if (V.size() >= 2) { // interval
      // m.element_for_interval(2, current_timestep-1, current_timestep, func2);
      m.element_for(2, lattice({
            local_domain.start(0), 
            local_domain.start(1), 
            // static_cast<size_t>(current_timestep - 1), 
            static_cast<size_t>(current_timestep),
          }, {
            local_domain.size(0), 
            local_domain.size(1), 
            1
          }),
          ftk::ELEMENT_SCOPE_INTERVAL, 
          func2, nthreads);
    }

    // m.element_for_ordinal(2, current_timestep, func2);
    m.element_for(2, lattice({ // ordinal
          local_domain.start(0), 
          local_domain.start(1), 
          static_cast<size_t>(current_timestep), 
        }, {
          local_domain.size(0), 
          local_domain.size(1), 
          1
        }), 
        ftk::ELEMENT_SCOPE_ORDINAL, 
        func2, nthreads);
  } else if (xl == FTK_XL_CUDA) {
#if FTK_HAVE_CUDA
    ftk::lattice domain3({
          domain.start(0), 
          domain.start(1), 
          0
        }, {
          domain.size(0)-1,
          domain.size(1)-1,
          std::numeric_limits<int>::max()
        });

    ftk::lattice ordinal_core({
          local_domain.start(0), 
          local_domain.start(1), 
          static_cast<size_t>(current_timestep), 
        }, {
          local_domain.size(0), 
          local_domain.size(1), 
          1
        });

    ftk::lattice interval_core({
          local_domain.start(0), 
          local_domain.start(1), 
          static_cast<size_t>(current_timestep-1), 
        }, {
          local_domain.size(0), 
          local_domain.size(1), 
          1
        });

    ftk::lattice ext({0, 0}, 
        {V[0].dim(1), V[0].dim(2)});

    if (V.size() >= 2) { // interval
      fprintf(stderr, "processing interval %d, %d\n", current_timestep - 1, current_timestep);
      auto results = extract_cp2dt_cuda(
          ELEMENT_SCOPE_INTERVAL, 
          current_timestep,
          domain3, 
          interval_core,
          ext,
          V[0].data(), // current
          V[1].data(), // last
          gradV[0].data(), 
          gradV[1].data(),
          scalar[0].data(),
          scalar[1].data(),
          use_explicit_coords, 
          coords.data()
        );
      fprintf(stderr, "interal_results#=%d\n", results.size());
      for (auto cp : results) {
        element_t e(3, 2);
        e.from_work_index(m, cp.tag, interval_core, ELEMENT_SCOPE_INTERVAL);
        discrete_critical_points[e] = cp;
      }
    }

    // ordinal
    auto results = extract_cp2dt_cuda(
        ELEMENT_SCOPE_ORDINAL, 
        current_timestep, 
        domain3,
        ordinal_core,
        ext,
        V[0].data(),
        V[0].data(),
        gradV[0].data(),
        gradV[0].data(),
        scalar[0].data(),
        scalar[0].data(),
        use_explicit_coords, 
        coords.data()
      );
    
    for (auto cp : results) {
      element_t e(3, 2);
      e.from_work_index(m, cp.tag, ordinal_core, ELEMENT_SCOPE_ORDINAL);
      discrete_critical_points[e] = cp;
    }
#else
    assert(false);
#endif
  }
}

inline void critical_point_tracker_2d_regular::trace_intersections()
{
  // scan 3-simplices to get connected components
  union_find<element_t> uf;
  for (const auto &kv : discrete_critical_points) 
    uf.add(kv.first);

  m.element_for(3, [&](const regular_simplex_mesh_element& f) {
    const auto sides = f.sides(m);
    std::set<element_t> intersected_sides;

    for (const auto& side : sides)
      if (uf.has(side)) 
        intersected_sides.insert(side);

    if (intersected_sides.size() > 1) {// the size of intersected_size should be 0 or 2
      for (auto it = std::next(intersected_sides.begin(), 1); it != intersected_sides.end(); it ++) {
        std::lock_guard<std::mutex> guard(mutex);
        uf.unite(*intersected_sides.begin(), *it);
      }
    }
  });
  uf.get_sets(connected_components);
}

inline void critical_point_tracker_2d_regular::trace_connected_components()
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
  for (const auto &kv : discrete_critical_points)
    elements.insert(kv.first);
  connected_components = extract_connected_components<element_t, std::set<element_t>>(
      neighbors, elements);

  for (const auto &component : connected_components) {
    std::vector<std::vector<double>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(component, neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      std::vector<critical_point_2dt_t> traj; 
      for (int k = 0; k < linear_graphs[j].size(); k ++)
        traj.push_back(discrete_critical_points[linear_graphs[j][k]]);
      traced_critical_points.emplace_back(traj);
    }
  }
}

template <typename I>
inline void critical_point_tracker_2d_regular::simplex_indices(
    const std::vector<std::vector<int>>& vertices, I indices[]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    indices[i] = m.get_lattice().to_integer(vertices[i]);
}

inline void critical_point_tracker_2d_regular::simplex_coordinates(
    const std::vector<std::vector<int>>& vertices, double X[][3]) const
{
  if (use_explicit_coords) {
    for (int i = 0; i < vertices.size(); i ++) {
      for (int j = 0; j < 2; j ++) 
        X[i][j] = coords(j, vertices[i][0], vertices[i][1]);
      X[i][2] = vertices[i][2];
    }
  } else {
    for (int i = 0; i < vertices.size(); i ++)
      for (int j = 0; j < 3; j ++)
        X[i][j] = vertices[i][j];
  }
}

template <typename T>
inline void critical_point_tracker_2d_regular::simplex_vectors(
    const std::vector<std::vector<int>>& vertices, T v[][2]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][2] == current_timestep ? 0 : 1;
    for (int j = 0; j < 2; j ++)
      v[i][j] = V[iv](j, 
          vertices[i][0] - local_array_domain.start(0), 
          vertices[i][1] - local_array_domain.start(1));
  }
}

inline void critical_point_tracker_2d_regular::simplex_scalars(
    const std::vector<std::vector<int>>& vertices, double values[]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][2] == current_timestep ? 0 : 1;
    values[i] = scalar[iv](
        vertices[i][0] - local_array_domain.start(0), 
        vertices[i][1] - local_array_domain.start(1));
  }
}

inline void critical_point_tracker_2d_regular::simplex_jacobians(
    const std::vector<std::vector<int>>& vertices, 
    double Js[][2][2]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][2] == 0 ? 0 : 1;
    for (int j = 0; j < 2; j ++) {
      for (int k = 0; k < 2; k ++) {
        Js[i][j][k] = gradV[iv](k, j, 
            vertices[i][0] - local_array_domain.start(0), 
            vertices[i][1] - local_array_domain.start(1));
      }
    }
  }
}

inline bool critical_point_tracker_2d_regular::check_simplex(
    const regular_simplex_mesh_element& e,
    critical_point_2dt_t& cp)
{
  typedef fixed_point<> fp_t;
  
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex
 
  double v[3][2]; // obtain vector values
  simplex_vectors(vertices, v);

#if 0 // working in progress, handling degeneracy cases
  fp_t vf[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++)
      vf[i][j] = v[i][j];
  
  const fp_t M[3][3] = {
    {vf[0][0], vf[1][0], vf[2][0]},
    {vf[0][1], vf[1][1], vf[2][1]},
    {1, 1, 1}
  };
  fp_t adjM[3][3];
  fp_t detM = ftk::det3(M);
  // std::cerr << "detM=" << detM << std::endl;
  adjugate3(M, adjM);
  // ftk::print3x3("M", M);
  // ftk::print3x3("adjM", adjM);
  // long long b[3] = {adjM[0][2]*factor, adjM[1][2]*factor, adjM[2][2]*factor};
  fp_t b[3] = {adjM[0][2], adjM[1][2], adjM[2][2]};

  int sign_detM = ftk::sign(detM);
  if (sign_detM < 0) {
    detM *= sign_detM;
    for (int k = 0; k < 3; k ++)
      b[k] *= sign_detM;
  }
  bool succ1 = (b[0] > 0 && b[0] < detM && b[1] > 0 && b[1] < detM && b[2] > 0 && b[2] < detM);
  if (!succ1) return false;
#endif

  // robust critical point test
  fp_t vf[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++)
      vf[i][j] = v[i][j];
  int indices[3];
  simplex_indices(vertices, indices);
  bool succ = robust_critical_point_in_simplex2(vf, indices);
  if (!succ) return false;

  double mu[3]; // check intersection
  bool succ2 = inverse_lerp_s2v2(v, mu);
  // if (!succ2) return false;

  double X[3][3]; // position
  simplex_coordinates(vertices, X);
  lerp_s2v3(X, mu, cp.x);

  if (scalar_field_source != SOURCE_NONE) {
    double values[3];
    simplex_scalars(vertices, values);
    cp.scalar = lerp_s2(values, mu);
  }

  double J[2][2] = {0}; // jacobian
  if (jacobian_field_source != SOURCE_NONE) { // lerp jacobian
    double Js[3][2][2];
    simplex_jacobians(vertices, Js);
    lerp_s2m2x2(Js, mu, J);
    ftk::make_symmetric2x2(J); // TODO
  } else {
    // TODO: jacobian is not given
  }
#if 0
  double X2[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++)
      X2[i][j] = X[i][j];
  jacobian_3dsimplex2(X2, v, J);
  ftk::make_symmetric2x2(J); // TODO
#endif
  cp.type = critical_point_type_2d(J, is_jacobian_field_symmetric);

  return true;
} 

inline bool critical_point_tracker_2d_regular::robust_check_simplex0(const element_t& e, critical_point_2dt_t& cp)
{
  typedef fixed_point<> fp_t;

  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex

  fp_t v[1][2]; // obtain vector values
  simplex_vectors<fp_t>(vertices, v);

  if (v[0][0] == 0 && v[0][1] == 0) {
    fprintf(stderr, "zero!\n");
  }
  return false;
#if 0 // TODO
  bool succ = inverse_lerp_s2v2(v, mu);
  if (!succ) return false;

  double X[3][3]; // position
  simplex_coordinates(vertices, X);
  lerp_s2v3(X, mu, cp.x);
#endif
}

inline bool critical_point_tracker_2d_regular::robust_check_simplex1(const element_t& e, critical_point_2dt_t& cp)
{
  typedef fixed_point<> fp_t;

  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex

  // fp_t v[2][2]; // obtain vector values
  // simplex_vectors<fp_t>(vertices, v);
 
  double V[2][2];
  long long iV[2][2];
  simplex_vectors(vertices, V);
  for (int i = 0; i < 2; i ++)
    for (int j = 0; j < 2; j ++)
      iV[i][j] = V[i][j] * 32768;

  bool succ = ftk::integer_inverse_lerp_s1v2(iV);
  if (succ) fprintf(stderr, "shoot\n");

  return false;
#if 0 // TODO
  bool succ = inverse_lerp_s2v2(v, mu);
  if (!succ) return false;

  double X[3][3]; // position
  simplex_coordinates(vertices, X);
  lerp_s2v3(X, mu, cp.x);
#endif
}

#if 0
void critical_point_tracker_2d_regular::robust_check_simplex2(const element_t& s, critical_point_2dt_t& cp)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex
 
  double v[3][2]; // obtain vector values
  simplex_vectors(vertices, v);
 
  double mu[3]; // check intersection
  bool succ = inverse_lerp_s2v2(v, mu);
  if (!succ) return false;

  double X[3][3]; // position
  simplex_coordinates(vertices, X);
  lerp_s2v3(X, mu, cp.x);
}
#endif

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> critical_point_tracker_2d_regular::get_traced_critical_points_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();

  for (const auto &curve : traced_critical_points)
    for (auto i = 0; i < curve.size(); i ++) {
      double p[3] = {curve[i][0], curve[i][1], curve[i][2]};
      points->InsertNextPoint(p);
    }

  size_t nv = 0;
  for (const auto &curve : traced_critical_points) {
    vtkSmartPointer<vtkPolyLine> polyLine = vtkPolyLine::New();
    polyLine->GetPointIds()->SetNumberOfIds(curve.size());
    for (int i = 0; i < curve.size(); i ++)
      polyLine->GetPointIds()->SetId(i, i+nv);

    cells->InsertNextCell(polyLine);
    nv += curve.size();
  }
  
  polyData->SetPoints(points);
  polyData->SetLines(cells);

  // point data for types
  if (type_filter) {
    vtkSmartPointer<vtkUnsignedIntArray> types = vtkSmartPointer<vtkUnsignedIntArray>::New();
    types->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &curve : traced_critical_points) {
      for (auto j = 0; j < curve.size(); j ++)
        types->SetValue(i ++, curve[j].type);
    }
    types->SetName("type");
    polyData->GetPointData()->AddArray(types);
  }

  if (1) { // ids
    vtkSmartPointer<vtkUnsignedIntArray> ids = vtkSmartPointer<vtkUnsignedIntArray>::New();
    ids->SetNumberOfValues(nv);
    size_t i = 0;
    for (auto k = 0; k < traced_critical_points.size(); k ++)
      for (auto j = 0; j < traced_critical_points[k].size(); j ++)
        ids->SetValue(i ++, k);
    ids->SetName("id");
    polyData->GetPointData()->AddArray(ids);

  }

  // point data for scalars
  // if (has_scalar_field) {
  if (1) { // scalar is 0 if no scalar field available
    vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
    scalars->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &curve : traced_critical_points) {
      for (auto j = 0; j < curve.size(); j ++)
        scalars->SetValue(i ++, curve[j].scalar);
    }
    scalars->SetName("scalar");
    polyData->GetPointData()->AddArray(scalars);
  }

  return polyData;
}

inline vtkSmartPointer<vtkPolyData> critical_point_tracker_2d_regular::get_discrete_critical_points_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  
  vtkIdType pid[1];
  for (const auto &kv : discrete_critical_points) {
    const auto &cp = kv.second;
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]};
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);

#if 0
  // point data for types
  vtkSmartPointer<vtkDoubleArray> types = vtkSmartPointer<vtkDoubleArray>::New();
  types->SetNumberOfValues(results.size());
  for (auto i = 0; i < results.size(); i ++) {
    types->SetValue(i, static_cast<double>(results[i].type));
  }
  types->SetName("type");
  polyData->GetPointData()->AddArray(types);
  
  // point data for scalars
  if (has_scalar_field) {
    vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
    scalars->SetNumberOfValues(results.size());
    for (auto i = 0; i < results.size(); i ++) {
      scalars->SetValue(i, static_cast<double>(results[i].scalar));
    }
    scalars->SetName("scalar");
    polyData->GetPointData()->AddArray(scalars);
  }
#endif
  return polyData;
}
#endif

inline void critical_point_tracker_2d_regular::write_discrete_critical_points(const std::string& filename) const
{
  diy::serializeToFile(discrete_critical_points, filename);
}

inline void critical_point_tracker_2d_regular::write_traced_critical_points(const std::string& filename) const 
{
  diy::serializeToFile(traced_critical_points, filename);
}

}

#endif
