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
#include <ftk/numeric/critical_point_type.hh>
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
    const double *Sl  // scalar of last timestep
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

  void advance_timestep();
  void update_timestep();
  
#if FTK_HAVE_VTK
  virtual vtkSmartPointer<vtkPolyData> get_traced_critical_points_vtk() const;
  virtual vtkSmartPointer<vtkPolyData> get_discrete_critical_points_vtk() const;
#endif

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

  virtual void simplex_positions(const std::vector<std::vector<int>>& vertices, double X[3][3]) const;
  virtual void simplex_vectors(const std::vector<std::vector<int>>& vertices, double v[3][2]) const;
  virtual void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[3]) const;
  virtual void simplex_jacobians(const std::vector<std::vector<int>>& vertices, 
      double Js[3][2][2]) const;
};


////////////////////
void critical_point_tracker_2d_regular::initialize()
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

void critical_point_tracker_2d_regular::finalize()
{
  diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, 0);

  if (comm.rank() == 0) {
    fprintf(stderr, "finalizing...\n");
    // trace_intersections();
    trace_connected_components();
  }
}

void critical_point_tracker_2d_regular::advance_timestep()
{
  update_timestep();

  const int nt = 2;
  if (scalar.size() > nt) scalar.pop_back();
  if (V.size() > nt) V.pop_back();
  if (gradV.size() > nt) gradV.pop_back();

  current_timestep ++;
}

void critical_point_tracker_2d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  // derive fields
  if (scalar_field_source == SOURCE_GIVEN) {
    if (vector_field_source == SOURCE_DERIVED) push_input_vector_field(gradient2D(scalar[0])); // 0 is the current timestep; 1 is the last timestep
    if (jacobian_field_source == SOURCE_DERIVED) push_input_jacobian_field(jacobian2D(V[0]));
  }

  // scan 2-simplices
  // fprintf(stderr, "tracking 2D critical points...\n");
  auto func2 = [=](element_t e) {
      critical_point_2dt_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        discrete_critical_points[e] = cp;
      }
    };

  if (xl == FTK_XL_NONE) {
    if (V.size() >= 2) { // interval
      // m.element_for_interval(2, current_timestep-1, current_timestep, func2);
      m.element_for(2, lattice({
            local_domain.start(0), 
            local_domain.start(1), 
            static_cast<size_t>(current_timestep - 1), 
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
          scalar[1].data()
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
        scalar[0].data()
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

void critical_point_tracker_2d_regular::trace_intersections()
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

void critical_point_tracker_2d_regular::trace_connected_components()
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

#if 1
  std::set<element_t> elements;
  for (const auto &kv : discrete_critical_points)
    elements.insert(kv.first);
  connected_components = extract_connected_components<element_t, std::set<element_t>>(
      neighbors, elements);
#endif

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

void critical_point_tracker_2d_regular::simplex_positions(
    const std::vector<std::vector<int>>& vertices, double X[3][3]) const
{
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j];
}

void critical_point_tracker_2d_regular::simplex_vectors(
    const std::vector<std::vector<int>>& vertices, double v[3][2]) const
{
  for (int i = 0; i < 3; i ++) {
    const int iv = vertices[i][2] == current_timestep ? 0 : 1;
    for (int j = 0; j < 2; j ++)
      v[i][j] = V[iv](j, 
          vertices[i][0] - local_array_domain.start(0), 
          vertices[i][1] - local_array_domain.start(1));
  }
}

void critical_point_tracker_2d_regular::simplex_scalars(
    const std::vector<std::vector<int>>& vertices, double values[3]) const
{
  for (int i = 0; i < 3; i ++) {
    const int iv = vertices[i][2] == current_timestep ? 0 : 1;
    values[i] = scalar[iv](
        vertices[i][0] - local_array_domain.start(0), 
        vertices[i][1] - local_array_domain.start(1));
  }
}

void critical_point_tracker_2d_regular::simplex_jacobians(
    const std::vector<std::vector<int>>& vertices, 
    double Js[3][2][2]) const
{
  for (int i = 0; i < 3; i ++) {
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

bool critical_point_tracker_2d_regular::check_simplex(
    const regular_simplex_mesh_element& e,
    critical_point_2dt_t& cp)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex
 
  double v[3][2]; // obtain vector values
  simplex_vectors(vertices, v);
 
  double mu[3]; // check intersection
  bool succ = inverse_lerp_s2v2(v, mu);
  if (!succ) return false;

  double X[3][3]; // position
  simplex_positions(vertices, X);
  lerp_s2v3(X, mu, cp.x);

  return true; // TODO: jacobians

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
  } else {
    // TODO: jacobian is not given
  }
  cp.type = critical_point_type_2d(J, is_jacobian_field_symmetric);
  if (cp.type & type_filter) return true;
  else return false;
} 

#if FTK_HAVE_VTK
vtkSmartPointer<vtkPolyData> critical_point_tracker_2d_regular::get_traced_critical_points_vtk() const
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

vtkSmartPointer<vtkPolyData> critical_point_tracker_2d_regular::get_discrete_critical_points_vtk() const
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
}


namespace diy {
  template <> struct Serialization<ftk::critical_point_2dt_t> {
    static void save(diy::BinaryBuffer& bb, const ftk::critical_point_2dt_t &cp) {
      diy::save(bb, cp.x[0]);
      diy::save(bb, cp.x[1]);
      diy::save(bb, cp.x[2]);
      diy::save(bb, cp.scalar);
      diy::save(bb, cp.type);
    }

    static void load(diy::BinaryBuffer& bb, ftk::critical_point_2dt_t& cp) {
      diy::load(bb, cp.x[0]);
      diy::load(bb, cp.x[1]);
      diy::load(bb, cp.x[2]);
      diy::load(bb, cp.scalar);
      diy::load(bb, cp.type);
    }
  };
}

#endif
