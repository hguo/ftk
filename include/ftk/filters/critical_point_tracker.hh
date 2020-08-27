#ifndef _FTK_CRITICAL_POINT_TRACKER_HH
#define _FTK_CRITICAL_POINT_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/critical_point.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy-ext/gather.hh>

#if FTK_HAVE_VTK
#include <vtkUnsignedIntArray.h>
#include <vtkVertex.h>
#include <vtkSmartPointer.h>
#endif

namespace ftk {

struct critical_point_tracker : public filter {
  critical_point_tracker() {}

  virtual void update() {}; 
  void reset() {
    field_data_snapshots.clear();
    traced_critical_points.clear();
  }
  
  virtual int cpdims() const = 0;

  void set_enable_streaming_trajectories(bool b) { enable_streaming_trajectories = b; }
  void set_enable_discarding_interval_points(bool b) { enable_discarding_interval_points = b; }
  void set_enable_discarding_degenerate_points(bool b) { enable_discarding_degenerate_points = b; }
  void set_enable_ignoring_degenerate_points(bool b) { enable_ignoring_degenerate_points = b; }

  void set_type_filter(unsigned int);

  void set_scalar_components(const std::vector<std::string>& c);
  int get_num_scalar_components() const {return scalar_components.size();}

  void update_traj_statistics();
  void select_traj(std::function<bool(const critical_point_traj_t& traj)>);

  void slice_traced_critical_points(unsigned int type_filter); // slice traces after finalization

public: // outputs
  const std::vector<critical_point_traj_t>& get_traced_critical_points() const {return traced_critical_points;}
  virtual std::vector<critical_point_t> get_critical_points() const = 0;
  const std::map<int, std::vector<std::tuple<critical_point_t, int>>>& get_sliced_critical_points() const {return sliced_critical_points;}

  void write_traced_critical_points_binary(const std::string& filename) const;
  void write_traced_critical_points_text(std::ostream& os) const;
  void write_traced_critical_points_text(const std::string& filename) const;
  void write_traced_critical_points_vtk(const std::string& filename) const;
#if FTK_HAVE_VTK
  // vtkSmartPointer<vtkPolyData> get_current_active_critical_points_vtk() const;
  vtkSmartPointer<vtkPolyData> get_traced_critical_points_vtk() const;
#endif

  void write_sliced_critical_points_text(int t, std::ostream& os) const;
  void write_sliced_critical_points_text(int t, const std::string& filename) const;
  void write_sliced_critical_points_vtk(int t, const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_sliced_critical_points_vtk(int t) const;
#endif

  void write_critical_points_binary(const std::string& filename) const;
  void write_critical_points_text(std::ostream& os) const;
  void write_critical_points_text(const std::string& filename) const;
  void write_critical_points_vtk(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_critical_points_vtk() const;
#endif

public: // inputs
  bool pop_field_data_snapshot();
  virtual void push_field_data_snapshot(
      const ndarray<double> &scalar, 
      const ndarray<double> &vector,
      const ndarray<double> &jacobian);
  virtual void push_scalar_field_snapshot(const ndarray<double> &scalar); // push scalar only

  virtual void push_field_data_spacetime(
      const ndarray<double> &scalars, 
      const ndarray<double> &vectors,
      const ndarray<double> &jacobians);
  void push_scalar_field_spacetime(const ndarray<double>& scalars);
  
  virtual void set_current_timestep(int t) {current_timestep = t;}
  int get_current_timestep() const {return current_timestep;}

protected:
  template <typename I> // mesh element type
  void grow_trajectories(
      std::vector<critical_point_traj_t> &trajectories,
      std::map<I, critical_point_t> &discrete_critical_poionts, // critical point tag needs to index mesh element ID.  Discrete critical points will be cleared after tracing
      std::function<std::set<I>(I)> neighbors, 
      std::function<I(unsigned long long)> tag_to_element
  );

	template <typename I> // mesh element type
	std::vector<critical_point_traj_t> trace_critical_points_offline(
		const std::map<I, critical_point_t> &discrete_critical_points,
		std::function<std::set<I>(I)> neighbors);

protected:
  struct field_data_snapshot_t {
    ndarray<double> scalar, vector, jacobian;
  };
  
  int current_timestep = 0;

  std::deque<field_data_snapshot_t> field_data_snapshots;
  
  // std::vector<std::vector<critical_point_t>> traced_critical_points;
  std::vector<critical_point_traj_t> traced_critical_points;
  std::map<int/*time*/, std::vector<std::tuple<critical_point_t, int/*id*/>>> sliced_critical_points;

  // type filter
  bool use_type_filter = false;
  unsigned int type_filter = 0;

  // scalar components
  std::vector<std::string> scalar_components = {"scalar"};

  bool enable_streaming_trajectories = false;
  bool enable_discarding_interval_points = false;
  bool enable_discarding_degenerate_points = false;
  bool enable_ignoring_degenerate_points = false;
};

///////

inline void critical_point_tracker::push_field_data_snapshot(
    const ndarray<double>& scalar,
    const ndarray<double>& vector,
    const ndarray<double>& jacobian)
{
  field_data_snapshot_t snapshot;
  snapshot.scalar = scalar;
  snapshot.vector = vector;
  snapshot.jacobian = jacobian;

  field_data_snapshots.emplace_back(snapshot);
}

inline void critical_point_tracker::push_scalar_field_snapshot(const ndarray<double>& scalar)
{
  field_data_snapshot_t snapshot;
  snapshot.scalar = scalar;

  field_data_snapshots.emplace_back(snapshot);
}

inline void critical_point_tracker::push_field_data_spacetime(
    const ndarray<double>& scalars,
    const ndarray<double>& vectors,
    const ndarray<double>& jacobians)
{
  for (size_t t = 0; t < scalars.shape(scalars.nd()-1); t ++) {
    auto scalar = scalars.slice_time(t);
    auto vector = vectors.slice_time(t);
    auto jacobian = jacobians.slice_time(t);

    push_field_data_snapshot(scalar, vector, jacobian);
  }
}

inline void critical_point_tracker::push_scalar_field_spacetime(const ndarray<double>& scalars)
{
  for (size_t t = 0; t < scalars.shape(scalars.nd()-1); t ++)
    push_scalar_field_snapshot( scalars.slice_time(t) );
}


inline bool critical_point_tracker::pop_field_data_snapshot()
{
  if (field_data_snapshots.size() > 0) {
    field_data_snapshots.pop_front();
    return true;
  } else return false;
}

//////
#if FTK_HAVE_VTK
inline void critical_point_tracker::write_traced_critical_points_vtk(const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_traced_critical_points_vtk();
    write_vtp(filename, poly);
  }
}

inline void critical_point_tracker::write_sliced_critical_points_vtk(int k, const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_sliced_critical_points_vtk(k);
    write_vtp(filename, poly);
  }
}

inline void critical_point_tracker::write_critical_points_vtk(const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_critical_points_vtk();
    write_vtp(filename, poly);
  }
}

inline vtkSmartPointer<vtkPolyData> critical_point_tracker::get_critical_points_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  
  vtkIdType pid[1];
  
  const auto critical_points = get_critical_points();
  for (const auto &cp : get_critical_points()) {
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]}; // TODO: time
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);

#if 0 // TODO
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

inline vtkSmartPointer<vtkPolyData> critical_point_tracker::get_traced_critical_points_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> lines = vtkCellArray::New();
  vtkSmartPointer<vtkCellArray> verts = vtkCellArray::New();

  if (cpdims() == 2) {
    for (const auto &curve : traced_critical_points)
      for (auto i = 0; i < curve.size(); i ++) {
        double p[3] = {curve[i][0], curve[i][1], curve[i].t};
        points->InsertNextPoint(p);
      }
  } else { // cpdims == 3
    for (const auto &curve : traced_critical_points)
      for (auto i = 0; i < curve.size(); i ++) {
        double p[3] = {curve[i][0], curve[i][1], curve[i][2]};
        points->InsertNextPoint(p);
      }
  }

  size_t nv = 0;
  for (const auto &curve : traced_critical_points) {
    if (curve.size() < 2) { // isolated vertex
      vtkSmartPointer<vtkVertex> obj = vtkVertex::New();
      obj->GetPointIds()->SetNumberOfIds(curve.size());
      for (int i = 0; i < curve.size(); i ++)
        obj->GetPointIds()->SetId(i, i+nv);
      verts->InsertNextCell(obj);
    } else { // lines
      vtkSmartPointer<vtkPolyLine> obj = vtkPolyLine::New();
      obj->GetPointIds()->SetNumberOfIds(curve.size());
      for (int i = 0; i < curve.size(); i ++)
        obj->GetPointIds()->SetId(i, i+nv);
      lines->InsertNextCell(obj);
    }
    nv += curve.size();
  }
 
  polyData->SetPoints(points);
  polyData->SetLines(lines);
  polyData->SetVerts(verts);

  // point data for types
  if (1) { // if (type_filter) {
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
  for (auto k = 0; k < scalar_components.size(); k ++) {
    vtkSmartPointer<vtkDoubleArray> scalar = vtkSmartPointer<vtkDoubleArray>::New();
    scalar->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &curve : traced_critical_points) {
      for (auto j = 0; j < curve.size(); j ++)
        scalar->SetValue(i ++, curve[j].scalar[k]);
    }
    scalar->SetName(scalar_components[k].c_str());
    polyData->GetPointData()->AddArray(scalar);
  }

  return polyData;
}
#else
inline void critical_point_tracker::write_traced_critical_points_vtk(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}

inline void critical_point_tracker::write_critical_points_vtk(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}

inline void critical_point_tracker::write_sliced_critical_points_vtk(int, const std::string&) const 
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

inline void critical_point_tracker::write_traced_critical_points_binary(const std::string& filename) const
{
  if (is_root_proc())
  	diy::serializeToFile(traced_critical_points, filename);
}

inline void critical_point_tracker::write_traced_critical_points_text(const std::string& filename) const
{
  if (is_root_proc()) {
    std::ofstream out(filename);
    write_traced_critical_points_text(out);
    out.close();
  }
}

inline void critical_point_tracker::write_traced_critical_points_text(std::ostream& os) const
{
  os << "#trajectories=" << traced_critical_points.size() << std::endl;
  for (int i = 0; i < traced_critical_points.size(); i ++) {
    const auto &curve = traced_critical_points[i];
    os << "--trajectory " << i << ", ";
   
    if (scalar_components.size() > 0) {
      os << "min=(";
      for (int k = 0; k < scalar_components.size(); k ++)
        if (k < scalar_components.size()-1) os << curve.min[k] << ", ";
        else os << curve.min[k] << "), ";
      
      os << "max=(";
      for (int k = 0; k < scalar_components.size(); k ++)
        if (k < scalar_components.size()-1) os << curve.max[k] << ", ";
        else os << curve.max[k] << "), ";
      
      os << "persistence=(";
      for (int k = 0; k < scalar_components.size(); k ++)
        if (k < scalar_components.size()-1) os << curve.persistence[k] << ", ";
        else os << curve.persistence[k] << "), ";
    }

    os << "bbmin=(";
    for (int k = 0; k < cpdims(); k ++)
      if (k < cpdims()) os << curve.bbmin[k] << ", ";
      else os << curve.bbmin[k] << "), ";
    
    os << "bbmax=(";
    for (int k = 0; k < cpdims(); k ++)
      if (k < cpdims()) os << curve.bbmax[k] << ", ";
      else os << curve.bbmax[k] << "), ";
    
    os << "tmin=" << curve.tmin << ", tmax=" << curve.tmax << ", ";

    os << "consistent_type=" << critical_point_type_to_string(cpdims(), curve.consistent_type, scalar_components.size()) << ", ";
    os << "loop=" << curve.loop;
    os << std::endl;

    for (int k = 0; k < curve.size(); k ++) {
      os << "---";
      curve[k].print(os, cpdims(), scalar_components) << std::endl;
    }
  }
}

inline void critical_point_tracker::write_critical_points_text(const std::string& filename) const
{
  if (is_root_proc()) {
    std::ofstream out(filename);
    write_critical_points_text(out);
    out.close();
  }
}

inline void critical_point_tracker::write_sliced_critical_points_text(int k, const std::string& filename) const
{
  if (is_root_proc()) {
    std::ofstream out(filename);
    write_sliced_critical_points_text(k, out);
    out.close();
  }
}

inline void critical_point_tracker::write_critical_points_text(std::ostream& os) const
{
  for (const auto &cp : get_critical_points())
    cp.print(os, cpdims(), scalar_components) << std::endl;
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> critical_point_tracker::get_sliced_critical_points_vtk(int t) const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();

  if (sliced_critical_points.find(t) == sliced_critical_points.end()) 
    return polyData;
  const auto &lcps = sliced_critical_points.at(t);
  // fprintf(stderr, "converting %d to vtk, n=%zu\n", t, lcps.size());
  
  vtkIdType pid[1];
  for (const auto &lcp : lcps) {
    const auto &cp = std::get<0>(lcp);
    double p[3] = {cp[0], cp[1], cp[2]};
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
  
  // types & ids
  vtkSmartPointer<vtkIntArray> ids = vtkSmartPointer<vtkIntArray>::New();
  ids->SetName("id");
  ids->SetNumberOfValues(lcps.size());
  
  vtkSmartPointer<vtkUnsignedIntArray> types = vtkSmartPointer<vtkUnsignedIntArray>::New();
  types->SetName("type");
  types->SetNumberOfValues(lcps.size());

  for (auto i = 0; i < lcps.size(); i ++) {
    ids->SetValue(i, std::get<1>(lcps[i]));
    types->SetValue(i, std::get<0>(lcps[i]).type);
  }

  polyData->GetPointData()->AddArray(ids);
  polyData->GetPointData()->AddArray(types);
  polyData->GetPointData()->SetActiveScalars("id");
  
  for (auto k = 0; k < scalar_components.size(); k ++) {
    vtkSmartPointer<vtkDoubleArray> scalar = vtkSmartPointer<vtkDoubleArray>::New();
    scalar->SetNumberOfValues(lcps.size());
    for (auto i = 0; i < lcps.size(); i ++)
      scalar->SetValue(i, std::get<0>(lcps[i]).scalar[k]);
    scalar->SetName(scalar_components[k].c_str());
    polyData->GetPointData()->AddArray(scalar);
  }

  return polyData;
}
#endif

inline void critical_point_tracker::write_sliced_critical_points_text(int t, std::ostream& os) const
{
  if (sliced_critical_points.find(t) == sliced_critical_points.end()) return;
  const auto &lcps = sliced_critical_points.at(t);

  for (const auto &lcp : lcps) {
    const auto &cp = std::get<0>(lcp);
    const int id = std::get<1>(lcp);
   
    os << "id=" << id << ", ";
    cp.print(os, cpdims(), scalar_components);
    os << std::endl;
  }
}

inline void critical_point_tracker::set_type_filter(unsigned int f)
{
  use_type_filter = true;
  type_filter = f;
}

inline void critical_point_tracker::set_scalar_components(const std::vector<std::string>& c) 
{
  assert(c.size() <= FTK_CP_MAX_NUM_VARS);
  scalar_components = c;
}

template <typename I>
void critical_point_tracker::grow_trajectories(
      std::vector<critical_point_traj_t> &trajectories,
      // std::vector<bool> &alive,
      std::map<I, critical_point_t> &discrete_critical_points, // will be cleared after tracing
      std::function<std::set<I>(I)> neighbors,
      std::function<I(unsigned long long)> tag_to_element)
{
  // 0. gather discrete trajectories
  diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, get_root_proc());
  if (!is_root_proc()) return;

  // 1. continue existing trajectories
  for (auto k = 0; k < trajectories.size(); k ++) {
    auto &traj = trajectories[k];
    if (traj.complete) continue;
    bool continued = false;

    // forward direction
    auto terminal = traj.back();
    auto current = tag_to_element(terminal.tag);
    while (1) {
      bool has_next = false;
      for (auto i : neighbors(current)) {
        if (discrete_critical_points.find(i) != discrete_critical_points.end()) 
        {
          current = i;
          const auto cp = discrete_critical_points[current];
          if (cp.ordinal)
            sliced_critical_points[cp.timestep].push_back(std::make_tuple(cp, k));
          traj.push_back(cp);
          discrete_critical_points.erase(current);
          has_next = true;
          continued = true;
          break;
        }
      }
      if (!has_next)
        break;
    }

    // backward direction
    terminal = traj.front();
    current = tag_to_element(terminal.tag);
    while (1) {
      bool has_next = false;
      for (auto i : neighbors(current)) {
        if (discrete_critical_points.find(i) != discrete_critical_points.end()) 
        {
          // fprintf(stderr, "tracing backwards!!\n");
          current = i;
          const auto cp = discrete_critical_points[current];
          if (cp.ordinal)
            sliced_critical_points[cp.timestep].push_back(std::make_tuple(cp, k));
          traj.insert(traj.begin(), cp); // TODO: improve performance by using list instead of vector
          discrete_critical_points.erase(current);
          has_next = true;
          continued = true;
          break;
        }
      }
      if (!has_next)
        break;
    }

    if (!continued) traj.complete = true;
  }

  // 2. generate new trajectories for the rest of discrete critical points
  std::set<I> elements;
  for (const auto &kv : discrete_critical_points)
    elements.insert(kv.first);
  auto connected_components = extract_connected_components<I, std::set<I>>(
      neighbors, elements);

  for (const auto &component : connected_components) 
  {
    auto linear_graphs = ftk::connected_component_to_linear_components<I>(component, neighbors);
    assert(linear_graphs.size() == 1);
    // fprintf(stderr, "size_component=%zu, size_linear_graph=%zu\n", component.size(), linear_graphs.size());
    for (int j = 0; j < linear_graphs.size(); j ++) {
      const auto &linear_graph = linear_graphs[j];
      critical_point_traj_t traj;
      traj.loop = is_loop(linear_graph, neighbors);

      const int new_id = trajectories.size();
      for (int k = 0; k < linear_graph.size(); k ++) {
        const auto &cp = discrete_critical_points[linear_graph[k]];
        traj.push_back(cp);
        if (cp.ordinal)
          sliced_critical_points[cp.timestep].push_back(std::make_tuple(cp, new_id));
      }
      trajectories.push_back(traj);
      // fprintf(stderr, "birth.\n");
    }
  }

  // 3. clear discrete critical points
  discrete_critical_points.clear();

  // write_sliced_critical_points_text(current_timestep, std::cerr);
}

inline void critical_point_tracker::update_traj_statistics()
{
  for (auto &traj : traced_critical_points)
    traj.update_statistics();
}

inline void critical_point_tracker::select_traj(std::function<bool(const critical_point_traj_t& traj)> f)
{
  std::vector<critical_point_traj_t> selected_traj;
  for (const auto &traj : traced_critical_points) 
    if (f(traj))
      selected_traj.push_back(traj);
  traced_critical_points = selected_traj;
}

template <typename element_t>
std::vector<critical_point_traj_t> critical_point_tracker::trace_critical_points_offline(
	const std::map<element_t, critical_point_t> &discrete_critical_points,
	std::function<std::set<element_t>(element_t)> neighbors)
{
  std::vector<critical_point_traj_t> traced_critical_points;

  std::set<element_t> elements;
  for (const auto &kv : discrete_critical_points)
    elements.insert(kv.first);
  auto connected_components = extract_connected_components<element_t, std::set<element_t>>(
      neighbors, elements);

  for (const auto &component : connected_components) {
    std::vector<std::vector<double>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(component, neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      critical_point_traj_t traj; 
      traj.loop = is_loop(linear_graphs[j], neighbors);
      for (int k = 0; k < linear_graphs[j].size(); k ++)
        traj.push_back(discrete_critical_points.at(linear_graphs[j][k]));
      traced_critical_points.emplace_back(traj);
    }
  }

  return traced_critical_points;
}

void critical_point_tracker::slice_traced_critical_points(unsigned int type_filter = 0xffffffff)
{
  for (auto i = 0; i < traced_critical_points.size(); i ++) {
    const auto &traj = traced_critical_points[i];
    for (auto j = 0; j < traj.size(); j ++) {
      const auto &cp = traj[j];
      if (cp.ordinal && (cp.type & type_filter)) //  && cp.type == CRITICAL_POINT_2D_MAXIMUM)
        sliced_critical_points[cp.timestep].push_back(
            std::make_tuple(cp, i));
    }
  }
}

}

#endif
