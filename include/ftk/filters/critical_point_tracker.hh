#ifndef _FTK_CRITICAL_POINT_TRACKER_HH
#define _FTK_CRITICAL_POINT_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/features/feature_point.hh>
#include <ftk/features/feature_curve.hh>
#include <ftk/features/feature_curve_set.hh>
#include <ftk/basic/duf.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/write_polydata.hh>
#include <ftk/utils/gather.hh>
#include <ftk/utils/redistribution.hh>
#include <iomanip>

namespace ftk {

enum {
  SOURCE_NONE, 
  SOURCE_GIVEN, // explicit
  SOURCE_DERIVED // implicit
};

struct critical_point_tracker : public virtual tracker {
  critical_point_tracker(diy::mpi::communicator comm) : tracker(comm) {}

  virtual void update() {}; 
  void reset() {
    field_data_snapshots.clear();
    traced_critical_points.clear();
  }

  void set_enable_robust_detection(bool b) { enable_robust_detection = b; }
  void set_enable_computing_degrees(bool b) { enable_computing_degrees = b; }
  void set_enable_streaming_trajectories(bool b) { enable_streaming_trajectories = b; }
  void set_enable_discarding_interval_points(bool b) { enable_discarding_interval_points = b; }
  void set_enable_discarding_degenerate_points(bool b) { enable_discarding_degenerate_points = b; }
  void set_enable_ignoring_degenerate_points(bool b) { enable_ignoring_degenerate_points = b; }

  void set_type_filter(unsigned int);

  void set_scalar_field_source(int s) {scalar_field_source = s;}
  void set_vector_field_source(int s) {vector_field_source = s;}
  void set_jacobian_field_source(int s) {jacobian_field_source = s;}
  void set_jacobian_symmetric(bool s) {is_jacobian_field_symmetric = s;}

  void set_scalar_components(const std::vector<std::string>& c);
  int get_num_scalar_components() const {return scalar_components.size();}

  void update_traj_statistics();

  void select_trajectories(std::function<bool(const feature_curve_t& traj)>);
  void select_sliced_critical_points(std::function<bool(const feature_point_t& cp)>);

  void slice_traced_critical_points(); // slice traces after finalization

public:
  // virtual void initialize() = 0;
  // virtual void finalize() = 0;

  bool advance_timestep();
  // virtual void update_timestep() = 0;

public: // i/o for traced critical points (trajectories)
  const feature_curve_set_t& get_traced_critical_points() const {return traced_critical_points;}
  feature_curve_set_t& get_traced_critical_points() {return traced_critical_points;}

  json get_traced_critical_points_json() const {return json(traced_critical_points);}
  void write_traced_critical_points_json(const std::string& filename, int indent=0) const;
  void read_traced_critical_points_json(const std::string& filename);
  void write_traced_critical_points_binary(const std::string& filename) const;
  void read_traced_critical_points_binary(const std::string& filename);
  void write_traced_critical_points_text(std::ostream& os) const;
  void write_traced_critical_points_text(const std::string& filename) const;
  void write_traced_critical_points_vtk(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_traced_critical_points_vtk() const {return traced_critical_points.to_vtp(scalar_components);}
#endif

public: // i/o for sliced critical points
  const std::map<int, std::vector<feature_point_t>>& get_sliced_critical_points() const {return sliced_critical_points;}
  json get_sliced_critical_points_json() const {return json(sliced_critical_points);}
  void write_sliced_critical_points_json(const std::string& filename) const {std::ofstream f(filename); f << get_sliced_critical_points_json(); f.close();}
  void read_sliced_critical_points_json(const std::string& filename) {std::ifstream f(filename); json j; f >> j; f.close(); sliced_critical_points = j.get<std::map<int, std::vector<feature_point_t>>>();}
  void write_sliced_critical_points_text(int t, std::ostream& os) const;
  void write_sliced_critical_points_text(int t, const std::string& filename) const;
  void write_sliced_critical_points_vtk(int t, const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_sliced_critical_points_vtk(int t) const;
#endif

public: // i/o for intercepted traced
  feature_curve_set_t get_intercepted_critical_point(int t0, int t1) const {return traced_critical_points.intercept(t0, t1);}
  void write_intercepted_critical_points_vtk(int t0, int t1, const std::string& filename) const;
  void write_intercepted_critical_points_text(int t0, int t1, const std::string& filename) const;
  void write_intercepted_critical_points_json(int t0, int t1, const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_intercepted_critical_points_vtk(int t0, int t1) const {return get_intercepted_critical_point(t0, t1).to_vtp(scalar_components);}
#endif

public: // i/o for discrete (untraced) critical points
  virtual std::vector<feature_point_t> get_critical_points() const = 0;
  json get_critical_points_json() const {return json(get_critical_points());}
  virtual void put_critical_points(const std::vector<feature_point_t>&) = 0;
  void write_critical_points_json(const std::string& filename) const {std::ofstream f(filename); f << get_critical_points_json(); f.close();}
  void read_critical_points_json(const std::string& filename);
  void write_critical_points_binary(const std::string& filename) const;
  void read_critical_points_binary(const std::string& filename);
  void write_critical_points_text(std::ostream& os) const;
  void write_critical_points_text(const std::string& filename) const;
  void write_critical_points_vtk(const std::string& filename) const;
  void write_critical_points(const std::string& filename) const; // automatically determine format
  void read_critical_points(const std::string& filename);
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_critical_points_vtk() const;
#endif

public: // post-processing and simplification
  // void foreach_trajectory(std::function<void(int, feature_curve_t&)> f) {for (int i = 0; i < traced_critical_points.size(); i ++) f(i, traced_critical_points[i]);}
  // void split_trajectories();

public: // inputs
  bool pop_field_data_snapshot();
  virtual void push_field_data_snapshot(
      const ndarray<double> &scalar, 
      const ndarray<double> &vector,
      const ndarray<double> &jacobian);
  virtual void push_scalar_field_snapshot(const ndarray<double> &scalar);
  virtual void push_vector_field_snapshot(const ndarray<double> &vector);

protected:
  bool filter_critical_point_type(const feature_point_t& cp);

  void update_vector_field_scaling_factor(int minbits=8, int maxbits=21);

protected:
  template <typename I> // mesh element type
  void trace_critical_points_online(
      feature_curve_set_t &trajectories,
      std::map<I, feature_point_t> &discrete_critical_poionts, // critical point tag needs to index mesh element ID.  Discrete critical points will be cleared after tracing
      std::function<std::set<I>(I)> neighbors, 
      std::function<I(unsigned long long)> tag_to_element,
      std::function<unsigned long long(I)> elementt_to_tag
  );

	template <typename I> // mesh element type
	std::vector<feature_curve_t> trace_critical_points_offline(
		std::map<I, feature_point_t> &discrete_critical_points, // id of each cp will be updated
		std::function<std::set<I>(I)> neighbors);

protected:
  struct field_data_snapshot_t {
    ndarray<double> scalar, vector, jacobian;
  };

  std::deque<field_data_snapshot_t> field_data_snapshots;
  
  // for robust detection
  double vector_field_resolution = std::numeric_limits<double>::max(); // min abs nonzero value of vector field.  for robust cp detection w/o gmp
  uint64_t vector_field_scaling_factor = 1;
  
  feature_curve_set_t traced_critical_points;
  std::map<int/*time*/, std::vector<feature_point_t>> sliced_critical_points;

  // type filter
  bool use_type_filter = false;
  unsigned int type_filter = 0;
  
  int scalar_field_source = SOURCE_NONE, 
      vector_field_source = SOURCE_NONE,
      jacobian_field_source = SOURCE_NONE;
  bool is_jacobian_field_symmetric = false;

  // scalar components
  std::vector<std::string> scalar_components = {"scalar"};

  bool enable_robust_detection = true;
  bool enable_computing_degrees = false;
  bool enable_streaming_trajectories = false;
  bool enable_discarding_interval_points = false;
  bool enable_discarding_degenerate_points = false;
  bool enable_ignoring_degenerate_points = false;
};

///////

inline bool critical_point_tracker::filter_critical_point_type(
    const feature_point_t& cp)
{
  // fprintf(stderr, "typefilter=%lu, type=%lu\n", 
  //     type_filter, cp.type);
  if (use_type_filter) {
    if (type_filter & cp.type) return true;
    else return false;
  }
  else return true;
}

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

inline void critical_point_tracker::push_vector_field_snapshot(const ndarray<double>& vector)
{
  field_data_snapshot_t snapshot;
  snapshot.vector = vector;

  field_data_snapshots.emplace_back(snapshot);
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
inline void critical_point_tracker::write_intercepted_critical_points_vtk(int t0, int t1, const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_intercepted_critical_points_vtk(t0, t1);
    write_polydata(filename, poly);
  }
}

inline void critical_point_tracker::write_traced_critical_points_vtk(const std::string& filename) const
{
  auto poly = traced_critical_points.to_vtp(scalar_components); 
  write_polydata(filename, poly, "auto", comm);
}

inline void critical_point_tracker::write_sliced_critical_points_vtk(int k, const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_sliced_critical_points_vtk(k);
    write_polydata(filename, poly);
  }
}

inline void critical_point_tracker::write_critical_points_vtk(const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_critical_points_vtk();
    write_polydata(filename, poly);
  }
}

inline vtkSmartPointer<vtkPolyData> critical_point_tracker::get_critical_points_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  
  vtkIdType pid[1];
  
  // const auto critical_points = get_critical_points();
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

inline void critical_point_tracker::write_intercepted_critical_points_vtk(int t0, int t1, const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

inline void critical_point_tracker::write_critical_points_binary(const std::string& filename) const 
{
  if (is_root_proc())
    diy::serializeToFile(get_critical_points(), filename);
}

inline void critical_point_tracker::read_critical_points_binary(const std::string& filename)
{
  if (is_root_proc()) {
    std::vector<feature_point_t> cps;
    diy::unserializeFromFile(filename, cps);
    put_critical_points(cps);
  }
}

inline void critical_point_tracker::write_traced_critical_points_binary(const std::string& filename) const
{
  if (is_root_proc()) 
    diy::serializeToFile(traced_critical_points, filename);
}

inline void critical_point_tracker::read_traced_critical_points_binary(const std::string& filename)
{
  if (is_root_proc()) 
    diy::unserializeFromFile(filename, traced_critical_points);
}

inline void critical_point_tracker::write_traced_critical_points_text(const std::string& filename) const
{
  if (is_root_proc()) {
    std::ofstream out(filename);
    traced_critical_points.write_text(out, scalar_components);
    out.close();
  }
}

inline void critical_point_tracker::write_traced_critical_points_text(std::ostream& os) const
{
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
    cp.print(os, scalar_components) << std::endl;
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> critical_point_tracker::get_sliced_critical_points_vtk(int t) const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();

  if (sliced_critical_points.find(t) == sliced_critical_points.end()) 
    return polyData;
  const auto &cps = sliced_critical_points.at(t);
  
  vtkIdType pid[1];
  for (const auto &cp : cps) {
    double p[3] = {cp[0], cp[1], cp[2]};
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
  
  // types & ids
  vtkSmartPointer<vtkIntArray> ids = vtkSmartPointer<vtkIntArray>::New();
  ids->SetName("id");
  ids->SetNumberOfValues(cps.size());
  
  vtkSmartPointer<vtkUnsignedIntArray> types = vtkSmartPointer<vtkUnsignedIntArray>::New();
  types->SetName("type");
  types->SetNumberOfValues(cps.size());
  
  vtkSmartPointer<vtkDoubleArray> velocities = vtkSmartPointer<vtkDoubleArray>::New();
  velocities->SetName("velocity");
  velocities->SetNumberOfComponents(3);
  velocities->SetNumberOfTuples(cps.size());

  int i = 0;
  for (const auto &cp : cps) {
    ids->SetValue(i, cp.id);
    types->SetValue(i, cp.type);
    velocities->SetTuple3(i, cp.v[0], cp.v[1], cp.v[2]);
    i ++;
  }

  polyData->GetPointData()->AddArray(ids);
  polyData->GetPointData()->AddArray(types);
  polyData->GetPointData()->AddArray(velocities);
  polyData->GetPointData()->SetActiveScalars("id");
  
  for (auto k = 0; k < scalar_components.size(); k ++) {
    vtkSmartPointer<vtkDoubleArray> scalar = vtkSmartPointer<vtkDoubleArray>::New();
    scalar->SetNumberOfValues(cps.size());
    int i = 0;
    for (const auto &cp : cps)
      scalar->SetValue(i ++, cp.scalar[k]);
    scalar->SetName(scalar_components[k].c_str());
    polyData->GetPointData()->AddArray(scalar);
  }

  return polyData;
}
#endif

inline void critical_point_tracker::write_sliced_critical_points_text(int t, std::ostream& os) const
{
  if (sliced_critical_points.find(t) == sliced_critical_points.end()) return;
  const auto &cps = sliced_critical_points.at(t);

  for (const auto &cp: cps) {
    cp.print(os, scalar_components);
    os << std::endl;
  }
}

inline void critical_point_tracker::write_traced_critical_points_json(const std::string& filename, int indent) const 
{
  if (is_root_proc()) {
    std::ofstream f(filename); 
    f << std::setw(indent) 
      << get_traced_critical_points_json(); 
    f.close();
  }
}

inline void critical_point_tracker::read_traced_critical_points_json(const std::string& filename) 
{
  if (is_root_proc()) {
    std::ifstream f(filename); 
    json j; 
    f >> j; 
    f.close(); 
    // traced_critical_points = j.get<std::map<int, feature_curve_t>>();
    traced_critical_points = j.get<feature_curve_set_t>();
  }
}

inline void critical_point_tracker::read_critical_points_json(const std::string& filename)
{
  if (is_root_proc()) {
    std::ifstream f(filename); 
    json j; 
    f >> j; 
    // fprintf(stderr, "#cps=%zu\n", j.size()); 
    f.close(); 
    put_critical_points(j.get<std::vector<feature_point_t>>());
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
void critical_point_tracker::trace_critical_points_online(
      feature_curve_set_t &trajectories,
      // std::vector<bool> &alive,
      std::map<I, feature_point_t> &discrete_critical_points, // will be cleared after tracing
      std::function<std::set<I>(I)> neighbors,
      std::function<I(unsigned long long)> tag_to_element,
      std::function<unsigned long long(I)> element_to_tag)
{
  // 0. gather discrete trajectories
  diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, get_root_proc());
  if (!is_root_proc()) return;

  // fprintf(stderr, "hey!!!, #traj=%zu, discrete points:\n", trajectories.size());
  // for (const auto kv : discrete_critical_points)
  //   std::cerr << "----" << kv.second.tag << ", " << kv.first << ", " << tag_to_element(kv.second.tag) << std::endl;

  // 1. continue existing trajectories
  trajectories.foreach([&](feature_curve_t& traj) {
    if (traj.complete) {
      // fprintf(stderr, "traj already complete.\n");
      return; // continue;
    }
    bool continued = false;

    // forward direction
    auto terminal = traj.back();
    auto current = tag_to_element(terminal.tag);
    // std::cerr << "forward terminal: " << terminal.tag << ", " << tag_to_element(terminal.tag) << std::endl;
    while (1) {
      bool has_next = false;
      for (auto i : neighbors(current)) {
        // std::cerr << "-------forward neighbor: " << element_to_tag(i) << ", " << i << std::endl;
        if (discrete_critical_points.find(i) != discrete_critical_points.end()) 
        {
          current = i;
          const auto cp = discrete_critical_points[current];
          // if (cp.ordinal) // TODO FIXME
          //   // sliced_critical_points[cp.timestep][k] = cp;
          //   sliced_critical_points[cp.timestep].push_back(cp);
          traj.push_back(cp);
          discrete_critical_points.erase(current);
          has_next = true;
          continued = true;
          // fprintf(stderr, "found next forward!!\n");
          break;
        }
      }
      if (!has_next)
        break;
    }

    // backward direction
    terminal = traj.front();
    current = tag_to_element(terminal.tag);
    // std::cerr << "backward terminal: " << terminal.tag << ", " << tag_to_element(terminal.tag) << std::endl;
    while (1) {
      bool has_next = false;
      for (auto i : neighbors(current)) {
        // std::cerr << "-------backward neighbor: " << element_to_tag(i) << ", " << i << std::endl;
        if (discrete_critical_points.find(i) != discrete_critical_points.end()) 
        {
          // fprintf(stderr, "tracing backwards!!\n");
          current = i;
          const auto cp = discrete_critical_points[current];
          // TODO FIXME
          // if (cp.ordinal)
          //   // sliced_critical_points[cp.timestep][k] = cp;
          //   sliced_critical_points[cp.timestep].push_back(cp);
          traj.insert(traj.begin(), cp); // TODO: improve performance by using list instead of vector
          discrete_critical_points.erase(current);
          has_next = true;
          continued = true;
          // fprintf(stderr, "found next backward!!\n");
          break;
        }
      }
      if (!has_next)
        break;
    }

    if (!continued) traj.complete = true;
  });

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
      feature_curve_t traj;
      traj.loop = is_loop(linear_graph, neighbors);

      // const int new_id = trajectories.size();
      for (int k = 0; k < linear_graph.size(); k ++) {
        const auto &cp = discrete_critical_points[linear_graph[k]];
        traj.push_back(cp);
        // if (cp.ordinal) // TODO FIXME
        //   sliced_critical_points[cp.timestep].push_back(cp);
        //   // sliced_critical_points[cp.timestep][new_id] = cp;
      }
      // trajectories.push_back(traj);
      trajectories.add(traj);
      // fprintf(stderr, "birth.\n");
    }
  }

  // 3. clear discrete critical points
  discrete_critical_points.clear();

  // write_sliced_critical_points_text(current_timestep, std::cerr);
}

inline void critical_point_tracker::update_traj_statistics()
{
  traced_critical_points.foreach([](feature_curve_t& t) {
      t.update_statistics();
  });
}

inline void critical_point_tracker::select_trajectories(std::function<bool(const feature_curve_t& traj)> f)
{
  traced_critical_points.filter(f);
}

inline void critical_point_tracker::select_sliced_critical_points(std::function<bool(const feature_point_t& cp)> f)
{
  for (auto &kv : sliced_critical_points) {
    auto &cps = kv.second;
    std::vector<feature_point_t> filtered_cps;
    for (auto i = 0; i < cps.size(); i ++) 
      if (f(cps[i]))
        filtered_cps.push_back(cps[i]);
    
    cps = filtered_cps;
  }
}

template <typename element_t>
std::vector<feature_curve_t> critical_point_tracker::trace_critical_points_offline(
	std::map<element_t, feature_point_t> &discrete_critical_points,
	std::function<std::set<element_t>(element_t)> neighbors)
{
  std::vector<feature_curve_t> traced_critical_points;
 
#if 0
	std::map<element_t, feature_point_t> all_discrete_critical_points;
  // diy::mpi::gather(comm, discrete_critical_points, all_discrete_critical_points, get_root_proc());
  diy::mpi::gather<std::map<element_t, feature_point_t>>(comm, discrete_critical_points, all_discrete_critical_points, get_root_proc());
  int ncps = discrete_critical_points.size(), nallcps;
  diy::mpi::all_reduce(comm, ncps, nallcps, std::plus<int>());
  fprintf(stderr, "#cps=%d, #sumcps=%d, #allcps=%zu\n", ncps, nallcps, all_discrete_critical_points.size());
  comm.barrier();
  // exit(1);
#endif

#if 1
  // fprintf(stderr, "dUF..\n");
  duf<element_t> uf(comm);
  diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, get_root_proc()); // TODO FIXME
  if (!is_root_proc()) discrete_critical_points.clear();
  
  for (const auto &kv : discrete_critical_points) {
    for (const auto &n : neighbors(kv.first))
      // if (true) // (discrete_critical_points.find(n) != discrete_critical_points.end())
      if (discrete_critical_points.find(n) != discrete_critical_points.end())
        uf.unite(kv.first, n);
  }
  // fprintf(stderr, "dUF sync...\n");
  uf.sync();
  // fprintf(stderr, "dUF done.\n");
  // fprintf(stderr, "dUF done., #pts=%zu, #roots=%zu\n", discrete_critical_points.size(), uf.get_roots().size());

  std::map<element_t/*root*/, std::map<element_t, feature_point_t>> ccs, rccs; // distributed cc
  for (const auto &kv : discrete_critical_points)
    ccs[ uf.find(kv.first) ].insert(kv);

  redistribute(comm, ccs, rccs);
  
  std::mutex my_mutex;
  // for (auto &cc : rccs) { 
  object::parallel_for_container<std::map<element_t, std::map<element_t, feature_point_t>>>
    (rccs, [&](typename std::map<element_t, std::map<element_t, feature_point_t>>::iterator icc) {
    std::set<element_t> component;
    for (const auto &kv : icc->second)
      component.insert(kv.first);

    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(component, neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      const unsigned int id = traced_critical_points.size(); 

      feature_curve_t traj; 
      traj.loop = is_loop(linear_graphs[j], neighbors);
      for (int k = 0; k < linear_graphs[j].size(); k ++) {
        // auto &cp = discrete_critical_points[linear_graphs[j][k]];
        auto cp = icc->second.at(linear_graphs[j][k]);
        cp.id = id;

        traj.push_back(cp);
        // sum1 ++;
      }
      
      {
        std::lock_guard<std::mutex> guard(my_mutex);
        traced_critical_points.emplace_back(traj);
      }
    }
  });

  fprintf(stderr, "rank=%d, #curves=%zu\n", comm.rank(), traced_critical_points.size());

  // FIXME: need to know earlier if the output trajectories will be written (into pvtp) later.  In this case, no gathering is needed
  diy::mpi::gather<std::vector<feature_curve_t>>(comm, traced_critical_points, traced_critical_points, get_root_proc(), 
      [](const std::vector<feature_curve_t>& in, std::vector<feature_curve_t>& out) {
        for (const auto& v : in)
          out.push_back(v);
      });

  if (comm.rank() == 0) 
    fprintf(stderr, "total curves: %zu\n", traced_critical_points.size());
#endif

#if 0
  fprintf(stderr, "computing cc..\n");
  std::set<element_t> elements;
  for (const auto &kv : discrete_critical_points)
    elements.insert(kv.first);
  auto connected_components = extract_connected_components<element_t, std::set<element_t>>(
      neighbors, elements);

  fprintf(stderr, "#cc=%zu\n", connected_components.size());
  // int sum = 0;
  // for (auto c : connected_components)
  //   sum += c.size();
  // fprintf(stderr, "#sum=%d, #dcp=%zu\n", sum, discrete_critical_points.size());

  // int sum1 = 0;
  std::mutex my_mutex;
  object::parallel_for(connected_components.size(), [&](int cid) {
    const auto &component = connected_components[cid];
    
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(component, neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      const unsigned int id = traced_critical_points.size(); 

      feature_curve_t traj; 
      traj.loop = is_loop(linear_graphs[j], neighbors);
      for (int k = 0; k < linear_graphs[j].size(); k ++) {
        auto &cp = discrete_critical_points[linear_graphs[j][k]];
        cp.id = id;

        traj.push_back(cp);
        // sum1 ++;
      }
      
      {
        std::lock_guard<std::mutex> guard(my_mutex);
        traced_critical_points.emplace_back(traj);
      }
    }
  });
#endif

#if 0
  for (const auto &component : connected_components) {
    std::vector<std::vector<double>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(component, neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      const unsigned int id = traced_critical_points.size(); 

      feature_curve_t traj; 
      traj.loop = is_loop(linear_graphs[j], neighbors);
      for (int k = 0; k < linear_graphs[j].size(); k ++) {
        auto &cp = discrete_critical_points[linear_graphs[j][k]];
        cp.id = id;

        traj.push_back(cp);
        // sum1 ++;
      }
      traced_critical_points.emplace_back(traj);
    }
  }
#endif

  // fprintf(stderr, "#sum1=%d\n", sum1);

  return traced_critical_points;
}

inline void critical_point_tracker::slice_traced_critical_points()
{
  int sum0 = 0;
  for (const auto &kv : traced_critical_points) {
    const auto &traj = kv.second;
    for (auto j = 0; j < traj.size(); j ++) {
      const auto &cp = traj[j];
      if (cp.ordinal) {
        sliced_critical_points[cp.timestep].push_back(cp); // [i] = cp;
        sum0 ++;
      }
    }
  }

  int sum = 0;
  for (const auto &kv : sliced_critical_points) {
    fprintf(stderr, "timestep=%d, n=%zu\n", kv.first, kv.second.size());
    sum += kv.second.size();
  }
  fprintf(stderr, "#sum_sliced=%d, %d\n", sum, sum0);
}

inline bool critical_point_tracker::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}
  
inline void critical_point_tracker::update_vector_field_scaling_factor(int minbits, int maxbits)
{
  // vector_field_resolution = std::numeric_limits<double>::max();
  for (const auto &s : field_data_snapshots)
    vector_field_resolution = std::min(vector_field_resolution, s.vector.resolution());
  
  int nbits = std::ceil(std::log2(1.0 / vector_field_resolution));
  nbits = std::max(minbits, std::min(nbits, maxbits));

  vector_field_scaling_factor = 1 << nbits;
 
  std::cerr << "resolution=" << vector_field_resolution 
    << ", factor=" << vector_field_scaling_factor 
    << ", nbits=" << nbits << std::endl;
}

}

#endif
