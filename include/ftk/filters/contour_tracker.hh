#ifndef _FTK_CONTOUR_TRACKER_HH
#define _FTK_CONTOUR_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/algorithms/cca.hh>
// #include <ftk/filters/contour.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/features/feature_point.hh>
#include <ftk/features/feature_surface.hh>
#include <ftk/features/feature_volume.hh>
#include <ftk/geometry/write_polydata.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/utils/gather.hh>
#include <iomanip>

namespace ftk {

/**
 * @brief Contour/isosurface tracker for scalar fields
 *
 * This class tracks contours (isosurfaces in 3D, isolines in 2D) of scalar
 * fields over time. It detects where the scalar field equals a threshold
 * value and tracks how these contours evolve, move, merge, or split across
 * timesteps.
 *
 * The tracker can:
 * - Extract isosurfaces at a specified threshold value
 * - Track intersections with mesh elements
 * - Output results in VTK formats for visualization
 * - Process time-varying scalar fields
 *
 * Usage example:
 * @code
 * contour_tracker tracker(comm);
 * tracker.set_threshold(0.5);
 * tracker.initialize();
 * for (int t = 0; t < num_timesteps; t++) {
 *   tracker.push_field_data_snapshot(scalar);
 *   tracker.update_timestep();
 *   tracker.advance_timestep();
 * }
 * tracker.finalize();
 * tracker.write_isovolume_vtu("contours.vtu");
 * @endcode
 */
struct contour_tracker : public virtual tracker {
  contour_tracker(diy::mpi::communicator comm) : tracker(comm) {}

  void update() override {};
  void reset() override {
    field_data_snapshots.clear();
    // traced_contours.clear();
  }

  /**
   * @brief Set names for scalar field components
   * @param c Vector of component names
   */
  void set_scalar_components(const std::vector<std::string>& c);

  /**
   * @brief Get the number of scalar components
   * @return Number of scalar field components
   */
  int get_num_scalar_components() const {return scalar_components.size();}

  /**
   * @brief Get the current threshold value
   * @return Threshold value for contour extraction
   */
  double get_threshold() const { return threshold; }

  /**
   * @brief Set the threshold value for contour extraction
   * @param t Threshold value (scalar field value to track)
   */
  void set_threshold(double t) {threshold = t;}

public:
  bool advance_timestep() override;

public: // inputs
  /**
   * @brief Remove the oldest field data snapshot
   * @return True if a snapshot was removed
   */
  bool pop_field_data_snapshot() override;

  /**
   * @brief Push a scalar field snapshot for the next timestep
   * @param scalar Scalar field data
   */
  virtual void push_field_data_snapshot(const ndarray<double> &scalar);

public:
  /**
   * @brief Get contour intersections with mesh elements
   * @return Vector of feature points representing contour intersections
   */
  virtual std::vector<feature_point_t> get_intersections() const = 0;

  /**
   * @brief Write contour intersections to file
   * @param filenames Output filename
   */
  void write_intersections(const std::string& filenames) const;

  /**
   * @brief Write contour intersections as VTK polydata
   * @param filenames Output .vtp filename
   */
  void write_intersections_vtp(const std::string& filenames) const;
#if FTK_HAVE_VTK
  /**
   * @brief Get contour intersections as VTK polydata
   * @return VTK polydata representation of intersections
   */
  vtkSmartPointer<vtkPolyData> get_intersections_vtp() const;
#endif

  /**
   * @brief Write isovolume as VTK unstructured grid
   * @param filename Output .vtu filename
   */
  virtual void write_isovolume_vtu(const std::string& filename) const = 0;

  /**
   * @brief Write time-sliced contours as VTK unstructured grids
   * @param pattern Filename pattern with timestep placeholder
   */
  virtual void write_sliced_vtu(const std::string& pattern) const {}

  /**
   * @brief Write time-sliced contours as VTK polydata
   * @param pattern Filename pattern with timestep placeholder
   */
  virtual void write_sliced_vtp(const std::string& pattern) const {}

protected:
  // virtual int cpdims() const = 0;

protected:
  struct field_data_snapshot_t {
    ndarray<double> scalar, gradient;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;
  
  double threshold = 0.0;
  
  // scalar components
  std::vector<std::string> scalar_components = {"scalar"};
};

///////

inline void contour_tracker::push_field_data_snapshot(const ndarray<double>& scalar)
{
  field_data_snapshot_t snapshot;
  snapshot.scalar = scalar;
  if (scalar.nd() == 2) 
    snapshot.gradient = gradient2D(scalar);
  else 
    snapshot.gradient = gradient3D(scalar);

  field_data_snapshots.emplace_back(snapshot);
}

inline bool contour_tracker::pop_field_data_snapshot()
{
  if (field_data_snapshots.size() > 0) {
    field_data_snapshots.pop_front();
    return true;
  } else return false;
}

inline bool contour_tracker::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> contour_tracker::get_intersections_vtp() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
 
  vtkIdType pid[1];
  for (const auto &cp : get_intersections()) {
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]}; 
    // if (cpdims() == 2) p[2] = cp.t;
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
  
  vtkSmartPointer<vtkDoubleArray> time_array = vtkSmartPointer<vtkDoubleArray>::New();
  time_array->SetNumberOfValues(get_intersections().size());
  int i = 0;
  for (const auto &cp : get_intersections())
    time_array->SetValue(i++, cp.t);
  time_array->SetName("time");
  polyData->GetPointData()->AddArray(time_array);

  return polyData;
}

inline void contour_tracker::write_intersections_vtp(const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_intersections_vtp();
    write_polydata(filename, poly);
  }
}
#else
inline void contour_tracker::write_intersections_vtp(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

}

#endif
