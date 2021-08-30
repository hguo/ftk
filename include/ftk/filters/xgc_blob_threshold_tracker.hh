#ifndef _FTK_XGC_BLOB_CCL_TRACKER_HH
#define _FTK_XGC_BLOB_CCL_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/filters/xgc_tracker.hh>
#include <ftk/tracking_graph/tracking_graph.hh>
#include <ftk/basic/duf.hh>

#if FTK_HAVE_VTK
#include <vtkThreshold.h>
#endif

namespace ftk {

struct xgc_blob_threshold_tracker : public xgc_tracker {
  xgc_blob_threshold_tracker(diy::mpi::communicator comm, 
      std::shared_ptr<simplicial_xgc_3d_mesh<>> m3) : xgc_tracker(comm, m3) {}
  virtual ~xgc_blob_threshold_tracker() {}
  
  void set_threshold(double t) { threshold = t; }

  // int cpdims() const { return 0; }
  
  void initialize() {}
  void update() {}
  void finalize();

  void update_timestep();

  void push_field_data_snapshot(const ndarray<double> &scalar);

public:
  ndarray<int> get_sliced(int t) const;
  void write_sliced(int t, const std::string& pattern) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> sliced_to_vtu_slices(int t) const;
  vtkSmartPointer<vtkUnstructuredGrid> sliced_to_vtu_solid(int t) const;
  // vtkSmartPointer<vtkUnstructuredGrid> sliced_to_vtu_partial_solid(int t) const;
#endif

protected:
  double threshold = 0.0;
  duf<int> uf;
};

/////

inline void xgc_blob_threshold_tracker::push_field_data_snapshot(const ndarray<double> &scalar)
{
  ndarray<double> grad, J; // no grad or jacobian needed
  
  fprintf(stderr, "threshold_2.5_sigma=%f\n", derive_threshold(scalar));

  if (m3->has_smoothing_kernel()) {
    ndarray<double> smoothed_scalar = m3->smooth_scalar(scalar);
    xgc_tracker::push_field_data_snapshot(smoothed_scalar, grad, J);
  } else 
    xgc_tracker::push_field_data_snapshot(scalar, grad, J);
}

inline void xgc_blob_threshold_tracker::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  const int m3n0 = m3->n(0);
  const auto &scalar = field_data_snapshots[0].scalar;
  // for (int i = 0; i < m3n0; i ++) {
  parallel_for(m3n0, [&](int i) {
    if (m3->interpolate(scalar, i) < threshold) return; // continue;
    const int ei = i + current_timestep * m3n0;
    // for (const auto j : m3->get_vertex_edge_vertex_nextnodes(i)) {
    for (const auto j : m3->get_vertex_edge_vertex(i)) {
      if (m3->interpolate(scalar, j) >= threshold) {
        const int ej = j + current_timestep * m3n0;
        uf.unite(ei, ej);
      }
    }
  });

  if (field_data_snapshots.size() >= 2) {
    const auto &scalar1 = field_data_snapshots[1].scalar;
    // for (int i = 0; i < m3n0; i ++) {
    parallel_for(m3n0, [&](int i) {
      // if (scalar[i] < threshold || scalar1[i] < threshold) continue;
      if (m3->interpolate(scalar, i) < threshold || m3->interpolate(scalar1, i) < threshold) return; // continue;
      else
        uf.unite(i + current_timestep * m3n0, i + (current_timestep+1) * m3n0);
    });
  }
  
  fprintf(stderr, "%zu\n", uf.size());
}

inline void xgc_blob_threshold_tracker::finalize()
{
}

void xgc_blob_threshold_tracker::write_sliced(int t, const std::string& pattern) const
{
  const std::string filename = series_filename(pattern, t);

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  // writer->SetInputData( sliced_to_vtu_slices(t) );
  writer->SetInputData( sliced_to_vtu_solid(t) );
  writer->Write();
#else
  fatal("FTK not compiled with VTK.");
#endif
}

ndarray<int> xgc_blob_threshold_tracker::get_sliced(int t) const
{
  ndarray<int> array({m3->n(0)}, -1);
  for (int i = 0; i < m3->n(0); i ++) {
    const int e = i + t * m3->n(0);
    if (uf.has(e)) 
      array[i] = uf.find(e);
  }
  return array;
}

#if FTK_HAVE_VTK
vtkSmartPointer<vtkUnstructuredGrid> xgc_blob_threshold_tracker::sliced_to_vtu_slices(int t) const
{
  auto grid = m3->to_vtu_slices();
  grid->GetPointData()->AddArray( get_sliced(t).to_vtk_data_array("id") );
  return grid;
}

vtkSmartPointer<vtkUnstructuredGrid> xgc_blob_threshold_tracker::sliced_to_vtu_solid(int t) const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();

  std::map<int, vtkIdType> idmap;
  for (int i = 0; i < m3->n(0); i ++) {
    const int e = i + t * m3->n(0);
    if (uf.has(e)) {
      double coords[3];
      m3->get_coords(i, coords);
      vtkIdType id = points->InsertNextPoint(coords[0], coords[1], coords[2]);
      idmap[i] = id;
    }
    // array[i] = uf.find(e);
  }
  grid->SetPoints(points);

  std::mutex mutex;
  for (int tid = 0; tid < m3->n(3); tid ++) {
  // parallel_for(m3->n(3), [&](int tid) {
    int tet[4];
    m3->get_simplex(3, tid, tet);

    std::vector<vtkIdType> ids;
    for (int k = 0; k < 4; k ++)
      if (idmap.find(tet[k]) != idmap.end())
        ids.push_back(idmap[tet[k]]);
      
    std::sort(ids.begin(), ids.end());
    {
      // std::lock_guard<std::mutex> guard(mutex);
      if (ids.size() == 4)
        grid->InsertNextCell(VTK_TETRA, 4, &ids[0]);
      else if (ids.size() == 3)
        grid->InsertNextCell(VTK_TRIANGLE, 3, &ids[0]);
      // else if (ids.size() == 2)
      //   grid->InsertNextCell(VTK_LINE, 2, &ids[0]);
    }
  }// );
  
  return grid;

#if 0
  vtkSmartPointer<vtkUnstructuredGrid> grid;
  if (mf3) grid = mf3->to_vtu_solid();
  else grid = m3->to_vtu_solid();

  auto ids = get_sliced(t).to_vtk_data_array("id");
  grid->GetPointData()->AddArray( ids );

  // auto scalar = field_data_snapshots[0].scalar.to_vtk_data_array("scalar");
  // grid->GetPointData()->AddArray( scalar );

  // grid->PrintSelf(std::cerr, vtkIndent(2));
  
  return grid;
#endif
}

// vtkSmartPointer<vtkUnstructuredGrid> xgc_blob_threshold_tracker::sliced_to_vtu_partial_solid(int t) const
// {
// }
#endif

}

#endif
