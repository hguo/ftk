#ifndef _FTK_XGC_BLOB_CCL_TRACKER_HH
#define _FTK_XGC_BLOB_CCL_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/filters/xgc_tracker.hh>
#include <ftk/tracking_graph/tracking_graph.hh>
#include <ftk/basic/duf.hh>

namespace ftk {

struct xgc_blob_threshold_tracker : public xgc_tracker {
  xgc_blob_threshold_tracker(diy::mpi::communicator comm, 
      std::shared_ptr<simplicial_xgc_3d_mesh<>> m3) : xgc_tracker(comm, m3) {}
  virtual ~xgc_blob_threshold_tracker() {}
  
  void set_threshold(double t) { threshold = t; }

  int cpdims() const { return 0; }
  
  void initialize() {}
  void update() {}
  void finalize();

  void update_timestep();

  void push_field_data_snapshot(const ndarray<double> &scalar);

public:
  void write_sliced(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const;
#endif

protected:
  double threshold = 0.0;
  duf<int> uf;
};

/////

inline void xgc_blob_threshold_tracker::push_field_data_snapshot(const ndarray<double> &scalar)
{
  ndarray<double> grad, J; // no grad or jacobian needed
  xgc_tracker::push_field_data_snapshot(scalar, grad, J);
}

inline void xgc_blob_threshold_tracker::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  const auto &scalar = field_data_snapshots[0].scalar;
  for (int i = 0; i < m3->n(0); i ++) {
    if (scalar[i] < threshold) continue;

    for (const auto j : m3->get_vertex_edge_vertex_nextnodes(i)) {
      if (scalar[j] >= threshold) 
        uf.unite(i, j);
    }
  }
  // m3->element_for(0, current_timestep, func, xl, nthreads, enable_set_affinity);
  
  fprintf(stderr, "%zu\n", uf.size());
}

inline void xgc_blob_threshold_tracker::finalize()
{
}

void xgc_blob_threshold_tracker::write_sliced(const std::string& filename) const
{
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( to_vtu() );
  writer->Write();
}

#if FTK_HAVE_VTK
vtkSmartPointer<vtkUnstructuredGrid> xgc_blob_threshold_tracker::to_vtu() const
{
  ndarray<int> array({m3->n(0)}, -1);
  for (int i = 0; i < m3->n(0); i ++) {
    if (uf.exists(i)) 
      array[i] = uf.find(i);
  }

  auto grid = m3->to_vtu_slices();
  grid->GetPointData()->AddArray( array.to_vtk_data_array() );

  return grid;
}
#endif

}

#endif
