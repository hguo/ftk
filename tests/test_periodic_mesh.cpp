#include <ftk/mesh/simplicial_unstructured_periodic_2d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_2d_mesh.hh>
#include <ftk/ndarray.hh>

diy::mpi::environment env;

int main(int argc, char **argv)
{
  using namespace ftk;
  auto m2 = simplicial_xgc_2d_mesh<>::from_xgc_mesh_file(argv[1]);
  auto m3 = simplicial_unstructured_3d_mesh<>::from_file(argv[2]);

  std::shared_ptr<simplicial_unstructured_periodic_2d_mesh<>> m;
  m.reset(new simplicial_unstructured_periodic_2d_mesh<>(
        std::dynamic_pointer_cast<simplicial_unstructured_2d_mesh<>>(m2), m3));
  
  const int nphi = 3;

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(m2->n(0) * (nphi + 1));
  
  for (int i=0; i < m2->n(0) * (nphi + 1); i++) {
    double coords[4];
    m->get_coords(i, coords);
    pts->SetPoint(i, coords[0], coords[1], coords[3]);
  }

  auto f3 = [&](int k) {
    int verts[4];
    m->get_simplex(3, k, verts);
    vtkIdType ids[4] = {verts[0], verts[1], verts[2], verts[3]};
    // fprintf(stderr, "%d: %d, %d, %d, %d\n", k, verts[0], verts[1], verts[2], verts[3]);
    grid->InsertNextCell(VTK_TETRA, 4, ids);
  };
  
  auto f2 = [&](int k) {
    int verts[3];
    m->get_simplex(2, k, verts);
    vtkIdType ids[3] = {verts[0], verts[1], verts[2]};
    // fprintf(stderr, "%d: %d, %d, %d, %d\n", k, verts[0], verts[1], verts[2], verts[3]);
    grid->InsertNextCell(VTK_TRIANGLE, 3, ids);
  };
  
  auto f1 = [&](int k) {
    int verts[2];
    m->get_simplex(1, k, verts);
    vtkIdType ids[2] = {verts[0], verts[1]};
    grid->InsertNextCell(VTK_LINE, 2, ids);
  };

#if 1
  for (int p = 0; p < nphi; p ++) {
    // m->element_for_ordinal(3, p, f3, FTK_XL_NONE, 1);
    m->element_for_ordinal(2, p, f2, FTK_XL_NONE, 1);
    m->element_for_interval(2, p, f2, FTK_XL_NONE, 1);
  }
  m->element_for_ordinal(2, nphi, f2, FTK_XL_NONE, 1);
#endif
  // m->element_for_ordinal(2, 0, f2, FTK_XL_NONE, 1);
  
  grid->SetPoints(pts);
  // grid->PrintSelf(std::cerr, vtkIndent(2));
  
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName("out.vtu");
  writer->SetInputData( grid );
  writer->Write();
#endif

  return 0;
}
