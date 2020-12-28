#ifndef _FTK_XGC_2D_MESH_HH
#define _FTK_XGC_2D_MESH_HH

#include <ftk/ftk_config.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_xgc_2d_mesh : public simplicial_unstructured_2d_mesh<I, F> {
  simplicial_xgc_2d_mesh(
      const ndarray<F>& coords, 
      const ndarray<I>& triangles,
      const ndarray<F>& psi,
      const ndarray<I>& nextnodes);

  static std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> from_xgc_mesh_h5(const std::string& filename);

  I nextnode(I i) const { return nextnodes[i]; }

public:
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_xgc_slices_3d_vtu(int nphi, int iphi) const;
  vtkSmartPointer<vtkUnstructuredGrid> scalar_to_xgc_slices_3d_vtu(const std::string& varname, const ndarray<F>& data, int nphi, int iphi) const;
  void to_xgc_slices_3d_vtu(const std::string&, int nphi, int iphi) const;
  void scalar_to_xgc_slices_3d_vtu(const std::string& filename, const std::string& varname, const ndarray<F>& data, int nphi, int iphi) const;
#endif

protected:
  ndarray<F> psi;
  ndarray<I> nextnodes;
};
/////////
  
template <typename I, typename F>
simplicial_xgc_2d_mesh<I, F>::simplicial_xgc_2d_mesh(
    const ndarray<F>& coords, 
    const ndarray<I>& triangles,
    const ndarray<F>& psi_,
    const ndarray<I>& nextnodes_) : 
  simplicial_unstructured_2d_mesh<I, F>(coords, triangles), 
  psi(psi_),
  nextnodes(nextnodes_)
{
}

template <typename I, typename F>
std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> simplicial_xgc_2d_mesh<I, F>::from_xgc_mesh_h5(const std::string& filename)
{
  ndarray<I> triangles;
  ndarray<F> coords;
  ndarray<I> nextnodes;
  ndarray<F> psi;

  triangles.from_h5(filename, "/cell_set[0]/node_connect_list");
  coords.from_h5(filename, "/coordinates/values");
  psi.from_h5(filename, "/psi");
  nextnodes.from_h5(filename, "/nextnode");

  return std::shared_ptr<simplicial_xgc_2d_mesh<I, F>>(
      new simplicial_xgc_2d_mesh<I, F>(coords, triangles, psi, nextnodes));

  // return std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>>(
  //     new simplicial_unstructured_2d_mesh<I, F>(coords, triangles));
}

#if FTK_HAVE_VTK
template <typename I, typename F>
void simplicial_xgc_2d_mesh<I, F>::
to_xgc_slices_3d_vtu(const std::string& filename, int nphi, int iphi) const
{
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( to_xgc_slices_3d_vtu(nphi, iphi) );
  writer->Write();
}

template <typename I, typename F>
void simplicial_xgc_2d_mesh<I, F>::
scalar_to_xgc_slices_3d_vtu(const std::string& filename, const std::string& varname, const ndarray<F>& scalar, int nphi, int iphi) const
{
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( scalar_to_xgc_slices_3d_vtu(varname, scalar, nphi, iphi) );
  writer->Write();
}

template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_xgc_2d_mesh<I, F>::
to_xgc_slices_3d_vtu(int nphi, int iphi) const
{
  const int np = nphi * iphi;
  
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(this->n(0) * np);
 
  for (int p = 0; p < np; p ++) {
    const F phi = p * 2 * M_PI / np;
    const vtkIdType offset = p * this->n(0);
    for (int i=0; i < this->n(0); i++)
      pts->SetPoint(offset + i, 
          this->vertex_coords[i*2] * cos(phi),
          this->vertex_coords[i*2] * sin(phi), 
          this->vertex_coords[i*2+1]);
  }
 
  for (int p = 0; p < np; p ++) {
    const vtkIdType offset = p * this->n(0);
    for (int i=0; i<this->n(2); i ++) {
      vtkIdType ids[3] = {
        offset + this->triangles[i*3], 
        offset + this->triangles[i*3+1], 
        offset + this->triangles[i*3+2]
      };
      grid->InsertNextCell(VTK_TRIANGLE, 3, ids);
    }
  }

  grid->SetPoints(pts);
  return grid;
}

template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_xgc_2d_mesh<I, F>::
scalar_to_xgc_slices_3d_vtu(const std::string& varname, const ndarray<F>& scalar, int nphi, int iphi) const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = to_xgc_slices_3d_vtu(nphi, iphi);

  vtkSmartPointer<vtkDataArray> array = vtkDoubleArray::New();
  array->SetName(varname.c_str());
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(this->n(0) * nphi * iphi);

  // std::cerr << scalar << std::endl;

  for (int i = 0; i < iphi; i ++) 
    for (int j = 0; j < nphi; j ++) 
      for (int k = 0; k < this->n(0); k ++)
        array->SetTuple1((i * nphi + j) * this->n(0) + k, scalar(k, j));

  grid->GetPointData()->AddArray(array);
  grid->GetPointData()->SetActiveScalars(varname.c_str());

  // grid->PrintSelf(std::cerr, vtkIndent(2));
  return grid;
}
#endif


} // namespace ftk

#endif
