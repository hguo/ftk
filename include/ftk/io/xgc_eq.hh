#ifndef _FTK_XGC_EQ_HH
#define _FTK_XGC_EQ_HH

#include <ftk/config.hh>
#include <ftk/utils/string.hh>
#include <fstream>
#include <string>
#include <ftk/external/json.hh>

#if FTK_HAVE_VTK
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkSmartPointer.h>
#endif

namespace ftk {

struct xgc_eq_t {
  void parse(const std::string& filename);
  void print() const;

  void eval_b(const double rz[], double b[]) const;

#if FTK_HAVE_HDF5
  void to_h5(const std::string filename) const;
#endif

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkRectilinearGrid> to_vtr() const;
  void to_vtr_file(const std::string filename) const;
#endif

  size_t mr, mz, mpsi;
  double min_r, max_r, min_z, max_z; 
  double axis_r, axis_z, axis_b;
  double x_psi, x_r, x_z;

  ndarray<double> psigrid, rgrid, zgrid;
  ndarray<double> I;
  ndarray<double> psi_rz; // 2D array
  ndarray<double> grad_psi_rz;
};

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkRectilinearGrid> xgc_eq_t::to_vtr() const
{
  vtkSmartPointer<vtkRectilinearGrid> grid = vtkRectilinearGrid::New();
  
  vtkSmartPointer<vtkDataArray> psi_rz = this->psi_rz.to_vtk_data_array("psi_rz");
  vtkSmartPointer<vtkDataArray> rgrid = this->rgrid.to_vtk_data_array("rgrid");
  vtkSmartPointer<vtkDataArray> zgrid = this->zgrid.to_vtk_data_array("zgrid");

  grid->SetDimensions(this->rgrid.dim(0), this->zgrid.dim(0), 1);
  grid->GetPointData()->SetScalars(psi_rz);
  grid->SetXCoordinates(rgrid);
  grid->SetYCoordinates(zgrid);

  return grid;
}

inline void xgc_eq_t::to_vtr_file(const std::string f) const
{
  vtkSmartPointer<vtkRectilinearGrid> grid = this->to_vtr();

  vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkXMLRectilinearGridWriter::New();
  writer->SetFileName(f.c_str());
  writer->SetInputData(grid);
  writer->Write();
}
#endif

inline void xgc_eq_t::parse(const std::string& filename)
{
  std::ifstream ifs(filename, std::ifstream::in);
  
  // ignore the first line
  ifs.ignore(1000, '\n');

  ifs >> mr >> mz >> mpsi;
  ifs >> min_r >> max_r >> min_z >> max_z;
  ifs >> axis_r >> axis_z >> axis_b;
  ifs >> x_psi >> x_r >> x_z;
  
  psigrid.reshape(mpsi);
  for (int i=0; i<mpsi; i++) 
    ifs >> psigrid[i];

  I.reshape(mpsi);
  for (int i=0; i<mpsi; i++) 
    ifs >> I[i];

  psi_rz.reshape(mr, mz); // column major order, the dim of mr changes faster
  // ndarray<double> psi_rz_column_major; 
  // psi_rz_column_major.reshape(mz, mr);
  for (int i=0; i<mr*mz; i++) 
    // ifs >> psi_rz_column_major[i];
    ifs >> psi_rz[i]; // _column_major[i];
  // for (int i=0; i<mr; i++) // transpose
  //   for (int j=0; j<mz; j++) 
  //     psi_rz[i*mz+j] = psi_rz_column_major[j*mr+i]; 

  int end_flag;
  ifs >> end_flag;
  assert(end_flag == -1);

  ifs.close();

  rgrid.reshape(mr);
  for (int i=0; i<mr; i++)
    rgrid[i] = min_r + (max_r - min_r) / (mr - 1) * i;

  zgrid.reshape(mz);
  for (int i=0; i<mz; i++) 
    zgrid[i] = min_z + (max_z - min_z) / (mz - 1) * i;

  // print();
}

inline void xgc_eq_t::print() const 
{
  fprintf(stderr, "mr=%zu, mz=%zu, mpsi=%zu\n", mr, mz, mpsi);
  fprintf(stderr, "min_r=%f, max_r=%f, min_z=%f, max_z=%f\n", 
      min_r, max_r, min_z, max_z);
  fprintf(stderr, "axis_r=%f, axis_z=%f, axis_b=%f\n", 
      axis_r, axis_z, axis_b);
  fprintf(stderr, "x_psi=%f, x_r=%f, x_z=%f\n", 
      x_psi, x_r, x_z);
}

} // namespace

#endif
