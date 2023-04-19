#ifndef _FTK_XGC_EQ_HH
#define _FTK_XGC_EQ_HH

#include <ftk/config.hh>
#include <ftk/utils/string.hh>
#include <fstream>
#include <string>
#include <ftk/external/json.hh>

#if FTK_HAVE_VTK
#include <vtkCellLocator.h>
#include <vtkTupleInterpolator.h>
#include <vtkLagrangeQuadrilateral.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#endif

namespace ftk {

struct xgc_eq_t {
  void parse(const std::string& filename);
  void print() const;

  void initialize(); // initialize shebangs for the interpolator
  void eval_b(const double rz[], double b[]) const;
  void eval_psi(const double rz[], double &psi, double &dpsi_dR, double &dpsi_dZ) const;
  double eval_I(double psi) const;
  
#if FTK_HAVE_ADIOS2
  void read_bp(const std::string filename, diy::mpi::communicator comm = MPI_COMM_WORLD);
  void read_bp(adios2::IO &io, adios2::Engine& reader);
#endif

#if FTK_HAVE_HDF5
  void to_h5(const std::string filename) const;
#endif

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkRectilinearGrid> to_vtr() const;
  void to_vtr_file(const std::string filename) const;
  
  vtkSmartPointer<vtkImageData> to_vti() const;
  void to_vti_file(const std::string filename) const;
  
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const;
  void to_vtu_file(const std::string filename) const;
#endif

  int mr, mz, mpsi;
  double min_r, max_r, min_z, max_z; 
  double axis_r, axis_z, axis_b;
  double x_psi, x_r, x_z;

  ndarray<double> psigrid, rgrid, zgrid;
  ndarray<double> I;
  ndarray<double> psi_rz; // 2D array
  ndarray<double> grad_psi_rz;

protected: // for the interpolator
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> vtu; // Lagrange quad grid
  vtkSmartPointer<vtkCellLocator> locator;
  vtkSmartPointer<vtkTupleInterpolator> interp; // 1D interpolator for I
#endif
};

inline void xgc_eq_t::initialize()
{
#if FTK_HAVE_VTK
  // 1D interpolation function for I(psi)
  this->interp = vtkTupleInterpolator::New();
  interp->SetNumberOfComponents(1);
  for (int i = 0; i < this->mpsi; i ++) {
    double tuple[1] = {this->I[i]};
    interp->AddTuple(this->psigrid[i], tuple);
  }

  // 2D interpolation grid for psi and its partial derivatives
  this->vtu = this->to_vtu();
  this->locator = vtkCellLocator::New();
  locator->SetDataSet(this->vtu);
  locator->BuildLocator();
#endif
}

inline double xgc_eq_t::eval_I(double psi) const
{
  double I[1] = {0.0};
#if FTK_HAVE_VTK
  interp->InterpolateTuple(psi, I);
#endif
  return I[0];
}
  
inline void xgc_eq_t::eval_b(const double rz[], double b[]) const
{
  double psi, dpsi_dR, dpsi_dZ;
  
  eval_psi(rz, psi, dpsi_dR, dpsi_dZ);
  const double I = eval_I(psi);
  const double R = rz[0];

  b[0] = -dpsi_dZ / R;
  b[1] =  dpsi_dR / R;
  b[2] =  I / R;
}

inline void xgc_eq_t::eval_psi(const double rz[], double &psi, double &dpsi_dR, double &dpsi_dZ) const
{
#if FTK_HAVE_VTK
  double closestPoint[3], pcoords[3], dist2;
  int subId;
  double weights[16], values[16], derivs[2];

  double rz_[2] = {rz[0], rz[1]}; // make findcell() happy
  auto cid = this->locator->FindCell(rz_);
  auto cell = this->vtu->GetCell(cid);
  cell->EvaluatePosition(rz, closestPoint, subId, pcoords, dist2, weights);

  psi = 0.0;
  for (int i = 0; i < 16; i ++) {
    auto id = cell->GetPointId(i);
    values[i] = this->psi_rz[id];
    psi += weights[i] * values[i];
  }
  cell->Derivatives(subId, pcoords, values, 1, derivs);
  
  dpsi_dR = derivs[0];
  dpsi_dZ = derivs[1];
#else
  assert(false); // not compiled with VTK
#endif
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkImageData> xgc_eq_t::to_vti() const
{
  vtkSmartPointer<vtkImageData> d = vtkImageData::New();
  vtkSmartPointer<vtkDataArray> psi_rz = this->psi_rz.to_vtk_data_array("psi_rz");

  d->SetDimensions(this->rgrid.dim(0), this->zgrid.dim(0), 1);
  d->GetPointData()->SetScalars(psi_rz);
  d->SetOrigin(this->min_r, this->min_z, 0.0);
  d->SetSpacing(
      (this->max_r - this->min_r) / (mr - 1), 
      (this->max_z - this->min_z) / (mz - 1),
      1.0);

  return d;
}

inline void xgc_eq_t::to_vti_file(const std::string f) const
{
  auto grid = this->to_vti();

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
  writer->SetFileName(f.c_str());
  writer->SetInputData(grid);
  writer->Write();
}

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

inline vtkSmartPointer<vtkUnstructuredGrid> xgc_eq_t::to_vtu() const 
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();

  pts->SetNumberOfPoints(this->mr * this->mz);

  for (int i = 0; i < this->mr; i ++) {
    for (int j = 0; j < this->mz; j ++) {
      pts->SetPoint(j * this->mr + i, 
          this->rgrid[i], this->zgrid[j], 0.0);
    }
  }
  grid->SetPoints(pts);

  for (int i = 0; i < this->mr / 3; i ++) {
    for (int j = 0; j < this->mz / 3; j ++) {
      vtkSmartPointer<vtkLagrangeQuadrilateral> quad = vtkLagrangeQuadrilateral::New();

      const int npts = 16;
      quad->GetPointIds()->SetNumberOfIds(npts);
      // quad->GetPoints()->SetNumberOfPoints(npts);
      quad->SetOrder(3, 3);
      quad->Initialize();

      for (int k = 0; k < 4; k ++) {
        for (int l = 0; l < 4; l ++) {
          int idx = quad->PointIndexFromIJK(k, l, 0);
          // fprintf(stderr, "order=%d\n", idx);
          quad->GetPointIds()->SetId(idx, (i*3+k)*this->mz + (j*3+l));
          // quad->GetPointIds()->SetId(k*4+l, (i+k)*this->mz + (j+l));
          // quad->GetPoints().SetPoint(k*4+l, (double)k, (double)l, 0.0);
        }
      }

      grid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());
    }
  }
  
  vtkSmartPointer<vtkDataArray> psi_rz = this->psi_rz.to_vtk_data_array("psi_rz");
  grid->GetPointData()->SetScalars(psi_rz);

  return grid;
}

inline void xgc_eq_t::to_vtu_file(const std::string f) const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = this->to_vtu();

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
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
  fprintf(stderr, "mr=%d, mz=%d, mpsi=%d\n", mr, mz, mpsi);
  fprintf(stderr, "min_r=%f, max_r=%f, min_z=%f, max_z=%f\n", 
      min_r, max_r, min_z, max_z);
  fprintf(stderr, "axis_r=%f, axis_z=%f, axis_b=%f\n", 
      axis_r, axis_z, axis_b);
  fprintf(stderr, "x_psi=%f, x_r=%f, x_z=%f\n", 
      x_psi, x_r, x_z);
}

#if FTK_HAVE_ADIOS2
inline void xgc_eq_t::read_bp(const std::string f, diy::mpi::communicator comm)
{
#if ADIOS2_USE_MPI
  adios2::ADIOS adios(comm);
#else 
  adios2::ADIOS adios;
#endif

  adios2::IO io = adios.DeclareIO("BPReader");
  adios2::Engine reader = io.Open(f, adios2::Mode::Read);

  this->read_bp(io, reader);

  reader.Close();
}

inline void xgc_eq_t::read_bp(adios2::IO &io, adios2::Engine& reader)
{
  this->I.read_bp(io, reader, "eq_I");
  this->psi_rz.read_bp(io, reader, "eq_psi_rz");
  this->psigrid.read_bp(io, reader, "eq_psi_grid");

  auto read_int = [&reader](const std::string v, int& p) { reader.Get<int>(v, p); };
  auto read_double = [&reader](const std::string v, double& p) { reader.Get<double>(v, p); };

  read_int("eq_mpsi", this->mpsi);
  read_int("eq_mr", this->mr);
  read_int("eq_mz", this->mz);

  read_double("eq_axis_r", this->axis_r);
  read_double("eq_axis_z", this->axis_z);
  read_double("eq_axis_b", this->axis_b);
  
  read_double("eq_min_r", this->min_r);
  read_double("eq_max_r", this->max_r);
  read_double("eq_min_z", this->min_z);
  read_double("eq_max_z", this->max_z);
  
  read_double("eq_x_psi", this->x_psi);
  read_double("eq_x_r", this->x_r);
  read_double("eq_x_z", this->x_z);
  
  rgrid.reshape(mr);
  for (int i=0; i<mr; i++)
    rgrid[i] = min_r + (max_r - min_r) / (mr - 1) * i;

  zgrid.reshape(mz);
  for (int i=0; i<mz; i++) 
    zgrid[i] = min_z + (max_z - min_z) / (mz - 1) * i;
}
#endif

#if FTK_HAVE_HDF5
inline void xgc_eq_t::to_h5(const std::string f) const
{
  auto fid = H5Fcreate(f.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  auto sid_scalar = H5Screate(H5S_SCALAR);

  { // psi grid
    auto sid_psi_grid = H5Screate(H5S_SIMPLE);
    hsize_t psi_grid_dims[] = {this->psigrid.size()};
    H5Sset_extent_simple(sid_psi_grid, 1, psi_grid_dims, psi_grid_dims);

    auto did_psi_grid = H5Dcreate1(fid, "eq_psi_grid", H5T_NATIVE_DOUBLE, sid_psi_grid, H5P_DEFAULT);
    H5Dwrite(did_psi_grid, H5T_NATIVE_DOUBLE, H5S_ALL, sid_psi_grid, H5P_DEFAULT, this->psigrid.data());
    H5Dclose(did_psi_grid);

    auto did_I = H5Dcreate1(fid, "eq_I", H5T_NATIVE_DOUBLE, sid_psi_grid, H5P_DEFAULT);
    H5Dwrite(did_I, H5T_NATIVE_DOUBLE, H5S_ALL, sid_psi_grid, H5P_DEFAULT, this->I.data());
    H5Dclose(did_I);
    
    H5Sclose(sid_psi_grid);
  }

  { // rz grid
    auto sid_rz_grid = H5Screate(H5S_SIMPLE);
    hsize_t rz_grid_dims[] = {this->psi_rz.dim(1), this->psi_rz.dim(0)};
    H5Sset_extent_simple(sid_rz_grid, 2, rz_grid_dims, rz_grid_dims);

    auto did_psi_rz = H5Dcreate1(fid, "eq_psi_rz", H5T_NATIVE_DOUBLE, sid_rz_grid, H5P_DEFAULT);
    H5Dwrite(did_psi_rz, H5T_NATIVE_DOUBLE, H5S_ALL, sid_rz_grid, H5P_DEFAULT, this->psi_rz.data());
    H5Dclose(did_psi_rz);

    H5Sclose(sid_rz_grid);
  }

  {
    auto did = H5Dcreate1(fid, "eq_mr", H5T_NATIVE_INT, sid_scalar, H5P_DEFAULT);
    int mr = this->rgrid.size();
    H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mr);
    H5Dclose(did);
  }
 
  {
    auto did = H5Dcreate1(fid, "eq_mz", H5T_NATIVE_INT, sid_scalar, H5P_DEFAULT);
    int mz = this->zgrid.size();
    H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mz);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_mpsi", H5T_NATIVE_INT, sid_scalar, H5P_DEFAULT);
    int mpsi = this->zgrid.size();
    H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mpsi);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_axis_r", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->axis_r);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_axis_z", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->axis_z);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_axis_b", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->axis_b);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_min_r", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->min_r);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_max_r", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->max_r);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_min_z", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->min_z);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_max_z", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->max_z);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_x_psi", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->x_psi);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_x_r", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->x_r);
    H5Dclose(did);
  }
  
  {
    auto did = H5Dcreate1(fid, "eq_x_z", H5T_NATIVE_DOUBLE, sid_scalar, H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->x_z);
    H5Dclose(did);
  }

  H5Sclose(sid_scalar);
  H5Fclose(fid);
}
#endif

} // namespace

#endif
