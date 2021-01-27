#include <ftk/config.hh>
#include <ftk/ndarray.hh>
#include "GLHeader.h"
#include "BDATReader.h"
#include "GLGPU_IO_Helper.h"
#include "glpp/GL_post_process.h"
#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

enum {
  GLGPU_ENDIAN_LITTLE = 0, 
  GLGPU_ENDIAN_BIG = 1
};

enum {
  GLGPU_TYPE_FLOAT = 0, 
  GLGPU_TYPE_DOUBLE = 1
};

static const int GLGPU_LEGACY_TAG_SIZE = 4;
static const char GLGPU_LEGACY_TAG[] = "CA02";

bool GLGPU_IO_Helper_ReadBDAT(
    const std::string& filename, 
    GLHeader &h,
    float **rho, float **phi, float **re, float **im, float **Jx, float **Jy, float **Jz,
    bool header_only, bool supercurrent)
{
  BDATReader *reader = new BDATReader(filename); 
  if (!reader->Valid()) {
    delete reader;
    return false;
  }

  std::string name, buf;
  while (1) {
    name = reader->ReadNextRecordInfo();
    if (name.size()==0) break;
    
    unsigned int type = reader->RecType(), 
                 recID = reader->RedID(); 
    float f; // temp var

    if (name != "psi")
      reader->ReadNextRecordData(&buf);
    else if (!header_only)
      reader->ReadNextRecordData(&buf);
    else break;
    void *p = (void*)buf.data();

    if (name == "dim") {
      assert(type == BDAT_INT32);
      memcpy(&h.ndims, p, sizeof(int));
      h.dims[0] = h.dims[1] = h.dims[2] = 1;
      assert(h.ndims == 2 || h.ndims == 3);
    } else if (name == "Nx") {
      assert(type == BDAT_INT32);
      memcpy(&h.dims[0], p, sizeof(int));
    } else if (name == "Ny") {
      assert(type == BDAT_INT32);
      memcpy(&h.dims[1], p, sizeof(int));
    } else if (name == "Nz") {
      assert(type == BDAT_INT32);
      memcpy(&h.dims[2], p, sizeof(int));
    } else if (name == "Lx") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.lengths[0] = f;
    } else if (name == "Ly") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.lengths[1] = f;
    } else if (name == "Lz") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.lengths[2] = f;
    } else if (name == "BC") {
      assert(type == BDAT_INT32);
      int btype; 
      memcpy(&btype, p, sizeof(int));
      h.pbc[0] = ((btype & 0x0000ff) == 0x01);
      h.pbc[1] = ((btype & 0x00ff00) == 0x0100);
      h.pbc[2] = ((btype & 0xff0000) == 0x010000); 
    } else if (name == "u") {
      // TODO
    } else if (name == "zaniso") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      // h.zaniso = p;
    } else if (name == "t") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.time = f;
    } else if (name == "Tf") {
      assert(type == BDAT_FLOAT);
    } else if (name == "Bx") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.B[0] = f;
    } else if (name == "By") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.B[1] = f;
    } else if (name == "Bz") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.B[2] = f;
    } else if (name == "Jxext") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.Jxext = f;
    } else if (name == "K") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.Kex = f;
    } else if (name == "V") {
      assert(type == BDAT_FLOAT);
      memcpy(&f, p, sizeof(float));
      h.V = f;
    } else if (name == "psi" && !header_only) {
      if (type == BDAT_FLOAT) {
        int count = buf.size()/sizeof(float)/2;
        int optype = recID == 2000 ? 0 : 1;
        float *data = (float*)p;

        *rho = (float*)malloc(sizeof(float)*count);
        *phi = (float*)malloc(sizeof(float)*count);
        *re = (float*)malloc(sizeof(float)*count);
        *im = (float*)malloc(sizeof(float)*count);

        if (optype == 0) { // re, im
#pragma omp parallel for
          for (int i=0; i<count; i++) {
            const float R = data[i*2], I = data[i*2+1];
            (*rho)[i] = sqrt(R*R + I*I);
            (*phi)[i] = atan2(I, R);
            (*re)[i] = R;
            (*im)[i] = I;
          }
        } else { // rho^2, phi
#pragma omp parallel for
          for (int i=0; i<count; i++) {
            const float Rho = sqrt(data[i*2]), Phi = data[i*2+1];
            (*rho)[i] = Rho; 
            (*phi)[i] = Phi;
            (*re)[i] = Rho * cos(Phi);
            (*im)[i] = Rho * sin(Phi);
          }
        }
      } else if (type == BDAT_DOUBLE) {
        // TODO
        assert(false);
      } else 
        assert(false);
    }
  }
  
  for (int i=0; i<3; i++) {
    h.origins[i] = -0.5 * h.lengths[i];
    if (h.pbc[i]) 
      h.cell_lengths[i] = h.lengths[i] / h.dims[i];
    else 
      h.cell_lengths[i] = h.lengths[i] / (h.dims[i] - 1);
  }

  if (supercurrent)
    GLGPU_IO_Helper_ComputeSupercurrent(h, *re, *im, Jx, Jy, Jz);

  delete reader;
  return true;
}

bool GLGPU_IO_Helper_ReadLegacy(
    const std::string& filename, 
    GLHeader& h,
    float **rho, float **phi, float **re, float **im, float **Jx, float **Jy, float **Jz,
    bool header_only, bool supercurrent)
{
  FILE *fp = fopen(filename.c_str(), "rb");
  if (!fp) return false;

  memset(&h, 0, sizeof(GLHeader));
  h.dims[0] = h.dims[1] = h.dims[2] = 1;

  // tag check
  char tag[GLGPU_LEGACY_TAG_SIZE+1] = {0};  
  fread(tag, 1, GLGPU_LEGACY_TAG_SIZE, fp);
  if (strcmp(tag, GLGPU_LEGACY_TAG) != 0) {
    fclose(fp);
    return false;
  }

  // endians
  int endian; 
  fread(&endian, sizeof(int), 1, fp); 

  // num_dims
  fread(&h.ndims, sizeof(int), 1, fp);

  // data type
  int size_real, datatype; 
  fread(&size_real, sizeof(int), 1, fp);
  if (size_real == 4) datatype = GLGPU_TYPE_FLOAT; 
  else if (size_real == 8) datatype = GLGPU_TYPE_DOUBLE; 
  else assert(false); 

  // dimensions 
  for (int i=0; i<h.ndims; i++) {
    fread(&h.dims[i], sizeof(int), 1, fp);
    if (datatype == GLGPU_TYPE_FLOAT) {
      float length; 
      fread(&length, sizeof(float), 1, fp);
      h.lengths[i] = length; 
    } else if (datatype == GLGPU_TYPE_DOUBLE) {
      fread(&h.lengths[i], sizeof(float), 1, fp); 
    }
  }

  // dummy
  int dummy; 
  fread(&dummy, sizeof(int), 1, fp);

  // time, fluctuation_amp, Bx, By, Bz, Jx
  if (datatype == GLGPU_TYPE_FLOAT) {
    float time_, fluctuation_amp_, B_[3], Jx_; 
    fread(&time_, sizeof(float), 1, fp);
    fread(&fluctuation_amp_, sizeof(float), 1, fp); 
    fread(&B_, sizeof(float), 3, fp);
    fread(&Jx_, sizeof(float), 1, fp); 
    h.time = time_;
    // _fluctuation_amp = fluctuation_amp;
    h.B[0] = B_[0]; 
    h.B[1] = B_[1]; 
    h.B[2] = B_[2];
    h.Jxext = Jx_;
  } else if (datatype == GLGPU_TYPE_DOUBLE) {
    float fluctuation_amp;
    fread(&h.time, sizeof(float), 1, fp); 
    fread(&fluctuation_amp, sizeof(float), 1, fp);
    fread(h.B, sizeof(float), 3, fp);
    fread(&h.Jxext, sizeof(float), 1, fp); 
  }

  // btype
  int btype; 
  fread(&btype, sizeof(int), 1, fp); 
  h.pbc[0] = btype & 0x0000ff;
  h.pbc[1] = btype & 0x00ff00;
  h.pbc[2] = btype & 0xff0000; 

  for (int i=0; i<3; i++) {
    h.origins[i] = -0.5 * h.lengths[i];
    if (h.pbc[i]) 
      h.cell_lengths[i] = h.lengths[i] / h.dims[i];
    else 
      h.cell_lengths[i] = h.lengths[i] / (h.dims[i] - 1);
  }

  // optype
  int optype; 
  fread(&optype, sizeof(int), 1, fp);
  if (datatype == GLGPU_TYPE_FLOAT) {
    float Kex_, Kex_dot_; 
    fread(&Kex_, sizeof(float), 1, fp);
    fread(&Kex_dot_, sizeof(float), 1, fp); 
    h.Kex = Kex_;
    h.Kex_dot = Kex_dot_;
  } else if (datatype == GLGPU_TYPE_DOUBLE) {
    float Kex_dot;
    fread(&h.Kex, sizeof(float), 1, fp);
    fread(&Kex_dot, sizeof(float), 1, fp); 
  }
 
  if (header_only) {
    fclose(fp);
    return true;
  }
  // read the actual data

  int count = 1; 
  for (int i=0; i<h.ndims; i++) 
    count *= h.dims[i]; 

  int offset = ftell(fp);

  // mem allocation 
  *rho = (float*)malloc(sizeof(float)*count);
  *phi = (float*)malloc(sizeof(float)*count);
  *re = (float*)malloc(sizeof(float)*count);
  *im = (float*)malloc(sizeof(float)*count);

  if (datatype == GLGPU_TYPE_FLOAT) {
    // raw data
    float *buf = (float*)malloc(sizeof(float)*count*2); // complex numbers
    fread(buf, sizeof(float), count*2, fp);
    
    if (optype == 0) { // re, im
#pragma omp parallel for
      for (int i=0; i<count; i++) {
        const float R = buf[i*2], I = buf[i*2+1]; 
        (*rho)[i] = sqrt(R*R + I*I);
        (*phi)[i] = atan2(I, R);
        (*re)[i] = R;
        (*im)[i] = I;
      }
    } else { // rho, phi
#pragma omp parallel for
      for (int i=0; i<count; i++) {
        const float Rho = buf[i*2], Phi = buf[i*2+1];
        (*rho)[i] = Rho; 
        (*phi)[i] = Phi;
        (*re)[i] = Rho * cos(Phi);
        (*im)[i] = Rho * sin(Phi);
      }
    }

    free(buf);
  } else if (datatype == GLGPU_TYPE_DOUBLE) {
    assert(false);
  }
  
  if (supercurrent)
    GLGPU_IO_Helper_ComputeSupercurrent(h, *re, *im, Jx, Jy, Jz);

  fclose(fp);
  return true;
}

bool GLGPU_IO_Helper_WriteRaw_rho_phi(
    const std::string& filename_rho, 
    const std::string& filename_phi, 
    GLHeader& h, 
    const float *rho, const float *phi)
{
  size_t num_elems = h.dims[2] * h.dims[1] * h.dims[0];

  FILE *fp_rho = fopen(filename_rho.c_str(), "wb");
  fwrite(rho, sizeof(float), num_elems, fp_rho);
  fclose(fp_rho);

  FILE *fp_phi = fopen(filename_phi.c_str(), "wb");
  fwrite(phi, sizeof(float), num_elems, fp_phi);
  fclose(fp_phi);

  return true;
}

bool GLGPU_IO_Helper_WriteNetCDF(
    const std::string& filename, 
    GLHeader& h,
    const float *rho, const float *phi,
    const float *re, const float *im, 
    const float *Jx, const float *Jy, const float *Jz)
{
#if FTK_HAVE_NETCDF
  int ncid; 
  int dimids[3]; 
  int varids[8];

  size_t starts[3] = {0, 0, 0}, 
         sizes[3]  = {(size_t)h.dims[2], (size_t)h.dims[1], (size_t)h.dims[0]};

#if 0
  const int cnt = sizes[0]*sizes[1]*sizes[2];
  float *rho = (float*)malloc(sizeof(float)*cnt), 
         *phi = (float*)malloc(sizeof(float)*cnt);
  for (int i=0; i<cnt; i++) {
    rho[i] = sqrt(re[i]*re[i] + im[i]*im[i]);
    phi[i] = atan2(im[i], re[i]);
  }
#endif

  fprintf(stderr, "netcdf filename=%s\n", filename.c_str());

  NC_SAFE_CALL( nc_create(filename.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &ncid) ); 
  NC_SAFE_CALL( nc_def_dim(ncid, "z", sizes[0], &dimids[0]) );
  NC_SAFE_CALL( nc_def_dim(ncid, "y", sizes[1], &dimids[1]) );
  NC_SAFE_CALL( nc_def_dim(ncid, "x", sizes[2], &dimids[2]) );
  NC_SAFE_CALL( nc_def_var(ncid, "rho", NC_FLOAT, 3, dimids, &varids[0]) );
  NC_SAFE_CALL( nc_def_var(ncid, "phi", NC_FLOAT, 3, dimids, &varids[1]) );
  NC_SAFE_CALL( nc_def_var(ncid, "re", NC_FLOAT, 3, dimids, &varids[2]) );
  NC_SAFE_CALL( nc_def_var(ncid, "im", NC_FLOAT, 3, dimids, &varids[3]) );
  NC_SAFE_CALL( nc_def_var(ncid, "Jx", NC_FLOAT, 3, dimids, &varids[4]) );
  NC_SAFE_CALL( nc_def_var(ncid, "Jy", NC_FLOAT, 3, dimids, &varids[5]) );
  NC_SAFE_CALL( nc_def_var(ncid, "Jz", NC_FLOAT, 3, dimids, &varids[6]) );
  NC_SAFE_CALL( nc_enddef(ncid) );

  NC_SAFE_CALL( nc_put_vara_float(ncid, varids[0], starts, sizes, rho) ); 
  NC_SAFE_CALL( nc_put_vara_float(ncid, varids[1], starts, sizes, phi) ); 
  NC_SAFE_CALL( nc_put_vara_float(ncid, varids[2], starts, sizes, re) ); 
  NC_SAFE_CALL( nc_put_vara_float(ncid, varids[3], starts, sizes, im) ); 
  NC_SAFE_CALL( nc_put_vara_float(ncid, varids[4], starts, sizes, Jx) ); 
  NC_SAFE_CALL( nc_put_vara_float(ncid, varids[5], starts, sizes, Jy) ); 
  NC_SAFE_CALL( nc_put_vara_float(ncid, varids[6], starts, sizes, Jz) ); 

  NC_SAFE_CALL( nc_close(ncid) );

  return true;
#else
  assert(false);
  return false;
#endif
}

void GLGPU_IO_Helper_ComputeSupercurrent(
    GLHeader &h, const float *re, const float *im, float **Jx, float **Jy, float **Jz)
{
  const int arraySize = h.dims[0] * h.dims[1] * h.dims[2];
  
  // GLPP
  GLPP *pp = new GLPP;
  // FIXME!
  pp->dim = h.ndims;
  pp->Nx = h.dims[0];
  pp->Ny = h.dims[1];
  pp->Nz = h.dims[2];
  pp->NN = arraySize;
  pp->btype = h.dtype;
  pp->Lx = h.lengths[0];
  pp->Ly = h.lengths[1];
  pp->Lz = h.lengths[2];
  pp->dx = h.cell_lengths[0];
  pp->dy = h.cell_lengths[1];
  pp->dz = h.cell_lengths[2];
  pp->Bx = h.B[0];
  pp->By = h.B[1]; 
  pp->Bz = h.B[2];
  pp->KEx = h.Kex;
  pp->psi = (COMPLEX*)malloc(sizeof(COMPLEX)*arraySize);
  for (int i=0; i<arraySize; i++) {
    pp->psi[i].re = re[i];
    pp->psi[i].im = im[i];
  }

  pp->calc_current();
  assert(pp->Jx != NULL);

  *Jx = (float*)malloc(sizeof(float)*arraySize);
  *Jy = (float*)malloc(sizeof(float)*arraySize);
  *Jz = (float*)malloc(sizeof(float)*arraySize);

  for (int i=0; i<arraySize; i++) {
    (*Jx)[i] = pp->Jx[i];
    (*Jy)[i] = pp->Jy[i];
    (*Jz)[i] = pp->Jz[i];
  }

  delete pp;
}
