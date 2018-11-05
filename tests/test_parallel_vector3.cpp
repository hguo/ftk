#include <iostream>
#include <vector>
#include <netcdf.h>
#include <ftk/numerics/parallel_vector.hh>
#include <ftk/numerics/cross_product.hh>
#include <ftk/numerics/norm.hh>
#include <ftk/numerics/barycentric_interpolation.hh>
#include <ftk/mesh_graph/MeshGraphRegular3DTets.h>

#define NC_SAFE_CALL(call) {\
  int retval = call;\
  if (retval != 0) {\
    fprintf(stderr, "[NetCDF Error] %s, in file '%s', line %i.\n", nc_strerror(retval), __FILE__, __LINE__); \
    exit(EXIT_FAILURE); \
  }\
}

template <typename T>
const T& texel3D(const std::vector<T> &p, int nv, int W, int H, int D, int v, int x, int y, int z)
{
  return p[(x + W*(y + H*z))*nv + v]; 
}

template <typename T>
T& texel3D(std::vector<T> &p, int nv, int W, int H, int D, int v, int x, int y, int z)
{
  return p[(x + W*(y + H*z))*nv + v]; 
}

const int W=128, H=128, D=128;
std::vector<float> VV, WW;

void open_tornado_nc(const std::string& filename, std::vector<float> &V)
{
  int ncid; 
  int varid[3]; 
  size_t starts[4] = {0, 0, 0, 0},
         sizes[4]  = {1, D, H, W};
 
  const std::string varname_u("U"), varname_v("V"), varname_w("W");
  const size_t n = W*H*D;
  std::vector<float> uu(n), vv(n), ww(n);
  V.resize(3*n);

  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) ); 
  NC_SAFE_CALL( nc_inq_varid(ncid, varname_u.c_str(), &varid[0]) ); 
  NC_SAFE_CALL( nc_inq_varid(ncid, varname_v.c_str(), &varid[1]) ); 
  NC_SAFE_CALL( nc_inq_varid(ncid, varname_w.c_str(), &varid[2]) ); 
  NC_SAFE_CALL( nc_get_vara_float(ncid, varid[0], starts, sizes, &uu[0]) ); 
  NC_SAFE_CALL( nc_get_vara_float(ncid, varid[1], starts, sizes, &vv[0]) ); 
  NC_SAFE_CALL( nc_get_vara_float(ncid, varid[2], starts, sizes, &ww[0]) ); 
  NC_SAFE_CALL( nc_close(ncid) );

  for (size_t i=0; i<n; i++) {
    V[i*3] = uu[i];
    V[i*3+1] = vv[i];
    V[i*3+2] = ww[i];
  }
}

void sweep_faces(const std::vector<float> &V, const std::vector<float> &V1) 
{
  int dims[3] = {W, H, D};
  bool pbc[3] = {0};
  ftk::MeshGraphRegular3DTets mg(dims, pbc);

  for (auto i = 0; i < mg.NFaces(); i++) {
    const auto f = mg.Face(i, true);
    if (f.Valid()) {
      const float VV[] = {
        V[f.nodes[0]*3], V[f.nodes[1]*3], V[f.nodes[2]*3], 
        V[f.nodes[0]*3+1], V[f.nodes[1]*3+1], V[f.nodes[2]*3+1], 
        V[f.nodes[0]*3+2], V[f.nodes[1]*3+2], V[f.nodes[2]*3+2]
      };
      
      const float WW[] = {
        V1[f.nodes[0]*3], V1[f.nodes[1]*3], V1[f.nodes[2]*3], 
        V1[f.nodes[0]*3+1], V1[f.nodes[1]*3+1], V1[f.nodes[2]*3+1], 
        V1[f.nodes[0]*3+2], V1[f.nodes[1]*3+2], V1[f.nodes[2]*3+2]
      };

      float lambda[3];
      auto b = ftk::parallel_vector(VV, WW, lambda);
      if (b) {
        float v[3], w[3], r[3], c[3], cn;
        ftk::barycentric_interpolation3(VV, lambda, v);
        ftk::barycentric_interpolation3(WW, lambda, w);
        ftk::cross_product(v, w, c);
        cn = ftk::norm2_3(c);

        for (int k=0; k<3; k++) 
          r[k] = v[k] / w[k];

        // if (cn <= 1e-3) { 
        if (1) {
          fprintf(stderr, "face={%llu, %llu, %llu}, lambda={%f, %f, %f}\n", 
              f.nodes[0], f.nodes[1], f.nodes[2],
              lambda[0], lambda[1], lambda[2]);
          fprintf(stderr, "lambda={%f, %f, %f}, v={%f, %f, %f}, w={%f, %f, %f}, ||v x w||=%f\n", 
              lambda[0], lambda[1], lambda[2],
              v[0], v[1], v[2], w[0], w[1], w[2], cn);
        }
      }
    }
  }
}

int main(int argc, char **argv)
{
  open_tornado_nc(argv[1], VV);
  open_tornado_nc(argv[2], WW);
  sweep_faces(VV, WW);

  // fprintf(stderr, "%f\n", texel3D<float>(gradV, 9, W, H, D, 0, 64, 64, 64));
  // fprintf(stderr, "%f\n", texel3D<float>(V, 3, W, H, D, 0, 64, 64, 64));
  return 0;
}
