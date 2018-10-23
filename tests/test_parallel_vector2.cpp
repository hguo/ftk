#include <iostream>
#include <vector>
#include <netcdf.h>
#include <ftk/numerics/parallel_vector.hh>
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
std::vector<float> V, gradV;

void open_tornado_nc(const std::string& filename)
{
  int ncid; 
  int varid[3]; 
  size_t starts[4] = {0, 0, 0, 0},
         sizes[4]  = {1, D, H, W};
 
  const std::string varname_u("U"), varname_v("V"), varname_w("W");
  const size_t n = W*H*D;
  std::vector<float> uu(n), vv(n), ww(n);
  V.resize(3*n);
  gradV.resize(9*n);

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

void derive_grad()
{
  for (int i=0; i<W; i++) {
    for (int j=0; j<H; j++) {
      for (int k=0; k<D; k++) {
        for (int v=0; v<3; v++) {
          texel3D<float>(gradV, 9, W, H, D, v*3, i, j, k) = 
            0.5 * (texel3D<float>(V, 3, W, H, D, 0, i+1, j, k) 
                - texel3D<float>(V, 3, W, H, D, 0, i-1, j, k));
          texel3D<float>(gradV, 9, W, H, D, v*3+1, i, j, k) = 
            0.5 * (texel3D<float>(V, 3, W, H, D, 0, i, j+1, k) -
            texel3D<float>(V, 3, W, H, D, 0, i, j-1, k));
          texel3D<float>(gradV, 9, W, H, D, v*3+2, i, j, k) = 
            0.5 * (texel3D<float>(V, 3, W, H, D, 0, i, j, k+1) -
            texel3D<float>(V, 3, W, H, D, 0, i, j, k-1));
        }
      }
    }
  }
}

void sweep_faces() 
{
  int dims[3] = {W, H, D};
  bool pbc[3] = {0};
  ftk::MeshGraphRegular3DTets mg(dims, pbc);

  for (auto i = 0; i < mg.NFaces(); i++) {
    const auto f = mg.Face(i, true);
    if (f.Valid()) {
      const float v[] = {
        V[f.nodes[0]*3], V[f.nodes[1]*3], V[f.nodes[2]*3], 
        V[f.nodes[0]*3+1], V[f.nodes[1]*3+1], V[f.nodes[2]*3+1], 
        V[f.nodes[0]*3+2], V[f.nodes[1]*3+2], V[f.nodes[2]*3+2]
      };

      float w[9]; // gradV * V
      for (int k = 0; k < 3; k ++) {
        float grad[9];
        for (int l = 0; l < 9; l ++) 
          grad[l] = gradV[f.nodes[k]*9 + l];

        ftk::mulmat3v(grad, v, w);
      }

      float lambda[3];
      auto b = ftk::parallel_vector(v, w, lambda);
      if (b) {
        fprintf(stderr, "face={%llu, %llu, %llu}, lambda={%f, %f, %f}\n", 
            f.nodes[0], f.nodes[1], f.nodes[2],
            lambda[0], lambda[1], lambda[2]);
      }
    }
  }
}

int main(int argc, char **argv)
{
  open_tornado_nc(argv[1]);
  derive_grad();
  sweep_faces();

  // fprintf(stderr, "%f\n", texel3D<float>(gradV, 9, W, H, D, 0, 64, 64, 64));
  // fprintf(stderr, "%f\n", texel3D<float>(V, 3, W, H, D, 0, 64, 64, 64));
  return 0;
}
