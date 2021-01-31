#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/mesh/simplicial_unstructured_extruded_3d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_3d_mesh.hh>
#include <ftk/ndarray.hh>

diy::mpi::environment env;

std::string xgc_data_path; 
std::string xgc_mesh_filename = "xgc.mesh.h5";

const int nphi = 16, iphi = 1;

#if FTK_HAVE_HDF5
TEST_CASE("xgc_mesh_3d_sides") {
  auto m2 = ftk::simplicial_xgc_2d_mesh<>::from_xgc_mesh_h5(xgc_mesh_filename);
  std::shared_ptr<ftk::simplicial_xgc_3d_mesh<>> mx(new ftk::simplicial_xgc_3d_mesh<>(m2, nphi, iphi));

  srand(0);
  for (int k = 0; k < 1000; k ++) {
    const int i = rand() % mx->n(0);
    int tri[4];
    mx->get_simplex(2, i, tri);
    fprintf(stderr, "tri=%d, tri=%d, %d, %d\n", 
        i, // mx->face_type(i),
        tri[0], tri[1], tri[2]);

    const auto tets = mx->side_of(2, i);
    for (const auto &tet : tets) {
      fprintf(stderr, "--tet=%d\n", tet);
      const auto sides = mx->sides(3, tet);
      for (const auto &side : sides) {
        fprintf(stderr, "----side=%d\n", side);
      }

      REQUIRE(sides.find(i) != sides.end());
      // if (sides.find(i) == sides.end()) {
      //   fprintf(stderr, "fatal error...\n");
      //   exit(1);
      // }
    }
  }
}

TEST_CASE("xgc_mesh_3d_find") {
  auto m2 = ftk::simplicial_xgc_2d_mesh<>::from_xgc_mesh_h5(xgc_mesh_filename);
  std::shared_ptr<ftk::simplicial_xgc_3d_mesh<>> mx(new ftk::simplicial_xgc_3d_mesh<>(m2, nphi, iphi));

  srand(0);
  for (int k = 0; k < 1000; k ++) {
    const int d = 1 + rand() % 3;
    const int i = rand() % mx->n(d);

    int verts[4] = {0};
    mx->get_simplex(d, i, verts);

    fprintf(stderr, "get %d-simplex %d: %d, %d, %d, %d\n", 
        d, i, verts[0], verts[1], verts[2], verts[3]);

    int i1;
    bool succ = mx->find_simplex(d, verts, i1);
    fprintf(stderr, "find %d-simplex %d: %d, %d, %d, %d, succ=%d\n",
        d, i1, verts[0], verts[1], verts[2], verts[3], succ);

    REQUIRE(succ);
    REQUIRE(i == i1);
  }
}
#endif

int main(int argc, char **argv)
{
  const char* path = std::getenv("FTK_XGC_TEST_DATA_PATH");
  if (path) {
    xgc_data_path = path;
    xgc_mesh_filename = xgc_data_path + "/" + xgc_mesh_filename;
  }
  // fprintf(stderr, "xgc_data_path=%s\n", xgc_data_dir.c_str());
  
  int requested = MPI_THREAD_FUNNELED, provided;
#if FTK_HAVE_MPI
  MPI_Init_thread(&argc, &argv, requested, &provided);
#endif

  Catch::Session session;
  int result = session.run(argc, argv);

#if FTK_HAVE_MPI
  MPI_Finalize();
#endif

  return result;
}
