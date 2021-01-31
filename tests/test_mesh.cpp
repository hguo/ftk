#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/mesh/simplicial_unstructured_extruded_3d_mesh.hh>
#include <ftk/ndarray.hh>

#if FTK_HAVE_VTK
TEST_CASE("mesh_extruded_3d_unstructured_pent_sides") {
  ftk::simplicial_unstructured_3d_mesh<> m;
  m.from_vtk_unstructured_grid_file("3d.vtu");

  ftk::simplicial_unstructured_extruded_3d_mesh<> m1(m);

  srand(0);
  for (int k = 0; k < 1000; k ++) {
    const int i = rand(); // 434318;
    int tet[4];
    m1.get_simplex(3, i, tet);
    fprintf(stderr, "tet=%d, type=%d, tet=%d, %d, %d, %d\n", i, m1.tet_type(i), tet[0], tet[1], tet[2], tet[3]);

    const auto pents = m1.side_of(3, i);
    for (const auto &pent : pents) {
      fprintf(stderr, "--pent=%d\n", pent);
      const auto sides = m1.sides(4, pent);
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

#if 1
TEST_CASE("mesh_extruded_3d_unstructured_tet_sides") {
  ftk::simplicial_unstructured_3d_mesh<> m;
  m.from_vtk_unstructured_grid_file("3d.vtu");

  ftk::simplicial_unstructured_extruded_3d_mesh<> m1(m);

  srand(0);
  for (int k = 0; k < 1000; k ++) {
    const int i = rand() % 100000000; // have to limit to prevent integer overflow
    int tri[4];
    m1.get_simplex(2, i, tri);
    fprintf(stderr, "tri=%d: %d, %d, %d, type=%d\n", 
        i, 
        tri[0], tri[1], tri[2], 
        m1.simplex_type(2, i));

    const auto tets = m1.side_of(2, i);
    for (const auto &tet : tets) {
      int tetverts[4];
      m1.get_simplex(3, tet, tetverts);
      fprintf(stderr, "--tet=%d: %d, %d, %d, %d, type=%d\n", tet, tetverts[0], tetverts[1], tetverts[2], tetverts[3], m1.simplex_type(3, tet));
      const auto sides = m1.sides(3, tet);
      for (const auto &side : sides) {
        int triverts[3];
        m1.get_simplex(2, side, triverts);
        fprintf(stderr, "----side=%d: %d, %d, %d, type=%d\n", side, triverts[0], triverts[1], triverts[2], m1.simplex_type(2, side));
      }

      REQUIRE(sides.find(i) != sides.end());
    }
  }
}
#endif

TEST_CASE("mesh_extruded_2d_unstructured") {
  ftk::simplicial_unstructured_2d_mesh<> m;
  m.from_vtk_unstructured_grid_file("1x1.vtu");

  ftk::simplicial_unstructured_extruded_2d_mesh<> m1(m);

  srand(0);
  for (int k = 0; k < 1000; k ++) {
    const int i = rand();
    int tri[4];
    m1.get_simplex(2, i, tri);
    fprintf(stderr, "tri=%d, type=%d, tri=%d, %d, %d\n", 
        i, m1.face_type(i),
        tri[0], tri[1], tri[2]);

    const auto tets = m1.side_of(2, i);
    for (const auto &tet : tets) {
      fprintf(stderr, "--tet=%d\n", tet);
      const auto sides = m1.sides(3, tet);
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

TEST_CASE("mesh_extruded_2d_unstructured_find") {
  ftk::simplicial_unstructured_2d_mesh<> m;
  m.from_vtk_unstructured_grid_file("1x1.vtu");

  ftk::simplicial_unstructured_extruded_2d_mesh<> m1(m);

  srand(0);
  for (int k = 0; k < 1000; k ++) {
    const int d = 1 + rand() % 3;
    const int i = rand() % m1.n(d);

    int verts[4] = {0};
    m1.get_simplex(d, i, verts);

    fprintf(stderr, "get %d-simplex %d: %d, %d, %d, %d\n", 
        d, i, verts[0], verts[1], verts[2], verts[3]);

    int i1;
    bool succ = m1.find_simplex(d, verts, i1);
    fprintf(stderr, "find %d-simplex %d: %d, %d, %d, %d, succ=%d\n",
        d, i1, verts[0], verts[1], verts[2], verts[3], succ);

    REQUIRE(succ);
    REQUIRE(i == i1);
  }
}

TEST_CASE("mesh_2d_unstructured_find") {
  ftk::simplicial_unstructured_2d_mesh<> m;
  m.from_vtk_unstructured_grid_file("1x1.vtu");

  srand(0);
  for (int k = 0; k < 1000; k ++) {
    const int d = 1 + rand() % 2;
    const int i = rand() % m.n(d);

    int verts[3] = {0};
    m.get_simplex(d, i, verts);

    fprintf(stderr, "get %d-simplex %d: %d, %d, %d\n", 
        d, i, verts[0], verts[1], verts[2]);

    int i1;
    bool succ = m.find_simplex(d, verts, i1);
    fprintf(stderr, "find %d-simplex %d: %d, %d, %d, succ=%d\n",
        d, i1, verts[0], verts[1], verts[2], succ);

    REQUIRE(succ);
    REQUIRE(i == i1);
  }
}
#endif

#include "main.hh"

#if 0
int main(int argc, char **argv)
{
  const std::string filename(argv[1]);

  ftk::ndarray<float> data;
  data.from_amira(filename);
  
  // data.reshape({data.dim(1), data.dim(2), data.dim(3)});
  data.to_binary_file(filename + ".bin");

  return 0;
}
#endif

#if 0
int main(int argc, char **argv)
{
  const std::string filename(argv[1]);

  ftk::ndarray<float> data;
  data.from_amira(filename);
  
  // data.reshape({data.dim(1), data.dim(2), data.dim(3)});
  data.to_vtk_image_data_file(filename + ".vti", true);

  return 0;

#if 0
  double su = 0, sv = 0;
  for (int k = 0; k < data.dim(3); k ++) 
    for (int j = 0; j < data.dim(2); j ++) 
      for (int i = 0; i < data.dim(1); i ++) {
        su += data(0, i, j, k);
        sv += data(1, i, j, k);
      }

  su /= data.dim(1) * data.dim(2) * data.dim(3);
  sv /= data.dim(1) * data.dim(2) * data.dim(3);
  fprintf(stderr, "su=%f, sv=%f\n", su, sv);
#endif // su=1.001280, sv=0.000033
  
  for (int k = 0; k < data.dim(3); k ++) 
    for (int j = 0; j < data.dim(2); j ++) 
      for (int i = 0; i < data.dim(1); i ++)
        data(0, i, j, k) -= 1.001280; 
  // std::cerr << data << std::endl;

#if 1
  auto d = data.slice_time();
  for (auto i = 0; i < d.size(); i ++) {
    // std::cerr << d[i] << std::endl;
    char my_filename[1024];
    sprintf(my_filename, "cylindar2D-%04d.bin", i);
    d[i].to_binary_file(my_filename); // , true);
    // sprintf(my_filename, "cylindar2D-%03d.vti", i);
    // d[i].to_vtk_image_data_file(my_filename, true);
  }
#endif

  return 0;
}
#endif

#if 0
int main(int argc, char **argv)
{
  ftk::simplicial_regular_mesh m(4);
  m.set_lb_ub({2, 2, 2, 0}, {18, 18, 18, 10});

  typedef ftk::simplicial_regular_mesh_element element_t; 
        
  auto neighbors = [&](element_t f) {
    std::set<element_t> neighbors;
    const auto cells = f.side_of(m);
    for (const auto c : cells) {
      const auto elements = c.sides(m);
      for (const auto f1 : elements)
        neighbors.insert(f1);
    }
    return neighbors;
  };

  while (1) {
    unsigned long long i;
    std::cin >> i;
    element_t e(m, 3, i); // atoi(argv[1])); // 1031755);
    std::cerr << e << std::endl << "-----\n";

    for (auto e1 : neighbors(e)) {
      std::cerr << e1 << "\t" << e1.to_integer(m) << std::endl;
    }
  }

  return 0;
}
#endif
