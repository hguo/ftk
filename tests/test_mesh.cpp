#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/ndarray.hh>

int main(int argc, char **argv)
{
  const std::string filename(argv[1]);

  ftk::ndarray<float> data;
  data.from_amira(filename);
  
  // data.reshape({data.dim(1), data.dim(2), data.dim(3)});
  // data.to_vtk_image_data_file(filename + ".vti", true);

  std::cerr << data << std::endl;

  auto d = data.slice_time();
  for (auto i = 0; i < d.size(); i ++) {
    // std::cerr << d[i] << std::endl;
    char my_filename[1024];
    sprintf(my_filename, "cylindar2D-%03d.bin", i);
    d[i].to_binary_file(my_filename); // , true);
  }

  return 0;
}

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
