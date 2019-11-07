#include <ftk/hypermesh/regular_simplex_mesh.hh>

int main(int argc, char **argv)
{
  ftk::lattice my_lattice({4, 4, 4, 4}, {128, 128, 128, 56});
  std::vector<int> idx = {13, 53, 7, 34};
  
  auto id = my_lattice.to_integer(idx);
  fprintf(stderr, "idx={%d, %d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], idx[3], id);

  idx = my_lattice.from_integer(id);
  fprintf(stderr, "idx={%d, %d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], idx[3], id);
  
  id = my_lattice.to_integer(idx);
  fprintf(stderr, "idx={%d, %d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], idx[3], id);
  
  idx = my_lattice.from_integer(id);
  fprintf(stderr, "idx={%d, %d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], idx[3], id);


#if 0
  const int nd = 4;
  ftk::regular_simplex_mesh m(nd);

  for (int i = 0; i < nd+1; i ++)
    fprintf(stderr, "%d, %d, %d\n", 
        m.ntypes(i), 
        m.ntypes(i, ftk::ELEMENT_SCOPE_ORDINAL), 
        m.ntypes(i, ftk::ELEMENT_SCOPE_INTERVAL));

#endif
  return 0;
}
