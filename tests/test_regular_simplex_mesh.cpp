#include <ftk/hypermesh/regular_simplex_mesh.hh>

int main(int argc, char **argv)
{
  const int nd = 4;
  ftk::regular_simplex_mesh m(nd);

  for (int i = 0; i < nd+1; i ++)
    fprintf(stderr, "%d, %d, %d\n", 
        m.ntypes(i), 
        m.ntypes(i, ftk::ELEMENT_SCOPE_ORDINAL), 
        m.ntypes(i, ftk::ELEMENT_SCOPE_INTERVAL));

  return 0;
}
