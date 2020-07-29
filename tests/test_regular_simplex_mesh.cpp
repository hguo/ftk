#include <ftk/basic/distributed_union_find.hh>
#include <ftk/basic/duf.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/external/diy-ext/gather.hh>

int main(int argc, char **argv)
{
  diy::mpi::environment env;

  ftk::simplicial_regular_mesh m(3);
  m.set_lb_ub({0, 0, 0}, {1, 1, 1});
  m.element_for(2, [&](ftk::simplicial_regular_mesh_element e) {
      e.print(std::cerr, m);
      std::cerr << std::endl;
    }, 1);

  return 0;

#if 0
  diy::mpi::communicator comm;
  ftk::duf<int> uf;
  uf.pid = [](int i) {
    if (i < 8) return 0; 
    else return 1;
  };

  if (comm.rank() == 0) {
    uf.unite(0, 1);
    uf.unite(0, 2);
    uf.unite(3, 4);
    uf.unite(2, 3);
    uf.unite(1, 2);
    uf.unite(1, 3);

    uf.unite(5, 6);
    uf.unite(5, 7);
    uf.unite(6, 7);
    uf.unite(6, 8);
    uf.unite(7, 8);
  } else {
    uf.unite(8, 6);
    uf.unite(8, 7);
    uf.unite(8, 9);

    uf.unite(10, 11);
    uf.unite(10, 12);
    uf.unite(11, 12);
  }

  uf.exchange();

  for (const auto &kv : uf.parents) {
    fprintf(stderr, "pid=%d, %d --> %d, is_root=%d\n", comm.rank(), kv.first, kv.second, kv.first == uf.find(kv.second));
  }
#endif  

#if 0
  MPI_Init(&argc, &argv);
  diy::mpi::communicator comm;

  std::map<int, int> mymap, mymap1;
  // if (comm.rank() != 2) 
    mymap[comm.rank()] = comm.rank();

  diy::mpi::gather(comm, mymap, mymap1, 0);

  if (comm.rank() == 0) {
    for (const auto kv : mymap1) 
      fprintf(stderr, "k=%d, v=%d\n", kv.first, kv.second);
  }
  MPI_Finalize();
#endif

#if 0
  ftk::lattice l({0, 0, 0}, {256, 256, 256});
  ftk::lattice_partitioner partitioner(l);

  // partitioner.partition(14);
  partitioner.partition(1024, {}, {2, 2, 2});

  std::cerr << partitioner <<  std::endl;
#endif

#if 0
  ftk::simplicial_regular_mesh m(3);
  m.set_lb_ub({0, 0, 0}, {3, 3, 3});

  m.element_for(2, [=](ftk::simplicial_regular_mesh_element e) {
      std::cerr << e << std::endl;
    }, 1);
#endif
 
#if 0
  m.element_for_ordinal(2, 1, [=](ftk::simplicial_regular_mesh_element e) {
      std::cerr << e << std::endl;
    }, 1);
#endif
 
#if 0
  m.element_for_interval(2, 1, 2, [=](ftk::simplicial_regular_mesh_element e) {
      std::cerr << e << std::endl;
    }, 1);
#endif

#if 0
  ftk::lattice my_lattice({1, 1, 1}, {4, 4, 1});
  my_lattice.print(std::cerr);
  
  std::vector<int> idx = {2, 2, 1};
  
  auto id = my_lattice.to_integer(idx);
  fprintf(stderr, "idx={%d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], id);

  idx = my_lattice.from_integer(id);
  fprintf(stderr, "idx={%d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], id);

  id = my_lattice.to_integer(idx);
  fprintf(stderr, "idx={%d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], id);

  idx = my_lattice.from_integer(id);
  fprintf(stderr, "idx={%d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], id);
#endif

#if 0
  // ftk::lattice my_lattice({4, 4, 4, 4}, {128, 128, 128, 56});
  ftk::lattice my_lattice({0, 0, 0, 0}, {128, 128, 128, 56});
  std::vector<int> idx = {13, 53, 7, 34};
  
  auto id = my_lattice.to_integer(idx);
  fprintf(stderr, "idx={%d, %d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], idx[3], id);

  idx = my_lattice.from_integer(id);
  fprintf(stderr, "idx={%d, %d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], idx[3], id);
  
  id = my_lattice.to_integer(idx);
  fprintf(stderr, "idx={%d, %d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], idx[3], id);
  
  idx = my_lattice.from_integer(id);
  fprintf(stderr, "idx={%d, %d, %d, %d}, id=%d\n", idx[0], idx[1], idx[2], idx[3], id);

  fprintf(stderr, "n=%d\n", my_lattice.n());
  fprintf(stderr, "prod=%d, %d, %d, %d\n", my_lattice.prod_[0], my_lattice.prod_[1], my_lattice.prod_[2], my_lattice.prod_[3]);
#endif

#if 0
  const int nd = 4;
  ftk::simplicial_regular_mesh m(nd);

  for (int i = 0; i < nd+1; i ++)
    fprintf(stderr, "%d, %d, %d\n", 
        m.ntypes(i), 
        m.ntypes(i, ftk::ELEMENT_SCOPE_ORDINAL), 
        m.ntypes(i, ftk::ELEMENT_SCOPE_INTERVAL));

#endif
  return 0;
}
