#ifndef _FTK_REGULAR_TRACKER_HH
#define _FTK_REGULAR_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/decomposition.hpp>

namespace ftk {

struct regular_tracker : public virtual tracker {
  regular_tracker(diy::mpi::communicator comm, int nd/*2 or 3*/) : tracker(comm), m(nd+1) {}
  virtual ~regular_tracker() {}
 
public:
  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void set_local_domain(const lattice&); // rank-specific "core" region of the block
  void set_local_array_domain(const lattice&); // rank-specific "ext" region of the block

  void set_coordinates(const ndarray<double>& coords_) {coords = coords_; use_explicit_coords = true;}

  void initialize();

protected:
  struct block_t {
    int gid;
    lattice local_domain, local_array_domain;
  };

protected:
  typedef simplicial_regular_mesh_element element_t;
  simplicial_regular_mesh m; // spacetime mesh
  
  lattice domain, array_domain, 
          local_domain, local_array_domain;

  bool use_explicit_coords = false;
  ndarray<double> coords;

protected: // internal use
  template <typename I=int> void simplex_indices(const std::vector<std::vector<int>>& vertices, I indices[]) const;

  void element_for_ordinal(int k, std::function<void(element_t)> f) { element_for(true, k, f); }
  void element_for_interval(int k, std::function<void(element_t)> f) { element_for(false, k, f); }
  void element_for(bool ordinal, int k, std::function<void(element_t)> f);
};

/////////////////////////////
inline void regular_tracker::initialize()
{
#if 0
  fprintf(stderr, "ordinal 1-simplices:\n");
  m.print_unit_simplices(1, ELEMENT_SCOPE_ORDINAL);
  
  fprintf(stderr, "interval 1-simplices:\n");
  m.print_unit_simplices(1, ELEMENT_SCOPE_INTERVAL);
#endif

  // initialize spacetime mesh
  {
    auto lb = domain.lower_bounds(), 
         ub = domain.upper_bounds();

    lb.push_back(start_timestep);
    ub.push_back(end_timestep);

    m.set_lb_ub(lb, ub);
  }
  
  // initialize spattial domain decomposittion
  if (use_default_domain_partition) {
    lattice_partitioner partitioner(domain);
   
    // a ghost size of 2 is necessary for jacobian derivaition; 
    // even if jacobian is not necessary, a ghost size of 1 is 
    // necessary for accessing values on boundaries
    std::vector<size_t> ghost;
    for (int i = 0; i < domain.nd(); i ++)
      ghost.push_back(2);
    
    partitioner.partition(comm.size(), {}, ghost);

    local_domain = partitioner.get_core(comm.rank());
    local_array_domain = partitioner.get_ext(comm.rank());
  }

  if (!is_input_array_partial)
    local_array_domain = array_domain;

#if 0 // experimental code to use DIY for multi-block handling
  // iniitalize partitions
  diy::RoundRobinAssigner assigner(comm.size(), nblocks);
  
  typedef diy::DiscreteBounds Bounds; 
  typedef diy::RegularGridLink RGLink;

  const int nd = m.nd() - 1;
  std::vector<bool> share_face, wrap;
  std::vector<int> ghosts(nd), divs;

  for (int i = 0; i < nd; i ++)
    ghosts[i] = 2;
  
  diy::RegularDecomposer<Bounds> decomposer(
      nd, domain.to_diy_bounds(), nblocks, 
      share_face, wrap, ghosts, divs);
  decomposer.decompose(comm.rank(), assigner, 
    [&](int gid, const Bounds &core, const Bounds &bounds, const Bounds&, const RGLink &link) {
      block_t *b = new block_t;
      RGLink *l = new RGLink(link);

      b->gid = gid;
      b->local_domain = core;
      b->local_array_domain = bounds;

      master.add(gid, b, l);
    });

  master.foreach([&](block_t *b, const diy::Master::ProxyWithLink&) {
    // fprintf(stderr, "hello %d\n", b->gid);
    // std::cerr << b->local_domain << std::endl;
    // std::cerr << b->local_array_domain << std::endl;
  });
#endif
}

template <typename I>
inline void regular_tracker::simplex_indices(
    const std::vector<std::vector<int>>& vertices, I indices[]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    indices[i] = m.get_lattice().to_integer(vertices[i]);
}

inline void regular_tracker::element_for(bool ordinal, int k, std::function<void(element_t)> f) 
{
  auto st = local_domain.starts(), sz = local_domain.sizes();
  st.push_back(current_timestep);
  sz.push_back(1);

  lattice local_spacetime_domain(st, sz);
  // std::cerr << local_spacetime_domain << std::endl;

  m.element_for(k, local_spacetime_domain, 
      ordinal ? ELEMENT_SCOPE_ORDINAL : ELEMENT_SCOPE_INTERVAL, 
      f, xl, nthreads, enable_set_affinity);
}

}


#endif
