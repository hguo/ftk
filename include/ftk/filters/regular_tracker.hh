#ifndef _FTK_REGULAR_TRACKER_HH
#define _FTK_REGULAR_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/filter.hh>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/decomposition.hpp>

namespace ftk {

struct regular_tracker : public virtual filter {
  regular_tracker(int nd/*2 or 3*/) : m(nd+1), master(comm) {}
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

  diy::Master master;

protected:
  simplicial_regular_mesh m; // spacetime mesh
  
  lattice domain, array_domain, 
          local_domain, local_array_domain;

  bool use_explicit_coords = false;
  ndarray<double> coords;

protected: // internal use
  template <typename I=int> void simplex_indices(const std::vector<std::vector<int>>& vertices, I indices[]) const;
};

/////////////////////////////
inline void regular_tracker::initialize()
{
  // initialize spacetime mesh
  {
    auto lb = domain.lower_bounds(), 
         ub = domain.upper_bounds();

    lb.push_back(start_timestep);
    ub.push_back(end_timestep);

    m.set_lb_ub(lb, ub);
  }

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
}

template <typename I>
inline void regular_tracker::simplex_indices(
    const std::vector<std::vector<int>>& vertices, I indices[]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    indices[i] = m.get_lattice().to_integer(vertices[i]);
}

}


#endif
