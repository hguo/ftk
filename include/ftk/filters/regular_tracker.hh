#ifndef _FTK_REGULAR_TRACKER_HH
#define _FTK_REGULAR_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/decomposition.hpp>

namespace ftk {

enum {
  REGULAR_COORDS_SIMPLE = 0, // uniform
  REGULAR_COORDS_BOUNDS = 1, // uniform with bounds
  REGULAR_COORDS_RECTILINEAR = 2, // x, y, z in separate arrays
  REGULAR_COORDS_EXPLICIT = 3// explicit (x, y, z)
};

struct regular_tracker : public virtual tracker {
  regular_tracker(diy::mpi::communicator comm, int nd/*2 or 3*/) : tracker(comm), m(nd+1) {}
  virtual ~regular_tracker() {}
 
public:
  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void set_local_domain(const lattice&); // rank-specific "core" region of the block
  void set_local_array_domain(const lattice&); // rank-specific "ext" region of the block

  void initialize();

  lattice get_local_array_domain() const { return local_array_domain; }
  
public: // physics coordinates
  // void set_coordinates(const ndarray<double>& coords_) {coords = coords_; use_explicit_coords = true;}
  std::array<double, 3> get_coords(const std::vector<int>&) const;

  void set_coords_bounds(const std::vector<double>& bounds) { bounds_coords = bounds; mode_phys_coords = REGULAR_COORDS_BOUNDS; }
  void set_coords_rectilinear(const std::vector<ndarray<double>> &rectilinear_coords_) { rectilinear_coords = rectilinear_coords_; mode_phys_coords = REGULAR_COORDS_RECTILINEAR; }
  void set_coords_explicit(const ndarray<double>& explicit_coords_) {  explicit_coords = explicit_coords_; mode_phys_coords = REGULAR_COORDS_EXPLICIT; }

#if FTK_HAVE_VTK
  void set_coords_bounds(vtkSmartPointer<vtkImageData>);
  void set_coords_rectilinear(vtkSmartPointer<vtkRectilinearGrid>);
  void set_coords_explicit(vtkSmartPointer<vtkStructuredGrid>);
#endif
  
protected:
  int mode_phys_coords = REGULAR_COORDS_SIMPLE; 
  ndarray<double> explicit_coords;
  std::vector<ndarray<double>> rectilinear_coords;
  std::vector<double> bounds_coords;

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
  ndarray<double> coords; // legacy

protected: // internal use
  template <typename I=int> void simplex_indices(const std::vector<std::vector<int>>& vertices, I indices[]) const;

  void element_for_ordinal(int k, std::function<void(element_t)> f) { element_for(true, k, f); }
  void element_for_interval(int k, std::function<void(element_t)> f) { element_for(false, k, f); }
  void element_for(bool ordinal, int k, std::function<void(element_t)> f);

  bool locate_spatial_quad(const double *x, int *corner, double *quad_coords) const;
};

/////////////////////////////

inline bool regular_tracker::locate_spatial_quad(const double *x, int *corner, double *quad_coords) const
{
  if (mode_phys_coords == REGULAR_COORDS_EXPLICIT) {
    fatal(FTK_ERR_NOT_IMPLEMENTED);
    return false;
  } else if (mode_phys_coords == REGULAR_COORDS_RECTILINEAR) {
    fatal(FTK_ERR_NOT_IMPLEMENTED);
    return false;
  } else if (mode_phys_coords == REGULAR_COORDS_BOUNDS) {
    fatal(FTK_ERR_NOT_IMPLEMENTED);
    return false;
  } else { // simple
    const int nd = m.nd() - 1;
    for (int i = 0; i < nd; i ++) {
      corner[i] = x[i];
      quad_coords[i] = x[i] - corner[i];
      if (corner[i] < 0 || corner[i] >= domain.upper_bound(i))
        return false;
    }
    return true;
  }
}

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
    lattice_partitioner partitioner(domain), partitioner_array(array_domain);
   
    // a ghost size of 2 is necessary for jacobian derivaition; 
    // even if jacobian is not necessary, a ghost size of 1 is 
    // necessary for accessing values on boundaries
    std::vector<size_t> ghost;
    for (int i = 0; i < domain.nd(); i ++)
      ghost.push_back(2);
    
    // partitioner.partition(comm.size(), {}, ghost);
    partitioner.partition(comm.size()); 

    local_domain = partitioner.get_core(comm.rank());
    // local_domain = partitioner.get_ext(comm.rank());
    local_array_domain = partitioner_array.add_ghost(local_domain, ghost, ghost);

    // std::cerr << "local_domain: " << local_domain << std::endl;
    // std::cerr << "local_array_domain: " << local_array_domain << std::endl;
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
  // fprintf(stderr, "st=%d, %d, sz=%d, %d\n", 
  //     st[0], st[1], sz[0], sz[1]);

  st.push_back(current_timestep);
  sz.push_back(1);

  lattice local_spacetime_domain(st, sz);
  // std::cerr << local_spacetime_domain << std::endl;

  m.element_for(k, local_spacetime_domain, 
      ordinal ? ELEMENT_SCOPE_ORDINAL : ELEMENT_SCOPE_INTERVAL, 
      f, xl, nthreads, enable_set_affinity);
}

#if FTK_HAVE_VTK
inline void regular_tracker::set_coords_bounds(vtkSmartPointer<vtkImageData> vti)
{
  bounds_coords.resize(6);
  vti->GetBounds(&bounds_coords[0]);
  mode_phys_coords = REGULAR_COORDS_BOUNDS;
}

inline void regular_tracker::set_coords_rectilinear(vtkSmartPointer<vtkRectilinearGrid> grid)
{
  const int nd = grid->GetDataDimension();
  const int *dims = grid->GetDimensions();

  // std::shared_ptr<simplicial_regular_mesh> m(
  //   new simplicial_regular_mesh(lattice(nd, dims)));

  std::vector<ndarray<double>> rectilinear_coords(nd);
  for (int i = 0; i < nd; i ++) {
    // rectilinear_coords[i].reshape(dims[i]);
    if (i == 0)
      rectilinear_coords[i].from_vtk_data_array(grid->GetXCoordinates());
    else if (i == 1)
      rectilinear_coords[i].from_vtk_data_array(grid->GetYCoordinates());
    else  // i == 2
      rectilinear_coords[i].from_vtk_data_array(grid->GetZCoordinates());
  }

  set_coords_rectilinear(rectilinear_coords);
}

inline void regular_tracker::set_coords_explicit(vtkSmartPointer<vtkStructuredGrid> grid)
{
  const int nd = grid->GetDataDimension();
  const int *dims = grid->GetDimensions();

  // std::shared_ptr<simplicial_regular_mesh> m(
  //   new simplicial_regular_mesh(lattice(nd, dims)));

  std::vector<size_t> shape_coords;
  shape_coords.push_back(3);
  for (int i = 0; i < nd; i ++)
    shape_coords.push_back(dims[i]);

  ndarray<double> coords;
  coords.reshape(shape_coords);
  const auto np = grid->GetNumberOfPoints();
  for (auto i = 0; i < np; i ++) {
    double *p = grid->GetPoint(i);
    for (int j = 0; j < 3; j ++)
      coords[j + 3*i] = p[j];
  }

  set_coords_explicit(coords);
}
#endif

}


#endif
