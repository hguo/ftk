#include <fstream>
#include <mutex>
#include <thread>
#include <cassert>
#include <cxxopts.hpp>

// #define NDEBUG // Disable assert()

#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
// #include <ftk/algorithms/cca.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <hypermesh/ndarray.hh>
#include <hypermesh/regular_simplex_mesh.hh>


#include <ftk/external/diy/mpi.hpp>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/assigner.hpp>
#include <ftk/external/diy/io/bov.hpp>
#include <ftk/external/diy/algorithms.hpp>

#include <ftk/basic/distributed_union_find_test.hh>
#include "connected_critical_point.hpp"

#if FTK_HAVE_QT5
#include "widget.h"
#include <QApplication>
#endif

#if FTK_HAVE_VTK
#include <ftk/geometry/curve2vtk.hh>
#include <vtkPolyDataMapper.h>
#include <vtkTubeFilter.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#endif

// for serialization
#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>

#define LOAD_BALANCING true

#define TIME_OF_STEPS true
#define MULTITHREAD false

#define PRINT_FEATURE_DENSITY false
#define PRINT_ELE_COUNT false

#define ALL_CRITICAL_POINT 0
#define MAXIMUM_POINT      1


int CRITICAL_POINT_TYPE; // By default, track maximum points

int nthreads;

int DW, DH; // the dimensionality of the data is DW*DH
int DT; // number of timesteps

double start, end; 

hypermesh::ndarray<float> scalar, grad, hess;
diy::DiscreteBounds data_box(3); // bounds of data box
std::vector<int> data_offset(3);

hypermesh::regular_simplex_mesh m(3); // the 3D space-time mesh

hypermesh::regular_simplex_mesh block_m(3); // the 3D space-time mesh of this block
hypermesh::regular_simplex_mesh block_m_ghost(3); // the 3D space-time mesh of this block with ghost cells

std::vector<std::tuple<hypermesh::regular_lattice, hypermesh::regular_lattice>> lattice_partitions;
float threshold; // threshold for trajectories. The max scalar on each trajectory should be larger than the threshold. 
int threshold_length; // threshold for trajectory length.  

std::mutex mutex;
 
Block_Critical_Point* b = new Block_Critical_Point(); 
int gid; // global id of this block / this process
std::map<std::string, intersection_t>* intersections;

// the output sets of connected elements
std::vector<std::set<std::string>> connected_components_str; // connected components 

// the output trajectories
std::vector<std::vector<float>> trajectories;

std::string filename_time_uf_w; // record time for each round of union-find

int scaling_factor; // the factor that controls the shape of the synthesize data, default value: 15

template <typename T> // the synthetic function
T f(T x, T y, T t) 
{
  return cos(x*cos(t)-y*sin(t))*sin(x*sin(t)+y*cos(t));
}

template <typename T>
hypermesh::ndarray<T> generate_synthetic_data(int DW, int DH, int DT)
{
  hypermesh::ndarray<T> scalar;
  // scalar.reshape(DW, DH, DT);
  scalar.reshape(data_box.max[0] - data_box.min[0] + 1, data_box.max[1] - data_box.min[1] + 1, data_box.max[2] - data_box.min[2] + 1);

  // for (int k = 0; k < DT; k ++) {
  //   for (int j = 0; j < DH; j ++) {
  //     for (int i = 0; i < DW; i ++) {

  // for (int k = std::max(0, block_m_ghost.lb(2)-2); k < std::min(block_m_ghost.ub(2)+3, DT); k ++) {
  //   for (int j = std::max(0, block_m_ghost.lb(1)-2); j < std::min(block_m_ghost.ub(1)+3, DH); j ++) {
  //     for (int i = std::max(0, block_m_ghost.lb(0)-2); i < std::min(block_m_ghost.ub(0)+3, DW); i ++) {

  for (int k = 0; k < data_box.max[2] + 1 - data_offset[2]; k ++) {
    for (int j = 0; j < data_box.max[1] + 1 - data_offset[1]; j ++) {
      for (int i = 0; i < data_box.max[0] + 1 - data_offset[0]; i ++) {
        const T x = ((T(i + data_offset[0]) / (DW-1)) - 0.5) * scaling_factor,
                y = ((T(j + data_offset[1]) / (DH-1)) - 0.5) * scaling_factor, 
                t = (T(k + data_offset[2]) / (DT-1)) + 1e-4;

        scalar(i, j, k) = f(x, y, t);
      }
    }
  }

  return scalar;
}

template <typename T>
hypermesh::ndarray<T> derive_gradients2(const hypermesh::ndarray<T>& scalar)
{
  hypermesh::ndarray<T> grad;
  grad.reshape(2, scalar.dim(0), scalar.dim(1), scalar.dim(2));
  
  // for (int k = 0; k < DT; k ++) {
  //   for (int j = 1; j < DH-1; j ++) {
  //     for (int i = 1; i < DW-1; i ++) {

  // for (int k = std::max(0, block_m_ghost.lb(2)-1); k < std::min(block_m_ghost.ub(2)+2, DT); k ++) {
  //   for (int j = std::max(1, block_m_ghost.lb(1)-1); j < std::min(block_m_ghost.ub(1)+2, DH-1); j ++) {
  //     for (int i = std::max(1, block_m_ghost.lb(0)-1); i < std::min(block_m_ghost.ub(0)+2, DW-1); i ++) {

  for (int k = std::max(0, block_m_ghost.lb(2)-1) - data_offset[2]; k < std::min(block_m_ghost.ub(2)+2, DT) - data_offset[2]; k ++) {
    for (int j = std::max(1, block_m_ghost.lb(1)-1) - data_offset[1]; j < std::min(block_m_ghost.ub(1)+2, DH-1) - data_offset[1]; j ++) {
      for (int i = std::max(1, block_m_ghost.lb(0)-1) - data_offset[0]; i < std::min(block_m_ghost.ub(0)+2, DW-1) - data_offset[0]; i ++) {

        grad(0, i, j, k) = 0.5 * (scalar(i+1, j, k) - scalar(i-1, j, k)) * (DW-1);
        grad(1, i, j, k) = 0.5 * (scalar(i, j+1, k) - scalar(i, j-1, k)) * (DH-1);
      }
    }
  }
  return grad;
}

template <typename T>
hypermesh::ndarray<T> derive_hessians2(const hypermesh::ndarray<T>& grad)
{
  hypermesh::ndarray<T> hess;
  hess.reshape(2, grad.dim(0), grad.dim(1), grad.dim(2), grad.dim(3));


  // for (int k = 0; k < DT; k ++) {
  //   for (int j = 2; j < DH-2; j ++) {
  //     for (int i = 2; i < DW-2; i ++) {

  // for (int k = std::max(0, block_m_ghost.lb(2)); k < std::min(block_m_ghost.ub(2)+1, DT); k ++) {
  //   for (int j = std::max(2, block_m_ghost.lb(1)); j < std::min(block_m_ghost.ub(1)+1, DH-2); j ++) {
  //     for (int i = std::max(2, block_m_ghost.lb(0)); i < std::min(block_m_ghost.ub(0)+1, DW-2); i ++) {

  for (int k = std::max(0, block_m_ghost.lb(2)) - data_offset[2]; k < std::min(block_m_ghost.ub(2)+1, DT) - data_offset[2]; k ++) {
    for (int j = std::max(2, block_m_ghost.lb(1)) - data_offset[1]; j < std::min(block_m_ghost.ub(1)+1, DH-2) - data_offset[1]; j ++) {
      for (int i = std::max(2, block_m_ghost.lb(0)) - data_offset[0]; i < std::min(block_m_ghost.ub(0)+1, DW-2) - data_offset[0]; i ++) {

        const T H00 = hess(0, 0, i, j, k) = // ddf/dx2
          0.5 * (grad(0, i+1, j, k) - grad(0, i-1, j, k)) * (DW-1);
        const T H01 = hess(0, 1, i, j, k) = // ddf/dxdy
          0.5 * (grad(0, i, j+1, k) - grad(0, i, j-1, k)) * (DH-1);
        const T H10 = hess(1, 0, i, j, k) = // ddf/dydx
          0.5 * (grad(1, i+1, j, k) - grad(1, i-1, j, k)) * (DW-1);
        const T H11 = hess(1, 1, i, j, k) = // ddf/dy2
          0.5 * (grad(1, i, j+1, k) - grad(1, i, j-1, k)) * (DH-1);
      }
    }
  }
  return hess;
}

void decompose_mesh(int nblocks) {
  std::vector<size_t> given = {0}; // partition the 2D spatial space and 1D timespace
  // std::vector<size_t> given = {0, 0, 1}; // Only partition the 2D spatial space
  // std::vector<size_t> given = {1, 1, 0}; // Only partition the 1D temporal space

  std::vector<size_t> ghost = {1, 1, 1}; // at least 1, larger is ok

  const hypermesh::regular_lattice& lattice = m.lattice(); 
  lattice.partition(nblocks, given, ghost, lattice_partitions); 

  auto& lattice_pair = lattice_partitions[gid]; 
  hypermesh::regular_lattice& lattice_p = std::get<0>(lattice_pair); 
  hypermesh::regular_lattice& lattice_ghost_p = std::get<1>(lattice_pair); 
  
  block_m.set_lb_ub(lattice_p);
  block_m_ghost.set_lb_ub(lattice_ghost_p);

  data_box.min[0] = std::max(0, block_m_ghost.lb(0)-2); data_box.max[0] = std::min(block_m_ghost.ub(0)+2, DW-1); 
  data_box.min[1] = std::max(0, block_m_ghost.lb(1)-2); data_box.max[1] = std::min(block_m_ghost.ub(1)+2, DH-1); 
  data_box.min[2] = std::max(0, block_m_ghost.lb(2)-2); data_box.max[2] = std::min(block_m_ghost.ub(2)+2, DT-1); 

  data_offset = {data_box.min[0], data_box.min[1], data_box.min[2]}; 
}

void check_simplex(const hypermesh::regular_simplex_mesh_element& f)
{
  if (!f.valid()) return; // check if the 2-simplex is valid
  const auto &vertices = f.vertices(); // obtain the vertices of the simplex
  float g[3][2], value[3];

  for (int i = 0; i < 3; i ++) {
    int _i = vertices[i][0] - data_offset[0];
    int _j = vertices[i][1] - data_offset[1];
    int _k = vertices[i][2] - data_offset[2];

    g[i][0] = grad(0, _i, _j, _k);
    g[i][1] = grad(1, _i, _j, _k);
    value[i] = scalar(_i, _j, _k);
  }
 
  float mu[3];
  bool succ = ftk::inverse_lerp_s2v2(g, mu);
  
  if (!succ) return;

  if(CRITICAL_POINT_TYPE == MAXIMUM_POINT) {

    float hessxx[3], hessxy[3], hessyy[3];
    for (int i = 0; i < vertices.size(); i ++) {
      int _i = vertices[i][0] - data_offset[0];
      int _j = vertices[i][1] - data_offset[1];
      int _k = vertices[i][2] - data_offset[2];

      hessxx[i] = hess(0, 0, _i, _j, _k);
      hessxy[i] = hess(0, 1, _i, _j, _k);
      hessyy[i] = hess(1, 1, _i, _j, _k);
    }
    float hxx = ftk::lerp_s2(hessxx, mu),
          hxy = ftk::lerp_s2(hessxy, mu), 
          hyy = ftk::lerp_s2(hessyy, mu);
    float eig[2];
    ftk::solve_eigenvalues_symmetric2x2(hxx, hxy, hyy, eig);
  
    if (eig[0] < 0 && eig[1] < 0) { 
      
    } else {
      return ;
    }

  }

  // A critical point is detected

  float X[3][3];
  for (int i = 0; i < vertices.size(); i ++)
    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j];

  float x[3];
  ftk::lerp_s2v3(X, mu, x);

  intersection_t I;
  I.eid = f.to_string();
  I.val = ftk::lerp_s2(value, mu);

  for(int i = 0; i < 3; ++i) {
    I.x.push_back(x[i]); 
    I.corner.push_back(f.corner[i]); 
  }

  {
    #if MULTITHREAD 
      std::lock_guard<std::mutex> guard(mutex);
    #endif

    intersections->insert(std::make_pair(I.eid, I)); 
    // fprintf(stderr, "x={%f, %f}, t=%f, val=%f\n", I.x[0], I.x[1], I.x[2], I.val);
  }

}

void scan_intersections() 
{
  block_m_ghost.element_for(2, check_simplex, nthreads); // iterate over all 2-simplices
}

bool is_in_mesh(const hypermesh::regular_simplex_mesh_element& f, const hypermesh::regular_lattice& _lattice) { 
  // If the corner of the face is contained by the core lattice _lattice, we consider the element belongs to _lattice

  for (int i = 0; i < f.corner.size(); ++i){
    if (f.corner[i] < _lattice.start(i) || f.corner[i] > _lattice.upper_bound(i)) {
      return false; 
    }
  }

  return true;
}

bool is_in_mesh(const std::string& eid, const hypermesh::regular_lattice& _lattice) { 
  hypermesh::regular_simplex_mesh_element f = hypermesh::regular_simplex_mesh_element(m, 2, eid); 
  
  return is_in_mesh(f, _lattice); 
}

  // If the intersection is contained in the lattice
bool is_in_mesh(const intersection_t& intersection, const hypermesh::regular_lattice& _lattice) { 
  return is_in_mesh(intersection.eid, _lattice); 
}


// Add union edges to the block
void add_unions(const hypermesh::regular_simplex_mesh_element& f) {
  if (!f.valid()) return; // check if the 3-simplex is valid

  const auto elements = f.sides();
  std::set<std::string> features;
  std::set<std::string> features_in_block;  
  std::map<std::string, hypermesh::regular_simplex_mesh_element> id2element; 

  for (const auto& ele : elements) {
    std::string eid = ele.to_string(); 
    if(intersections->find(eid) != intersections->end()) {
      features.insert(eid); 
      id2element.insert(std::make_pair(eid, ele)); 

      if(is_in_mesh(ele, block_m.lattice())) {
        features_in_block.insert(eid); 
      }
      
    }
  }

  if(features.size()  > 1) {
    if(features_in_block.size() > 1) {
      // When features are local, we just need to relate to the first feature element

      #if MULTITHREAD
        std::lock_guard<std::mutex> guard(mutex);
      #endif

      std::string first_feature = *(features_in_block.begin()); 
      for(std::set<std::string>::iterator ite_i = std::next(features_in_block.begin(), 1); ite_i != features_in_block.end(); ++ite_i) {
        std::string curr_feature = *ite_i; 

        if(first_feature < curr_feature) { // Since only when the id of related_ele < ele, then the related_ele can be the parent of ele
          intersections->find(curr_feature)->second.related_elements.insert(first_feature); 
        } else {
          intersections->find(first_feature)->second.related_elements.insert(curr_feature); 
        }
        
      }
    }

    if(features_in_block.size() == 0 || features.size() == features_in_block.size()) {
      return ;
    }
  
    // When features are across processors, we need to relate all local feature elements to all remote feature elements
    for(auto& feature: features) {
      if(features_in_block.find(feature) == features_in_block.end()) { // if the feature is not in the block
        #if MULTITHREAD 
          std::lock_guard<std::mutex> guard(mutex); // Use a lock for thread-save. 
        #endif

        for(auto& feature_in_block : features_in_block) {

          if(feature < feature_in_block) { // When across processes, also, since only when the id of related_ele < ele, then the related_ele can be the parent of ele
            intersections->find(feature_in_block)->second.related_elements.insert(feature); 
          }

        }
      }
    }
  }
}

void add_related_elements_to_intersections() {
    // std::cout<<"Start Adding Elements to Blocks: "<<world.rank()<<std::endl; 

  std::vector<hypermesh::regular_simplex_mesh_element> eles_with_intersections;
  for(auto& pair : *intersections) {
    auto& eid = pair.first;
    // hypermesh::regular_simplex_mesh_element f = hypermesh::regular_simplex_mesh_element(m, 2, eid); 
    auto&f = eles_with_intersections.emplace_back(m, 2, eid); 
  }

  // std::cout<<"Finish Adding Elements to Blocks: "<<world.rank()<<std::endl; 

  // Connected Component Labeling by using union-find. 

  // std::cout<<"Start Adding Union Operations of Elements to Blocks: "<<world.rank()<<std::endl; 

  // // Method one:
  //   // For dense critical points
  //   // Enumerate each 3-d element to connect 2-d faces that contain critical points  
  // _m_ghost.element_for(3, add_unions, nthreads); 

  // Method two:
    // For sparse critical points
    // Enumerate all critical points, find their higher-order geometry; to connect critical points in this higher-order geometry
  std::set<std::string> visited_hypercells;
  for(auto& e : eles_with_intersections) { 
    const auto hypercells = e.side_of();
    for(auto& hypercell : hypercells) {
      std::string id_hypercell = hypercell.to_string(); 
      if(visited_hypercells.find(id_hypercell) == visited_hypercells.end()) {
        visited_hypercells.insert(id_hypercell); 
        add_unions(hypercell); 
      }
    }
  }

  // std::cout<<"Finish Adding Union Operations of Elements to Blocks: "<<world.rank()<<std::endl; 
}

void init_block_without_load_balancing() {

  intersections->clear(); // *
  for(auto p : b->points) {
    intersections->insert(std::make_pair(p.eid, p)); // *

    b->add(p.eid); 
    b->set_gid(p.eid, gid);
  }

  for(auto& p : b->points) {
    // std::cout<<p.related_elements.size()<<std::endl;

    for(auto& related_ele : p.related_elements) {
      // std::cout<<related_ele<<std::endl;

      b->add_related_element(p.eid, related_ele); 

      hypermesh::regular_simplex_mesh_element f = hypermesh::regular_simplex_mesh_element(m, 2, related_ele);

      if(!b->has_gid(related_ele)) { // If the block id of this feature is unknown, search the block id of this feature
        for(int i = 0; i < lattice_partitions.size(); ++i) {
          int rgid = i;
          if(rgid == gid) { // We know the feature is not in this partition
            continue;
          }

          auto& _lattice = std::get<0>(lattice_partitions[i]); 
          if(is_in_mesh(f, _lattice)) { // the feature is in mith partition
            b->set_gid(related_ele, rgid); // Set gid of this feature to mi  
          }

        }
      }

    }
  }
}

void init_block_after_load_balancing(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner) {
  // , diy::RegularContinuousLink* link

  intersections->clear(); 
  for(auto p : b->points) {
    intersections->insert(std::make_pair(p.eid, p)); 

    b->add(p.eid); 
    b->set_gid(p.eid, gid);
  }

  for(auto& p : b->points) {
    // std::cout<<p.related_elements.size()<<std::endl;

    for(auto& related_ele : p.related_elements) {
      // std::cout<<related_ele<<std::endl;

      b->add_related_element(p.eid, related_ele); 

      hypermesh::regular_simplex_mesh_element f = hypermesh::regular_simplex_mesh_element(m, 2, related_ele); 

      if(!b->has_gid(related_ele)) { // If the block id of this feature is unknown, search the block id of this feature

        for(int i = 0; i < b->block_bounds.size(); ++i) {
          int rgid = i;
          if(rgid == gid) {
            continue;
          }

          bool flag = true;
          for(int j = 0; j < 3; ++j) {
            if(f.corner[j] < b->block_bounds[i].min[j] || f.corner[j] > b->block_bounds[i].max[j]) {
              flag = false;
              break;
            }
          }
          if(flag) {
            if(b->has_gid(related_ele)) { // if this point is within multiple processes' ranges; probably on a boundary
              // std::cout<<"Multiple gids! "<<std::endl;

              int _rgid = b->get_gid(related_ele); 
              if(_rgid < rgid) {
                b->set_gid(related_ele, rgid);
              }

            } else {
              b->set_gid(related_ele, rgid);
            }

            // break ;
          }
        }

        assert(b->has_gid(related_ele)); // If not, Error! Cannot find the gid of the related element!
        // if(!b->has_gid(related_ele)) {
        //   std::cout<<"Error! Cannot find the gid of the related element! "<<std::endl;
        //   exit(0);
        // }

      }

    }
  }
}


void load_balancing(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner) {
  bool wrap = false; 
  int hist = 128; //32; 512

  hist = std::min({hist, DW, DH, DT}); 

  // DIM = 3

  diy::ContinuousBounds domain(3);

  // m.set_lb_ub({2, 2, 0}, {DW-3, DH-3, DT-1});
  domain.min[0] = 2; domain.max[0] = DW-3; 
  domain.min[1] = 2; domain.max[1] = DH-3; 
  domain.min[2] = 0; domain.max[2] = DT-1; 

  diy::RegularContinuousLink* link = new diy::RegularContinuousLink(3, domain, domain);
  master.add(gid, b, link); 

  diy::kdtree<Block_Critical_Point, intersection_t>(master, assigner, 3, domain, &Block_Critical_Point::points, 2*hist, wrap);
      // For weighted kdtree, look at kdtree.hpp diy::detail::KDTreePartition<Block,Point>::compute_local_histogram, pass and weights along with particles

  // Everybody sends their bounds to everybody else
  diy::all_to_all(master, assigner, [&](void* _b, const diy::ReduceProxy& srp) {
    Block_Critical_Point* b = static_cast<Block_Critical_Point*>(_b);
    if (srp.round() == 0) {
      diy::RegularContinuousLink* link = static_cast<diy::RegularContinuousLink*>(srp.master()->link(srp.master()->lid(srp.gid())));
      for (int i = 0; i < world.size(); ++i) {
        srp.enqueue(srp.out_link().target(i), link->bounds());
      }
    } else {
      b->block_bounds.resize(srp.in_link().size());
      for (int i = 0; i < srp.in_link().size(); ++i) {
        int _gid = srp.in_link().target(i).gid;

        assert(i == _gid);

        srp.dequeue(_gid, b->block_bounds[_gid]);
      }
    }
  });
}


void unite_disjoint_sets(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& components_str)
{
  // Initialization
    // Init union-find blocks
  std::vector<Block_Union_Find*> local_blocks;
  local_blocks.push_back(b); 

  // std::cout<<"Start Distributed Union-Find: "<<world.rank()<<std::endl; 

  // get_connected_components
  bool is_iexchange = true; // true false
  exec_distributed_union_find(world, master, assigner, local_blocks, is_iexchange, filename_time_uf_w); 

  #ifdef FTK_HAVE_MPI
    #if TIME_OF_STEPS
      MPI_Barrier(world);
      end = MPI_Wtime();
      if(world.rank() == 0) {
        std::cout << "UF: Unions of Disjoint Sets: " << end - start << " seconds. " << std::endl;
      }
      start = end; 
    #endif
  #endif

  // std::cout<<"Finish Distributed Union-Find: "<<world.rank()<<std::endl; 
}

void trace_intersections(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner)
{
  typedef hypermesh::regular_simplex_mesh_element element_t; 

  // std::cout<<"Start Extracting Connected Components: "<<world.rank()<<std::endl; 

  unite_disjoint_sets(world, master, assigner, connected_components_str);

  // std::cout<<"Finish Extracting Connected Components: "<<world.rank()<<std::endl; 

  // Obtain geometries

  // Get disjoint sets of element IDs
  // b = static_cast<Block_Critical_Point*> (master.get(0)); // load block with local id 0
  b->get_sets(world, master, assigner, connected_components_str); 

  // #ifdef FTK_HAVE_MPI
  //   #if TIME_OF_STEPS
  //     MPI_Barrier(world);
  //     end = MPI_Wtime();
  //     if(world.rank() == 0) {
  //       std::cout << "Gather Connected Components: " << end - start << " seconds. " << std::endl;
  //     }
  //     start = end; 
  //   #endif
  // #endif 

  // Convert connected components to geometries

  // if(world.rank() == 0) {
  if(connected_components_str.size() > 0) {
    
    std::vector<std::set<element_t>> cc; // connected components 
    // Convert element IDs to elements
    for(auto& comp_str : connected_components_str) {
      std::set<element_t>& comp = cc.emplace_back(); 
      for(auto& eid : comp_str) {
        comp.insert(hypermesh::regular_simplex_mesh_element(m, 2, eid)); 
      }
    }

    auto neighbors = [](element_t f) {
      std::set<element_t> neighbors;
      const auto cells = f.side_of();
      for (const auto c : cells) {
        const auto elements = c.sides();
        for (const auto f1 : elements)
          neighbors.insert(f1);
      }
      return neighbors;
    };
    
    for (int i = 0; i < cc.size(); i ++) {
      std::vector<std::vector<float>> mycurves;
      auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(cc[i], neighbors);
      
      for (int j = 0; j < linear_graphs.size(); j ++) {
        int _size = linear_graphs[j].size();
        if(_size > threshold_length) {
          std::vector<float> mycurve, mycolors;
          float max_value = std::numeric_limits<float>::min();
          for (int k = 0; k < linear_graphs[j].size(); k ++) {
            auto p = intersections->at(linear_graphs[j][k].to_string());
            mycurve.push_back(p.x[0]); //  / (DW-1));
            mycurve.push_back(p.x[1]); //  / (DH-1));
            mycurve.push_back(p.x[2]); //  / (DT-1));
            mycurve.push_back(p.val);
            max_value = std::max(max_value, p.val);
          }
          if (max_value > threshold) {
            trajectories.emplace_back(mycurve);
          }

        }

      }
    }
  }

  // #ifdef FTK_HAVE_MPI
  //   #if TIME_OF_STEPS
  //     MPI_Barrier(world);
  //     end = MPI_Wtime();
  //     if(world.rank() == 0) {
  //       std::cout << "Generate trajectories: " << end - start << " seconds. " << std::endl;
  //     }
  //     start = end; 
  //   #endif
  // #endif

  #ifdef FTK_HAVE_MPI
    #if TIME_OF_STEPS
      MPI_Barrier(world);
      end = MPI_Wtime();
      if(world.rank() == 0) {
        std::cout << "Obtain trajectories: " << end - start << " seconds. " << std::endl;
      }
      start = end; 
    #endif
  #endif
}

void print_trajectories()
{
  printf("We found %lu trajectories:\n", trajectories.size());
  for (int i = 0; i < trajectories.size(); i ++) {
    printf("--Curve %d:\n", i);
    const auto &curve = trajectories[i];
    for (int k = 0; k < curve.size()/4; k ++) {
      printf("---x=(%f, %f), t=%f, val=%f\n", curve[k*4], curve[k*4+1], curve[k*4+2], curve[k*4+3]);
    }
  }
}

void read_traj_file(const std::string& f)
{
  std::ifstream ifs(f, std::ios::in | std::ios::binary);

  while(!ifs.eof()){
    float ncoord;
    ifs.read(reinterpret_cast<char*>(&ncoord), sizeof(float)); 

    // std::cout<<ncoord<<std::endl;

    auto& traj = trajectories.emplace_back(); 
    for(int j = 0; j < (int) ncoord; ++j) {
      float coord;
      ifs.read(reinterpret_cast<char*>(&coord), sizeof(float)); 

      traj.push_back(coord); 
    }
  }

  ifs.close();
}

void write_traj_file(diy::mpi::communicator& world, const std::string& f)
{
  std::vector<float> buf; 

  for(auto& traj : trajectories) {
    buf.push_back(traj.size()); 

    // std::cout<<traj.size()<<std::endl;

    for(auto& coord : traj) {
      buf.push_back(coord);  
    }
  }

  MPI_Status status;
  MPI_File fh;

  MPI_File_open(world, f.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  // MPI_File_write_shared(fh, &buf[0], buf.size(), MPI_FLOAT, &status);
  MPI_File_write_ordered(fh, &buf[0], buf.size(), MPI_FLOAT, &status); // Collective version of write_shared // refer to: https://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node283.htm

  MPI_File_close(&fh);
}

#if FTK_HAVE_VTK
void write_traj_vtk_file(const std::string& f) {
  ftk::write_curves_vtk(trajectories, f, 4); 
}
#endif

void read_dump_file(const std::string& f)
{
  std::vector<intersection_t> vector;

  std::ifstream ifs(f);
  cereal::JSONInputArchive ar(ifs);
  ar(vector);
  ifs.close();

  for (const auto &i : vector) {
    hypermesh::regular_simplex_mesh_element e(m, 2, i.eid);
    intersections->at(e.to_string()) = i;
  }
}

void write_dump_file(const std::string& f)
{
  std::vector<intersection_t> vector;
  for (const auto &i : *intersections)
    vector.push_back(i.second);
  
  std::ofstream ofs(f);
  cereal::JSONOutputArchive ar(ofs);
  ar(vector);
  ofs.close();
}

void write_element_sets_file(diy::mpi::communicator& world, const std::string& f)
{
  // diy::mpi::io::file out(world, f, diy::mpi::io::file::wronly | diy::mpi::io::file::create | diy::mpi::io::file::sequential | diy::mpi::io::file::append);

  std::stringstream ss;
  for(auto& comp_str : connected_components_str) {
    std::vector comp_str_vec(comp_str.begin(), comp_str.end()); 
    std::sort(comp_str_vec.begin(), comp_str_vec.end()); 

    for(auto& ele_id : comp_str_vec) {
      ss<<ele_id<<" "; 
    }
    ss<<std::endl; 
  }

  const std::string buf = ss.str();
  // out.write_at_all(0, buf.c_str(), buf.size()); 

  MPI_Status status;
  MPI_File fh;

  MPI_File_open(world, f.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  // MPI_File_write_shared(fh, buf.c_str(), buf.length(), MPI_CHAR, &status);
  MPI_File_write_ordered(fh, buf.c_str(), buf.length(), MPI_CHAR, &status); // Collective version of write_shared // refer to: https://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node283.htm

  MPI_File_close(&fh);

  // out.close();
}

#if FTK_HAVE_VTK
void start_vtk_window()
{
  auto vtkcurves = ftk::curves2vtk(trajectories, 4);
  vtkcurves->Print(std::cerr);

  vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
  tubeFilter->SetInputData(vtkcurves);
  tubeFilter->SetRadius(1);
  tubeFilter->SetNumberOfSides(50);
  tubeFilter->Update();
  
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  // mapper->SetInputData(vtkcurves);
  mapper->SetInputConnection(tubeFilter->GetOutputPort());

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // a renderer and render window
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // add the actors to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(1, 1, 1); // Background color white

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = 
      vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle( style );

  renderWindowInteractor->Start();
}
#endif

int main(int argc, char **argv)
{
  diy::mpi::environment     env(0, 0); // env(NULL, NULL)
  diy::mpi::communicator    world;

  #ifdef FTK_HAVE_MPI
    #if TIME_OF_STEPS
      if(world.rank() == 0) {
        std::cout << "Start! " << std::endl;
      }
    #endif
  #endif

  int nblocks = world.size(); 

  nthreads = 1; 
  
  #if MULTITHREAD
    if(world.size() == 1) {
      nthreads = std::thread::hardware_concurrency(); 
      std::cout<<"Use MULTITHREAD! "<<std::endl; 
    }
  #endif

  diy::Master               master(world, nthreads);
  diy::ContiguousAssigner   assigner(world.size(), nblocks);

  std::vector<int> gids; // global ids of local blocks
  assigner.local_gids(world.rank(), gids);
  gid = gids[0]; // We just assign one block for each process

  std::string pattern, format;
  std::string filename_dump_r, filename_dump_w;
  std::string filename_traj_r, filename_traj_w, filename_traj_vtk_w;
  std::string filename_sets_w; 

  bool show_qt = false, show_vtk = false;

  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "input file name pattern", cxxopts::value<std::string>(pattern))
    ("f,format", "input file format", cxxopts::value<std::string>(format)->default_value("float32"))
    ("read-dump", "read dump file", cxxopts::value<std::string>(filename_dump_r))
    ("write-dump", "write dump file", cxxopts::value<std::string>(filename_dump_w))
    ("read-traj", "read traj file", cxxopts::value<std::string>(filename_traj_r))
    ("write-traj", "write traj file", cxxopts::value<std::string>(filename_traj_w))
    ("write-traj-vtk", "write traj file with vtk format", cxxopts::value<std::string>(filename_traj_vtk_w))
    ("write-sets", "write sets of connected elements", cxxopts::value<std::string>(filename_sets_w))
    ("write-time-union-find", "write time for each round of union-find", cxxopts::value<std::string>(filename_time_uf_w))
    ("w,width", "width", cxxopts::value<int>(DW)->default_value("128"))
    ("h,height", "height", cxxopts::value<int>(DH)->default_value("128"))
    ("t,timesteps", "timesteps", cxxopts::value<int>(DT)->default_value("10"))
    ("critical-point-type", "Track which type of critical points", cxxopts::value<int>(CRITICAL_POINT_TYPE)->default_value(std::to_string(MAXIMUM_POINT)))
    ("scaling-factor", "scaling factor for synthetic data", cxxopts::value<int>(scaling_factor)->default_value("15"))
    ("threshold", "threshold", cxxopts::value<float>(threshold)->default_value("0"))
    ("threshold-length", "threshold for trajectory length", cxxopts::value<int>(threshold_length)->default_value("-1"))
    ("vtk", "visualization with vtk", cxxopts::value<bool>(show_vtk))
    ("qt", "visualization with qt", cxxopts::value<bool>(show_qt))
    ("d,debug", "enable debugging");
  auto results = options.parse(argc, argv);

  #ifdef FTK_HAVE_MPI
    start = MPI_Wtime();
  #endif

  // Decompose mesh
  // ========================================

  m.set_lb_ub({2, 2, 0}, {DW-3, DH-3, DT-1}); // update the mesh; set the lower and upper bounds of the mesh

  decompose_mesh(nblocks); 

  intersections = &b->intersections; 

  // Load data
  // ========================================

  if (pattern.empty()) { // if the input data is not given, generate a synthetic data for the demo
    scalar = generate_synthetic_data<float>(DW, DH, DT);
  } else { // load the binary data

    diy::mpi::io::file in(world, pattern, diy::mpi::io::file::rdonly);
    
    diy::DiscreteBounds diy_box(3);
    diy_box.min[0] = data_box.min[2]; diy_box.max[0] = data_box.max[2]; 
    diy_box.min[1] = data_box.min[1]; diy_box.max[1] = data_box.max[1]; 
    diy_box.min[2] = data_box.min[0]; diy_box.max[2] = data_box.max[0]; 

    std::vector<unsigned> shape;
    shape.push_back(DT);
    shape.push_back(DH);
    shape.push_back(DW);

    diy::io::BOV reader(in, shape);

    // int box_size = (box.max[0] - box.min[0] + 1) * (box.max[1] - box.min[1] + 1) * (box.max[2] - box.min[2] + 1); 
    // std::vector<float> _data(box_size);
    // reader.read(box, &_data[0], true);

    // int ite = 0; 
    // for (int k = box.min[0]; k < box.max[0]+1; k ++) {
    //   for (int j = box.min[1]; j < box.max[1]+1; j ++) {
    //     for (int i = box.min[2]; i < box.max[2]+1; i ++) {
    //       scalar(i, j, k) = _data[ite++];
    //     }
    //   }
    // }

    scalar.reshape(data_box.max[0] - data_box.min[0] + 1, data_box.max[1] - data_box.min[1] + 1, data_box.max[2] - data_box.min[2] + 1);

    reader.read(diy_box, scalar.data(), true);
    // reader.read(diy_box, scalar.data());
  }

  #ifdef FTK_HAVE_MPI
    #if TIME_OF_STEPS
      MPI_Barrier(world);
      end = MPI_Wtime();
      if(world.rank() == 0) {
        std::cout << "Init Data and Partition Mesh: " << end - start << " seconds. " << std::endl;
      }
      start = end;  
    #endif
  #endif

// ========================================

  //// For debug
  
  // if(world.rank() == 0) {
  //   for (auto& _m_pair : ms) {
  //     hypermesh::regular_simplex_mesh& _m = std::get<0>(_m_pair); 
  //     hypermesh::regular_simplex_mesh& _m_ghost = std::get<1>(_m_pair); 

  //     auto sizes = _m_ghost.sizes(); 
  //     std::cout << sizes[0] << " " << sizes[1] << " " << sizes[2] << std::endl; 
  //   }
  // }
  // exit(0); 

// ========================================
  
  if (!filename_traj_r.empty()) { // if the trajectory file is given, skip all the analysis and visualize/print the trajectories
    if(world.rank() == 0) {
      read_traj_file(filename_traj_r);

      if(threshold_length > -1) {
        auto ite = trajectories.begin(); 
        while(ite != trajectories.end()) {
          if(ite->size() / 4 <= threshold_length) {
            ite = trajectories.erase(ite); 
          } else {
            ++ite ;
          }
        }
      }

      if(threshold > 0) {

        auto ite = trajectories.begin(); 
        while(ite != trajectories.end()) {
          float max_value = std::numeric_limits<float>::min();
          
          for(int i = 3; i < ite->size(); i += 4) {
            max_value = std::max(max_value, ite->at(i));
          }

          if(max_value <= threshold) {
            ite = trajectories.erase(ite); 
          } else {
            ++ite ;
          }

        } 
      }

    }
  } else { // otherwise do the analysis

    if (!filename_dump_r.empty()) { // if the dump file is given, skill the sweep step; otherwise do sweep-and-trace
      read_dump_file(filename_dump_r);
    } else { // derive gradients and do the sweep
      grad = derive_gradients2(scalar);
      hess = derive_hessians2(grad);

      // #ifdef FTK_HAVE_MPI
      //   #if TIME_OF_STEPS
      //     MPI_Barrier(world);
      //     end = MPI_Wtime();
      //     if(world.rank() == 0) {
      //       std::cout << "Derive gradients: " << end - start << " seconds. " << std::endl;
      //     }
      //     start = end; 
      //   #endif
      // #endif

      // std::cout<<"Start scanning: "<<world.rank()<<std::endl; 

      scan_intersections();

      // #ifdef FTK_HAVE_MPI
      //   #if TIME_OF_STEPS
      //     end = MPI_Wtime();
      //     std::cout << end - start <<std::endl;
      //     // std::cout << "Scan Critical Points: " << end - start << " seconds. " << gid <<std::endl;
      //     start = end; 
      //   #endif
      // #endif

      // #ifdef FTK_HAVE_MPI
      //   #if TIME_OF_STEPS
      //     MPI_Barrier(world);
      //     end = MPI_Wtime();
      //     if(world.rank() == 0) {
      //       std::cout << "Scan for Critical Points: " << end - start << " seconds. " <<std::endl;
      //     }
      //     start = end; 
      //   #endif
      // #endif


      // std::cout<<"Finish scanning: "<<world.rank()<<std::endl; 
    }

    add_related_elements_to_intersections(); 

    // #ifdef FTK_HAVE_MPI
    //   #if TIME_OF_STEPS
    //     MPI_Barrier(world);
    //     end = MPI_Wtime();
    //     if(world.rank() == 0) {
    //       std::cout << "Get Edges: " << end - start << " seconds. " << std::endl;
    //     }
    //     start = end; 
    //   #endif
    // #endif

    // std::vector<intersection_t> points; 
    for(auto intersection : *intersections) {
      if(is_in_mesh(intersection.second, block_m.lattice())) {
        b->points.emplace_back(intersection.second);  
      }
    }

    // #ifdef FTK_HAVE_MPI
    //   #if TIME_OF_STEPS
    //     MPI_Barrier(world);
    //     end = MPI_Wtime();
    //     if(world.rank() == 0) {
    //       std::cout << "Add critical points into data block: " << end - start << " seconds. " <<std::endl;
    //     }
    //     start = end; 
    //   #endif
    // #endif

    #ifdef FTK_HAVE_MPI
      #if TIME_OF_STEPS
        MPI_Barrier(world);
        end = MPI_Wtime();
        if(world.rank() == 0) {
          std::cout << "Detect Features and Their Relationships: " << end - start << " seconds. " <<std::endl;
        }
        start = end; 
      #endif
    #endif

    // ===============================

    #if PRINT_ELE_COUNT || PRINT_FEATURE_DENSITY
      int feature_ele_cnt = 0;

      if (world.rank() == 0) {
        diy::mpi::reduce<int>(world, b->points.size(), feature_ele_cnt, 0, std::plus<int>());
      } else {
        diy::mpi::reduce<int>(world, b->points.size(), 0, std::plus<int>());
      }
    #endif

    #if PRINT_ELE_COUNT
      if (world.rank() == 0) {
        std::cout<<"Feature Element Count is " << feature_ele_cnt << std::endl; 
      }
    #endif

    #if PRINT_FEATURE_DENSITY
      int element_cnt = 0; 

      m.element_for(2, [&](const hypermesh::regular_simplex_mesh_element& f){
        element_cnt++ ;
      }, nthreads);

      if (world.rank() == 0) {
        std::cout<<"Feature Density: "<< feature_ele_cnt / (float)element_cnt << std::endl; 
      }
    #endif

    #if PRINT_ELE_COUNT || PRINT_FEATURE_DENSITY
      MPI_Barrier(world);
      end = MPI_Wtime();
      start = end; 
    #endif

    // ===============================


    #if LOAD_BALANCING
      // std::cout << gid << " : " << b->points.size() << std::endl;
      if(world.size() > 1) {
        load_balancing(world, master, assigner); 
      }
      // std::cout << gid << " : " << b->points.size() << std::endl;
    #endif

    #ifdef FTK_HAVE_MPI
      #if TIME_OF_STEPS
        MPI_Barrier(world);
        end = MPI_Wtime();
        if(world.rank() == 0) {
            std::cout << "UF: Load balancing: " << end - start << " seconds. " << std::endl;
        }
        start = end; 
      #endif
    #endif

    #if LOAD_BALANCING
      if(world.size() > 1) {
        init_block_after_load_balancing(world, master, assigner); 
        master.clear();  // clear the added block, since we will add block into master again. 
      } else {
        init_block_without_load_balancing();  
      }
    #else
      init_block_without_load_balancing();
    #endif

    #ifdef FTK_HAVE_MPI
      #if TIME_OF_STEPS
        MPI_Barrier(world);
        end = MPI_Wtime();
        if(world.rank() == 0) {
            std::cout << "UF: Init Data Structure of Union-Find: " << end - start << " seconds. " << std::endl;
        }
        start = end; 
      #endif
    #endif

    if (!filename_dump_w.empty())
      write_dump_file(filename_dump_w);

    // std::cout<<"Start tracing: "<<world.rank()<<std::endl; 

    trace_intersections(world, master, assigner);

    // std::cout<<"Finish tracing: "<<world.rank()<<std::endl; 

    #ifdef FTK_HAVE_MPI
      if (!filename_traj_w.empty()) {
        write_traj_file(world, filename_traj_w);

        #ifdef FTK_HAVE_MPI
          #if TIME_OF_STEPS
            MPI_Barrier(world);
            end = MPI_Wtime();
            if(world.rank() == 0) {
              std::cout << "Output trajectories: " << end - start << " seconds. " << std::endl;
            }
            start = end; 
          #endif
        #endif
      }
    #endif

    #ifdef FTK_HAVE_MPI
      // if(world.rank() == 0) {
        if (!filename_sets_w.empty())
          write_element_sets_file(world, filename_sets_w);
      // }
    #endif
  }

#if FTK_HAVE_VTK
  if(world.rank() == 0) {
    if (!filename_traj_vtk_w.empty())
      write_traj_vtk_file(filename_traj_vtk_w);
  }
#endif

  if(world.rank() == 0) {
    if (show_qt) {
  #if FTK_HAVE_QT5
      QApplication app(argc, argv);
      QGLFormat fmt = QGLFormat::defaultFormat();
      fmt.setSampleBuffers(true);
      fmt.setSamples(16);
      QGLFormat::setDefaultFormat(fmt);

      // scalar data with full domain
      hypermesh::ndarray<float> full_scalar; full_scalar.reshape(DW, DH, DT);
      
      for (int k = data_box.min[2]; k < data_box.max[2]+1; k ++) {
        for (int j = data_box.min[1]; j < data_box.max[1]+1; j ++) {
          for (int i = data_box.min[0]; i < data_box.max[0]+1; i ++) {
            full_scalar(i, j, k) = scalar(i - data_offset[0], j - data_offset[1], k - data_offset[2]); 
          }
        }
      }

      CGLWidget *widget = new CGLWidget(full_scalar);
      widget->set_trajectories(trajectories, threshold);
      widget->show();
      return app.exec();
  #else
      fprintf(stderr, "Error: the executable is not compiled with Qt\n");
  #endif
    } else if (show_vtk) {
  #if FTK_HAVE_VTK
      start_vtk_window();
      // ftk::write_curves_vtk(trajectories, "trajectories.vtp");
  #else
      fprintf(stderr, "Error: the executable is not compiled with VTK\n");
  #endif
    } else {
      // print_trajectories();
    }
  }

  return 0;
}
