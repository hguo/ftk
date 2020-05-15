#define DIM 3
#define FEATURE_DIM 0

#include <fstream>
#include <mutex>
#include <thread>
#include <cassert>
#include <cxxopts.hpp>

#define NDEBUG // Disable assert()

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

// #include <ftk/basic/distributed_union_find.hh>
#include "connected_component_labeling.hpp"

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

#define LOAD_BALANCING false

#define TIME_OF_STEPS true
#define MULTITHREAD false

#define PRINT_FEATURE_DENSITY false
#define PRINT_ELE_COUNT false
#define PRINT_EDGE_COUNT false

#define PRINT_SIMPLEX_COUNT false

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
float threshold; // threshold for super level sets. 
// int threshold_length; // threshold for level_set length.  

std::mutex mutex;
 
Block_Feature* b = new Block_Feature(); 
int gid; // global id of this block / this process
std::map<std::string, intersection_t>* intersections;

// the output sets of connected elements
std::vector<std::set<std::string>> connected_components_str; // connected components 

// the output level_sets
std::vector<std::vector<float>> level_sets;

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

  for (int k = 0; k < data_box.max[2] + 1 - data_offset[2]; k ++) {
    for (int j = 0; j < data_box.max[1] + 1 - data_offset[1]; j ++) {
      for (int i = 0; i < data_box.max[0] + 1 - data_offset[0]; i ++) {
        
        // const T x = ((T(i + data_offset[0]) / (DW-1)) - 0.5) * scaling_factor,
        //         y = ((T(j + data_offset[1]) / (DH-1)) - 0.5) * scaling_factor, 
        //         t = (T(k + data_offset[2]) / (DT-1)) + 1e-4;

        // For test of weak scaling
        const T x = ((T(i + data_offset[0]) / (DW-1)) - 0.5) * scaling_factor,
                y = ((T(j + data_offset[1]) / (DH-1)) - 0.5) * scaling_factor, 
                t = (T(k + data_offset[2]) / (128-1)) + 1e-4;

        scalar(i, j, k) = f(x, y, t);
      }
    }
  }

  return scalar;
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

  data_box.min[0] = block_m_ghost.lb(0); data_box.max[0] = block_m_ghost.ub(0); 
  data_box.min[1] = block_m_ghost.lb(1); data_box.max[1] = block_m_ghost.ub(1); 
  data_box.min[2] = block_m_ghost.lb(2); data_box.max[2] = block_m_ghost.ub(2); 

  data_offset = {data_box.min[0], data_box.min[1], data_box.min[2]}; 
}

void check_simplex(const hypermesh::regular_simplex_mesh_element& f)
{
  if (!f.valid()) return; // check if the 0-simplex is valid
  const auto &vertices = f.vertices(); // obtain the vertices of the simplex

  float value; 

  {
    int _i = vertices[0][0] - data_offset[0];
    int _j = vertices[0][1] - data_offset[1];
    int _k = vertices[0][2] - data_offset[2];

    value = scalar(_i, _j, _k);
  }

  if(value < threshold) {
    return ;
  }

  intersection_t I;
  I.eid = f.to_string();
  I.val = value; 
  // I.x.resize(DIM); 

  for(int i = 0; i < 3; ++i) {
    I.x[i] = f.corner[i]; 
    // I.corner.push_back(f.corner[i]); 
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
  block_m_ghost.element_for(0, check_simplex, nthreads); // iterate over all 0-simplices
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


// Vertex to string
// std::string vertex_to_string(const std::vector<int>& corner) {
//   std::stringstream res;

//   std::copy(corner.begin(), corner.end(), std::ostream_iterator<int>(res, ","));

//   return res.str(); 
// }


void trace_intersections(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner)
{
  typedef hypermesh::regular_simplex_mesh_element element_t; 

  // std::cout<<"Start Extracting Connected Components: "<<world.rank()<<std::endl; 

  unite_disjoint_sets(world, master, assigner, connected_components_str);

  // std::cout<<"Finish Extracting Connected Components: "<<world.rank()<<std::endl; 

  // Obtain geometries

  // Get disjoint sets of element IDs
  // b = static_cast<Block_Feature*> (master.get(0)); // load block with local id 0
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

  if(connected_components_str.size() > 0) {

    // std::vector<std::set<element_t>> cc; // connected components 

    // Convert element IDs to elements
    for(auto& comp_str : connected_components_str) {
      // std::set<element_t>& comp = cc.emplace_back(); 
      std::vector<element_t> eles; 

      for(auto& eid : comp_str) {
        // comp.insert(hypermesh::regular_simplex_mesh_element(m, 0, eid)); 

        hypermesh::regular_simplex_mesh_element f(m, 0, eid);
        eles.push_back(f); 
      }

      std::sort(eles.begin(), eles.end(), [&](auto&a, auto& b) {
        return (a.corner[2] < b.corner[2]) || (a.corner[2] == b.corner[2] && a.corner[0] < b.corner[0]) || (a.corner[2] == b.corner[2] && a.corner[0] == b.corner[0] && a.corner[1] < b.corner[1]); 
      }); 

      std::vector<float>& level_set = level_sets.emplace_back();
      for (int k = 0; k < eles.size(); k ++) {
        auto& p = intersections->at(eles[k].to_string());

        level_set.push_back(p.x[0]); //  / (DW-1));
        level_set.push_back(p.x[1]); //  / (DH-1));
        level_set.push_back(p.x[2]); //  / (DT-1));
        level_set.push_back(p.val);
      } 
    }

  }

  // #ifdef FTK_HAVE_MPI
  //   #if TIME_OF_STEPS
  //     MPI_Barrier(world);
  //     end = MPI_Wtime();
  //     if(world.rank() == 0) {
  //       std::cout << "Generate level_sets: " << end - start << " seconds. " << std::endl;
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

void print_level_sets()
{
  printf("We found %lu level_sets:\n", level_sets.size());
  for (int i = 0; i < level_sets.size(); i ++) {
    printf("--Level-set %d:\n", i);
    const auto &level_set = level_sets[i];
    for (int k = 0; k < level_set.size()/4; k ++) {
      printf("---x=(%f, %f), t=%f, val=%f\n", level_set[k*4], level_set[k*4+1], level_set[k*4+2], level_set[k*4+3]);
    }
  }
}

void read_level_set_file(const std::string& f)
{
  std::ifstream ifs(f, std::ios::in | std::ios::binary);

  while(!ifs.eof()){
    float ncoord;
    ifs.read(reinterpret_cast<char*>(&ncoord), sizeof(float)); 

    // std::cout<<ncoord<<std::endl;

    auto& level_set = level_sets.emplace_back(); 
    for(int j = 0; j < (int) ncoord; ++j) {
      float coord;
      ifs.read(reinterpret_cast<char*>(&coord), sizeof(float)); 

      level_set.push_back(coord); 
    }
  }

  ifs.close();
}

void write_level_set_file(diy::mpi::communicator& world, const std::string& f)
{
  std::vector<float> buf; 

  for(auto& level_set : level_sets) {
    buf.push_back(level_set.size()); 

    // std::cout<<level_set.size()<<std::endl;

    for(auto& coord : level_set) {
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

void read_dump_file(const std::string& f)
{
  std::vector<intersection_t> vector;

  std::ifstream ifs(f);
  cereal::JSONInputArchive ar(ifs);
  ar(vector);
  ifs.close();

  for (const auto &i : vector) {
    hypermesh::regular_simplex_mesh_element e(m, 0, i.eid);
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
  std::string filename_level_set_r, filename_level_set_w;
  std::string filename_sets_w; 

  bool show_qt = false, show_vtk = false;

  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "input file name pattern", cxxopts::value<std::string>(pattern))
    ("f,format", "input file format", cxxopts::value<std::string>(format)->default_value("float32"))
    ("read-dump", "read dump file", cxxopts::value<std::string>(filename_dump_r))
    ("write-dump", "write dump file", cxxopts::value<std::string>(filename_dump_w))
    ("read-level-set", "read level_set file", cxxopts::value<std::string>(filename_level_set_r))
    ("write-level-set", "write level_set file", cxxopts::value<std::string>(filename_level_set_w))
    ("write-sets", "write sets of connected elements", cxxopts::value<std::string>(filename_sets_w))
    ("write-time-union-find", "write time for each round of union-find", cxxopts::value<std::string>(filename_time_uf_w))
    ("w,width", "width", cxxopts::value<int>(DW)->default_value("128"))
    ("h,height", "height", cxxopts::value<int>(DH)->default_value("128"))
    ("t,timesteps", "timesteps", cxxopts::value<int>(DT)->default_value("10"))
    ("critical-point-type", "Track which type of critical points", cxxopts::value<int>(CRITICAL_POINT_TYPE)->default_value(std::to_string(MAXIMUM_POINT)))
    ("scaling-factor", "scaling factor for synthetic data", cxxopts::value<int>(scaling_factor)->default_value("15"))
    ("threshold", "threshold", cxxopts::value<float>(threshold)->default_value("0"))
    // ("threshold-length", "threshold for level_set length", cxxopts::value<int>(threshold_length)->default_value("-1"))
    ("vtk", "visualization with vtk", cxxopts::value<bool>(show_vtk))
    ("qt", "visualization with qt", cxxopts::value<bool>(show_qt))
    ("d,debug", "enable debugging");
  auto results = options.parse(argc, argv);

  #ifdef FTK_HAVE_MPI
    start = MPI_Wtime();
  #endif

  // Decompose mesh
  // ========================================

  m.set_lb_ub({0, 0, 0}, {DW-1, DH-1, DT-1}); // update the mesh; set the lower and upper bounds of the mesh

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

  #if PRINT_FEATURE_DENSITY || PRINT_SIMPLEX_COUNT
    int element_cnt = 0; 

    m.element_for(0, [&](const hypermesh::regular_simplex_mesh_element& f){
      element_cnt++ ;
    }, nthreads);
  #endif

  #if PRINT_SIMPLEX_COUNT
    if (world.rank() == 0) {
      std::cout<<"Simplex Count: "<< element_cnt << std::endl; 
      exit(0); 
    }

    MPI_Barrier(world);
    end = MPI_Wtime();
    start = end; 
  #endif
  
  if (!filename_level_set_r.empty()) { // if the level_set file is given, skip all the analysis and visualize/print the level_sets
    if(world.rank() == 0) {
      read_level_set_file(filename_level_set_r);
    }
  } else { // otherwise do the analysis

    if (!filename_dump_r.empty()) { // if the dump file is given, skill the sweep step; otherwise do sweep-and-trace
      read_dump_file(filename_dump_r);
    } else { // derive gradients and do the sweep

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
      //       std::cout << "Scan for Satisfied Points: " << end - start << " seconds. " <<std::endl;
      //     }
      //     start = end; 
      //   #endif
      // #endif

      // std::cout<<"Finish scanning: "<<world.rank()<<std::endl; 
    }

    add_related_elements_to_intersections(intersections, m, block_m, FEATURE_DIM); 

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

    add_points_to_block(m, block_m, b, FEATURE_DIM); 

    // #ifdef FTK_HAVE_MPI
    //   #if TIME_OF_STEPS
    //     MPI_Barrier(world);
    //     end = MPI_Wtime();
    //     if(world.rank() == 0) {
    //       std::cout << "Add points into data block: " << end - start << " seconds. " <<std::endl;
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
      if(world.size() > 1) {
        load_balancing_resize_bounds(world, master, assigner, m, gid, b); 
      }
    #endif

    #ifdef FTK_HAVE_MPI
      #if TIME_OF_STEPS
        MPI_Barrier(world);
        end = MPI_Wtime();
        if(world.rank() == 0) {
            std::cout << "UF: Load balancing - Resize Bounds: " << end - start << " seconds. " << std::endl;
        }
        start = end; 
      #endif
    #endif

    #if LOAD_BALANCING
      if(world.size() > 1) {
        load_balancing_redistribute_data(world, master, assigner, m, gid, b, FEATURE_DIM); 
      }
    #endif

    #ifdef FTK_HAVE_MPI
      #if TIME_OF_STEPS
        MPI_Barrier(world);
        end = MPI_Wtime();
        if(world.rank() == 0) {
            std::cout << "UF: Load balancing - Redistribute Data: " << end - start << " seconds. " << std::endl;
        }
        start = end; 
      #endif
    #endif

    #if PRINT_EDGE_COUNT
      int edge_cnt = 0;

      int local_edge_cnt = 0; 
      for(auto& feature : b->features) {
        local_edge_cnt += feature.related_elements.size();
      }

      if (world.rank() == 0) {
        diy::mpi::reduce<int>(world, local_edge_cnt, edge_cnt, 0, std::plus<int>());
      } else {
        diy::mpi::reduce<int>(world, local_edge_cnt, 0, std::plus<int>());
      }

      if (world.rank() == 0) {
        std::cout<<"Edge Count is " << edge_cnt << std::endl; 
      }
      
      MPI_Barrier(world);
      end = MPI_Wtime();
      start = end; 
    #endif

    #if LOAD_BALANCING
      if(world.size() > 1) {
        init_block_after_load_balancing(world, master, assigner, m, gid, b, FEATURE_DIM); 
        master.clear();  // clear the added block, since we will add block into master again. 
      } else {
        init_block_without_load_balancing(lattice_partitions, m, gid, b, FEATURE_DIM);  
      }
    #else
      init_block_without_load_balancing(lattice_partitions, m, gid, b, FEATURE_DIM);
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
      if (!filename_level_set_w.empty()) {
        write_level_set_file(world, filename_level_set_w);

        #ifdef FTK_HAVE_MPI
          #if TIME_OF_STEPS
            MPI_Barrier(world);
            end = MPI_Wtime();
            if(world.rank() == 0) {
              std::cout << "Output level_sets: " << end - start << " seconds. " << std::endl;
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

  return 0;
}
