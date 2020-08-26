#define DIM 4
#define FEATURE_DIM DIM-1

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

#define LOAD_BALANCING true

#define TIME_OF_STEPS true
#define MULTITHREAD false

#define PRINT_FEATURE_DENSITY false
#define PRINT_ELE_COUNT true

#define PRINT_ROUND_COUNT_AND_PEAK_MEMORY true

#define ALL_CRITICAL_POINT 0
#define MAXIMUM_POINT      1


int CRITICAL_POINT_TYPE; // By default, track maximum points

int nthreads;

int DW, DH, DD; // the dimensionality of the data is DW*DH*DD
int DT; // number of timesteps
std::vector<int> D_sizes(DIM);

double start, end; 

hypermesh::ndarray<float> scalar, grad, hess;
diy::DiscreteBounds data_box(DIM); // bounds of data box
std::vector<int> data_offset(DIM);

hypermesh::regular_simplex_mesh m(DIM); // the 3D space-time mesh

hypermesh::regular_simplex_mesh block_m(DIM); // the 3D space-time mesh of this block
hypermesh::regular_simplex_mesh block_m_ghost(DIM); // the 3D space-time mesh of this block with ghost cells

std::vector<std::tuple<hypermesh::regular_lattice, hypermesh::regular_lattice>> lattice_partitions;
float threshold; // threshold for trajectories. The max scalar on each trajectory should be larger than the threshold. 
int threshold_length; // threshold for trajectory length.  

std::mutex mutex;
 
Block_Feature* b = new Block_Feature(); 
int gid; // global id of this block / this process
std::map<std::string, intersection_t>* intersections;

// the output sets of connected elements
std::vector<std::set<std::string>> connected_components_str; // connected components 

// the output trajectories
std::vector<std::vector<float>> trajectories;

std::string filename_time_uf_w; // record time for each round of union-find

template <typename T>
hypermesh::ndarray<T> derive_gradients3(const hypermesh::ndarray<T>& scalar)
{
  hypermesh::ndarray<T> grad;
  std::vector<size_t> reshape(DIM+1);
  reshape[0] = DIM-1;
  for(int i = 0; i < DIM; ++i) {
    reshape[i+1] = scalar.dim(i); 
  }
  grad.reshape(reshape);
  
  for (int l = std::max(0, block_m_ghost.lb(3)-1) - data_offset[3]; l < std::min(block_m_ghost.ub(3)+2, DT) - data_offset[3]; l ++) {
    for (int k = std::max(0, block_m_ghost.lb(2)-1) - data_offset[2]; k < std::min(block_m_ghost.ub(2)+2, DD-1) - data_offset[2]; k ++) {
      for (int j = std::max(1, block_m_ghost.lb(1)-1) - data_offset[1]; j < std::min(block_m_ghost.ub(1)+2, DH-1) - data_offset[1]; j ++) {
        for (int i = std::max(1, block_m_ghost.lb(0)-1) - data_offset[0]; i < std::min(block_m_ghost.ub(0)+2, DW-1) - data_offset[0]; i ++) {

          grad(0, i, j, k, l) = 0.5 * (scalar(i+1, j, k, l) - scalar(i-1, j, k, l)) * (DW-1);
          grad(1, i, j, k, l) = 0.5 * (scalar(i, j+1, k, l) - scalar(i, j-1, k, l)) * (DH-1);
          grad(2, i, j, k, l) = 0.5 * (scalar(i, j, k+1, l) - scalar(i, j, k-1, l)) * (DD-1);
        }
      }
    }
  }

  return grad;
}

template <typename T>
hypermesh::ndarray<T> derive_hessians3(const hypermesh::ndarray<T>& grad)
{
  hypermesh::ndarray<T> hess;
  std::vector<size_t> reshape(DIM+2);
  reshape[0] = DIM-1;
  for(int i = 0; i < DIM+1; ++i) {
    reshape[i+1] = grad.dim(i); 
  }

  for (int l = std::max(0, block_m_ghost.lb(3)) - data_offset[3]; l < std::min(block_m_ghost.ub(3)+1, DT) - data_offset[3]; l ++) {
    for (int k = std::max(2, block_m_ghost.lb(2)) - data_offset[2]; k < std::min(block_m_ghost.ub(2)+1, DD-2) - data_offset[2]; k ++) {
      for (int j = std::max(2, block_m_ghost.lb(1)) - data_offset[1]; j < std::min(block_m_ghost.ub(1)+1, DH-2) - data_offset[1]; j ++) {
        for (int i = std::max(2, block_m_ghost.lb(0)) - data_offset[0]; i < std::min(block_m_ghost.ub(0)+1, DW-2) - data_offset[0]; i ++) {

          const T H00 = hess(0, 0, i, j, k) = // ddf/dx2
            0.5 * (grad(0, i+1, j, k) - grad(0, i-1, j, k)) * (DW-1);
          const T H01 = hess(0, 1, i, j, k) = // ddf/dxdy
            0.5 * (grad(0, i, j+1, k) - grad(0, i, j-1, k)) * (DH-1);
          const T H02 = hess(0, 2, i, j, k) = // ddf/dxdz
            0.5 * (grad(0, i, j, k+1) - grad(0, i, j, k-1)) * (DD-1);

          const T H10 = hess(1, 0, i, j, k) = // ddf/dydx
            0.5 * (grad(1, i+1, j, k) - grad(1, i-1, j, k)) * (DW-1);
          const T H11 = hess(1, 1, i, j, k) = // ddf/dy2
            0.5 * (grad(1, i, j+1, k) - grad(1, i, j-1, k)) * (DH-1);
          const T H12 = hess(1, 2, i, j, k) = // ddf/dydz
            0.5 * (grad(1, i, j, k+1) - grad(1, i, j, k-1)) * (DD-1);

          const T H20 = hess(2, 0, i, j, k) = // ddf/dzdx
            0.5 * (grad(2, i+1, j, k) - grad(2, i-1, j, k)) * (DW-1);
          const T H21 = hess(2, 1, i, j, k) = // ddf/dzdy
            0.5 * (grad(2, i, j+1, k) - grad(2, i, j-1, k)) * (DH-1);
          const T H22 = hess(2, 2, i, j, k) = // ddf/dz2
            0.5 * (grad(2, i, j, k+1) - grad(2, i, j, k-1)) * (DD-1);

        }
      }
    }
  }

  return hess;
}

void decompose_mesh(int nblocks) {
  std::vector<size_t> given = {0}; // partition the 2D spatial space and 1D timespace
  // std::vector<size_t> given = {0, 0, 1}; // Only partition the 2D spatial space
  // std::vector<size_t> given = {1, 1, 0}; // Only partition the 1D temporal space

  std::vector<size_t> ghost(DIM); // at least 1, larger is ok
  for(int i = 0; i < DIM; ++i) {
    ghost[i] = 1; 
  }

  const hypermesh::regular_lattice& lattice = m.lattice(); 
  lattice.partition(nblocks, given, ghost, lattice_partitions); 

  auto& lattice_pair = lattice_partitions[gid]; 
  hypermesh::regular_lattice& lattice_p = std::get<0>(lattice_pair); 
  hypermesh::regular_lattice& lattice_ghost_p = std::get<1>(lattice_pair); 
  
  block_m.set_lb_ub(lattice_p);
  block_m_ghost.set_lb_ub(lattice_ghost_p);

  for(int i = 0; i < DIM; ++i) {
    data_box.min[i] = std::max(0, block_m_ghost.lb(i)-2); data_box.max[i] = std::min(block_m_ghost.ub(i)+2, D_sizes[i]-1);  
  }

  for(int i = 0; i < DIM; ++i) {
    data_offset[i] = data_box.min[i]; 
  }
}

void check_simplex(const hypermesh::regular_simplex_mesh_element& f)
{
  if (!f.valid()) return; // check if the 3-simplex is valid
  const auto &vertices = f.vertices(); // obtain the vertices of the simplex
  float g[4][3], value[4];

  for (int i = 0; i < 4; i ++) {
    int _i = vertices[i][0] - data_offset[0];
    int _j = vertices[i][1] - data_offset[1];
    int _k = vertices[i][2] - data_offset[2];
    int _l = vertices[i][3] - data_offset[3];

    g[i][0] = grad(0, _i, _j, _k, _l);
    g[i][1] = grad(1, _i, _j, _k, _l);
    g[i][2] = grad(2, _i, _j, _k, _l);
    value[i] = scalar(_i, _j, _k, _l);
  }
 
  // check intersection
  float mu[4];
  bool succ = ftk::inverse_lerp_s3v3(g, mu);
  
  if (!succ) return;

  if(CRITICAL_POINT_TYPE == MAXIMUM_POINT) {

    // check hessian
    float H[4][3][3], h[3][3];
    for (int i = 0; i < 4; i ++)
      for (int j = 0; j < 3; j ++)
        for (int k = 0; k < 3; k ++)
          H[i][j][k] = hess(j, k, vertices[i][0], vertices[i][1], vertices[i][2], vertices[i][3]);
    ftk::lerp_s3m3x3(H, mu, h);

    float eig[3];
    ftk::solve_eigenvalues_symmetric3x3(h, eig);
  
    if (eig[0] < 0 && eig[1] < 0 && eig[2] < 0) { // local maxima
      
    } else {
      return ;
    }

  }

  // A critical point is detected

  float X[4][4];
  for (int i = 0; i < vertices.size(); i ++)
    for (int j = 0; j < DIM; j ++)
      X[i][j] = vertices[i][j];

  float x[4];
  ftk::lerp_s3v4(X, mu, x);

  intersection_t I;
  I.eid = f.to_string();
  I.val = ftk::lerp_s3(value, mu);
  // I.x.resize(DIM); 

  for(int i = 0; i < DIM; ++i) {
    I.x[i] = x[i]; 
    // I.corner.push_back(f.corner[i]); 
  }

  {
    #if MULTITHREAD 
      std::lock_guard<std::mutex> guard(mutex);
    #endif

    intersections->insert(std::make_pair(I.eid, I)); 
    // fprintf(stderr, "x={%f, %f}, t=%f, val=%f\n", I.x[0], I.x[1], I.x[2], I.x[3], I.val);
  }

}

void scan_intersections() 
{
  block_m_ghost.element_for(DIM-1, check_simplex, nthreads); // iterate over all 2-simplices
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
      start = MPI_Wtime(); 
    #endif

    #if PRINT_ROUND_COUNT_AND_PEAK_MEMORY
      int total_rount_cnt = 0; 
      MPI_Reduce(&b->nrounds, &total_rount_cnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 

      double total_peak_memory = 0; 
      MPI_Reduce(&b->peak_memory, &total_peak_memory, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 

      if(world.rank() == 0) {
        std::cout << "\t Number of Iterations per Process: " << total_rount_cnt / (double) world.size() << std::endl;
        std::cout << "\t Peak Memory per Process: " << total_peak_memory / (double) world.size() << std::endl;
      }

      MPI_Barrier(world);
      start = MPI_Wtime();
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

  // if(world.rank() == 0) {
  if(connected_components_str.size() > 0) {
    
    std::vector<std::set<element_t>> cc; // connected components 
    // Convert element IDs to elements
    for(auto& comp_str : connected_components_str) {
      std::set<element_t>& comp = cc.emplace_back(); 
      for(auto& eid : comp_str) {
        comp.insert(hypermesh::regular_simplex_mesh_element(m, DIM-1, eid)); 
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
            
            for(int l = 0; l < DIM; ++l) {
              mycurve.push_back(p.x[l]); 
            }

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
    for (int k = 0; k < curve.size()/(DIM+1); k ++) {
      // printf("---x=(%f, %f), t=%f, val=%f\n", curve[k*4], curve[k*4+1], curve[k*4+2], curve[k*4+3]);
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
  ftk::write_curves_vtk(trajectories, f, DIM+1); 
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
    hypermesh::regular_simplex_mesh_element e(m, DIM-1, i.eid);
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
  auto vtkcurves = ftk::curves2vtk(trajectories, DIM+1);
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
    ("d,depth", "depth", cxxopts::value<int>(DD)->default_value("128"))
    ("t,timesteps", "timesteps", cxxopts::value<int>(DT)->default_value("10"))
    ("critical-point-type", "Track which type of critical points", cxxopts::value<int>(CRITICAL_POINT_TYPE)->default_value(std::to_string(MAXIMUM_POINT)))
    ("threshold", "threshold", cxxopts::value<float>(threshold)->default_value("0"))
    ("threshold-length", "threshold for trajectory length", cxxopts::value<int>(threshold_length)->default_value("-1"))
    ("vtk", "visualization with vtk", cxxopts::value<bool>(show_vtk))
    ("qt", "visualization with qt", cxxopts::value<bool>(show_qt)); 
  auto results = options.parse(argc, argv);

  #ifdef FTK_HAVE_MPI
    start = MPI_Wtime();
  #endif

  // Decompose mesh
  // ========================================

  // ?????
  D_sizes = {DW, DH, DD, DT}; 
  m.set_lb_ub({0, 0, 0, 0}, {DW-1, DH-1, DD-1, DT-1}); // update the mesh; set the lower and upper bounds of the mesh

  decompose_mesh(nblocks); 

  intersections = &b->intersections; 

  // Load data
  // ========================================

  if (pattern.empty()) { // if the input data is not given, generate a synthetic data for the demo
    std::cout<<"Error! No Data! "<<std::endl; 
    exit(0);
  } else { // load the binary data

    diy::mpi::io::file in(world, pattern, diy::mpi::io::file::rdonly);
    
    diy::DiscreteBounds diy_box(DIM);
    for(int i = 0; i < DIM; ++i) {
      diy_box.min[i] = data_box.min[DIM-1-i]; diy_box.max[i] = data_box.max[DIM-1-i];  
    }

    std::vector<unsigned> shape(DIM);
    for(int i = 0; i < DIM; ++i) {
      shape[i] = D_sizes[DIM-1-i]; 
    }

    diy::io::BOV reader(in, shape);

    std::vector<size_t> reshape(DIM);
    for(int i = 0; i < DIM; ++i) {
      reshape[i] = data_box.max[i] - data_box.min[i] + 1; 
    }
    scalar.reshape(reshape);

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
          if(ite->size() / (DIM+1) <= threshold_length) {
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
          
          for(int i = DIM; i < ite->size(); i += (DIM+1)) {
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
      grad = derive_gradients3(scalar);
      hess = derive_hessians3(grad);

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
        std::cout<<"\t Feature Element Count is " << feature_ele_cnt << std::endl; 
      }
    #endif

    #if PRINT_FEATURE_DENSITY
      int element_cnt = 0; 

      m.element_for(DIM-1, [&](const hypermesh::regular_simplex_mesh_element& f){
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
      hypermesh::ndarray<float> full_scalar; full_scalar.reshape(DW, DH, DD, DT);
      
      for (int l = data_box.min[3]; l < data_box.max[3]+1; l ++) {
        for (int k = data_box.min[2]; k < data_box.max[2]+1; k ++) {
          for (int j = data_box.min[1]; j < data_box.max[1]+1; j ++) {
            for (int i = data_box.min[0]; i < data_box.max[0]+1; i ++) {
              full_scalar(i, j, k, l) = scalar(i - data_offset[0], j - data_offset[1], k - data_offset[2], l - data_offset[3]); 
            }
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
