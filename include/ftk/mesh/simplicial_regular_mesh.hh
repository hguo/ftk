#ifndef _HYPERMESH_simplicial_regular_MESH_HH
#define _HYPERMESH_simplicial_regular_MESH_HH

#include <ftk/ftk_config.hh>
#include <iostream>
#include <sstream> 
#include <vector>
#include <tuple>
#include <set>
#include <string>
#include <algorithm>
#include <numeric>
#include <thread>
#include <functional>
#include <cassert>
#include <iterator>
#include <functional>
#include <ftk/mesh/lattice.hh>
#include <ftk/external/diy/serialization.hpp>

#if FTK_HAVE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#if FTK_HAVE_TBB
#include <tbb/tbb.h>
#endif

namespace ftk {

enum {
  ELEMENT_SCOPE_ALL = 0,
  ELEMENT_SCOPE_ORDINAL = 1, 
  ELEMENT_SCOPE_INTERVAL = 2
};

struct simplicial_regular_mesh;

struct simplicial_regular_mesh_element {
  friend class simplicial_regular_mesh;

  simplicial_regular_mesh_element() : dim(0), type(0) {}
  simplicial_regular_mesh_element(const simplicial_regular_mesh_element& e) : corner(e.corner), dim(e.dim), type(e.type) {}
  simplicial_regular_mesh_element(int nd, int d);
  simplicial_regular_mesh_element(const std::vector<int>& corner, int d, int type); // nd is encoded in corner
  // simplicial_regular_mesh_element(const simplicial_regular_mesh &m, int d, size_t index);
  simplicial_regular_mesh_element(const simplicial_regular_mesh &m, int d, size_t work_index, 
      const lattice& l, int scope = ELEMENT_SCOPE_ALL);
  template <typename uint=uint64_t> simplicial_regular_mesh_element(const simplicial_regular_mesh& m, int d, uint i);
  simplicial_regular_mesh_element(const std::string &i);

  simplicial_regular_mesh_element& operator=(const simplicial_regular_mesh_element& e);
  bool operator!=(const simplicial_regular_mesh_element& e) const {return !(*this == e);}
  bool operator<(const simplicial_regular_mesh_element& e) const;
  bool operator==(const simplicial_regular_mesh_element& e) const;
  // simplicial_regular_mesh_element& operator++();
  friend std::ostream& operator<<(std::ostream& os, const simplicial_regular_mesh_element&);
  void print(std::ostream& os, const simplicial_regular_mesh& m) const;

  std::vector<std::vector<int> > vertices(const simplicial_regular_mesh&) const;

  // void increase_corner(int d=0);
  bool valid(const simplicial_regular_mesh& m) const; // validate the element in the context of the given mesh

  // template <int nd> std::tuple<std::array<int, nd>, int> to_index() const;
  // template <int nd> void from_index(const std::tuple<std::array<int, nd>, int>&);

  size_t to_work_index(const simplicial_regular_mesh& m, const lattice& l, int scope = ELEMENT_SCOPE_ALL) const;
  void from_work_index(const simplicial_regular_mesh& m, size_t, const lattice& l, int scope = ELEMENT_SCOPE_ALL);

  template <typename uint = uint64_t> uint to_integer(const simplicial_regular_mesh& m) const;
  template <typename uint = uint64_t> void from_integer(const simplicial_regular_mesh& m, uint i);

  // std::string to_string() const;
  // void from_string(const std::string& index);

  template <typename Archive> void serialize(Archive &ar) {ar(dim, type, corner);}

  std::vector<simplicial_regular_mesh_element> sides(const simplicial_regular_mesh& m) const;
  std::vector<simplicial_regular_mesh_element> side_of(const simplicial_regular_mesh& m) const;
  
  bool is_ordinal(const simplicial_regular_mesh& m) const;

  // const simplicial_regular_mesh &m; // avoid the ref to mesh to ease (de)serialization
  std::vector<int> corner;
  int dim, type;
};


struct simplicial_regular_mesh {
  friend class simplicial_regular_mesh_element;
  typedef simplicial_regular_mesh_element iterator;

  simplicial_regular_mesh(int n) : nd_(n), lattice_(n) {
    for (int i = 0; i < n; i ++) {
      lb_.push_back(0);
      ub_.push_back(0);
    }

    initialize_subdivision();

    lattice_.reshape(lb_, sizes()); 
  }

  simplicial_regular_mesh(int n, lattice& _lattice) : simplicial_regular_mesh(n) {
    std::vector<int> _lb, _ub;

    for (int i = 0; i < n; i ++) {
      int _start = _lattice.start(i); 
      int _size = _lattice.size(i); 

      _lb.push_back(_start);
      _ub.push_back(_start + _size - 1);
    }

    // lattice_.set_unlimited_time(_lattice.unlimited_time()); 
    set_lb_ub(_lb, _ub); 
  }

  // Dimensionality of the mesh
  int nd() const {return nd_;}

  // Number of unique types of d-simplices
  int ntypes(int d, int scope = ELEMENT_SCOPE_ALL) const;

  // Number of d-dimensional elements
  size_t n(int d, int scope = ELEMENT_SCOPE_ALL) const;

  // number of d-dimensional ordinal and interval elements
  size_t n_ordinal(int d) const;
  size_t n_interval(int d) const;

  // Returns d+1 vertices that build up the d-dimensional simplex of the given type
  std::vector<std::vector<int>> unit_simplex(int d, int t) const {return unit_simplices[d][t];}

  void get_lb_ub(std::vector<int>& lb, std::vector<int>& ub) {lb = lb_; ub = ub_;}
  void set_lb_ub(const std::vector<int>& lb, const std::vector<int>& ub);
  void set_lb_ub(const lattice& lattice); 
  int lb(int d) const {return lb_[d];}
  int ub(int d) const {return ub_[d];}
  const std::vector<int>& lb() const {return lb_;}
  const std::vector<int>& ub() const {return ub_;}
  std::vector<int> sizes();  

  const lattice& get_lattice() const {return lattice_; }

  iterator element_begin(int d, int scope = ELEMENT_SCOPE_ALL);
  iterator element_end(int d, int scope = ELEMENT_SCOPE_ALL);

  // Iterate all d-simplices with the given f() function using the given number of threads.
  // There are three different modes: all, fixed time, and fixed interval.  Please
  // see the comments in the next two functions for more details.
  void element_for(int d, std::function<void(simplicial_regular_mesh_element)> f, 
      int nthreads=std::thread::hardware_concurrency()); // , int scope = ELEMENT_SCOPE_ALL);

  void element_for_ordinal(int d, int t, std::function<void(simplicial_regular_mesh_element)> f,
      int nthreads=std::thread::hardware_concurrency());

  void element_for_interval(int d, int t0, int t1, std::function<void(simplicial_regular_mesh_element)> f,
      int nthreads=std::thread::hardware_concurrency());
  
  void element_for(int d, const lattice& subdomain, int scope, 
      std::function<void(simplicial_regular_mesh_element)> f,
      int nthreads=std::thread::hardware_concurrency());

#if 0
public: // partitioning
  void partition(int np, std::vector<std::tuple<simplicial_regular_mesh, simplicial_regular_mesh>>& partitions);  
  void partition(int np, const std::vector<size_t> &given, std::vector<std::tuple<simplicial_regular_mesh, simplicial_regular_mesh>>& partitions);  
  void partition(int np, const std::vector<size_t> &given, const std::vector<size_t> &ghost, std::vector<std::tuple<simplicial_regular_mesh, simplicial_regular_mesh>>& partitions);  
  void partition(int np, const std::vector<size_t> &given, const std::vector<size_t> &ghost_low, const std::vector<size_t> &ghost_high, std::vector<std::tuple<simplicial_regular_mesh, simplicial_regular_mesh>>& partitions);  
#endif

public:
  void print_unit_simplices(int nd, int d) const;

private: // initialization functions
  void initialize_subdivision();

  // Generate spatial subdivision for the unit n-cube
  // input: dimensionality n
  // output: a vector (n!) of n-simplices; each simplex is a vector of vertices; 
  std::vector<std::vector<std::vector<int>>> subdivide_unit_cube(int n); 

  // Reduce a simplex in the unit cube.  The reduction means to re-encode the simplex with 
  // the corner coordinate (offset) and the node coordinates of the reduced simplex.  
  // For example in a 2-cube, simplex (01,11) can be represented by the cube at (0,1) and 
  // with vertices (00,10)
  // input: a k-dimensional simplex inside the unit cube
  // output: 1) simplex: the sample simplex that is encoded by the 
  //            original cube or another cube. 
  //         2) offset: the offset is all-zero if the simplex can be 
  //            encoded by the original cube, otherwise is the coordinate 
  //            of the cube that uniquely encodes the simplex
  std::tuple<std::vector<std::vector<int>>, std::vector<int>> reduce_unit_simplex(const std::vector<std::vector<int>> &simplex) const;

  // Enumerate all types reduced k-simplices in the n-cube
  std::vector<std::vector<std::vector<int>>> enumerate_unit_simplices(int n, int k);

  // Enumerate all sides of the specific type of k-simplex
  std::vector<std::tuple<int, std::vector<int>>> enumerate_unit_simplex_sides(int k, int type);
 
  // Enumerate all element that contains the specific type of k-simplex
  std::vector<std::tuple<int, std::vector<int>>> enumerate_unit_simplex_side_of(int k, int type);

  void derive_ordinal_and_interval_simplices();

  // bool is_simplex_identical(const std::vector<std::string>&, const std::vector<std::string>&) const;

private:
  const int nd_;
  std::vector<int> lb_, ub_; // lower and upper bounds of each dimension
  std::vector<int> ntypes_, ntypes_ordinal_, ntypes_interval_; // number of types for k-simplex
  std::vector<int> dimprod_;

  struct lattice lattice_; 

  // list of k-simplices types; each simplex contains k vertices
  // unit_simplices[d][type] retunrs d+1 vertices that build up the simplex
  std::vector<std::vector<std::vector<std::vector<int>>>> unit_simplices;

  std::vector<std::vector<int>> unit_ordinal_simplex_types, 
                                unit_interval_simplex_types;
  
  std::vector<std::vector<bool>> is_unit_simplex_type_ordinal;
  
  // (dim,type) --> vector of (type,offset)
  std::vector<std::vector<std::vector<std::tuple<int, std::vector<int>>>>> unit_simplex_sides;
  std::vector<std::vector<std::vector<std::tuple<int, std::vector<int>>>>> unit_simplex_side_of;
};



//////////////////////////////////
inline simplicial_regular_mesh_element::simplicial_regular_mesh_element(int nd, int d_) 
  : dim(d_), type(0) 
{
  corner.resize(nd);
}

inline simplicial_regular_mesh_element::simplicial_regular_mesh_element(
    const std::vector<int>& corner_, int d_, int type_)
  : corner(corner_), dim(d_), type(type_)
{
}

template <typename uint>
inline simplicial_regular_mesh_element::simplicial_regular_mesh_element(const simplicial_regular_mesh& m, int d, uint i)
{
  corner.resize(m.nd()); dim = d; from_integer(m, i);
}

inline simplicial_regular_mesh_element::simplicial_regular_mesh_element(
    const simplicial_regular_mesh &m, int d_, size_t i, const lattice& l, int scope)
{
  dim = d_;
  corner.resize(m.nd());
  from_work_index(m, i, l, scope);
}

#if 0
inline simplicial_regular_mesh_element::simplicial_regular_mesh_element(const std::string &i)
{
  from_string(i);
}
#endif

inline simplicial_regular_mesh_element& simplicial_regular_mesh_element::operator=(const simplicial_regular_mesh_element& e)
{
  corner = e.corner;
  dim = e.dim;
  type = e.type;
  return *this;
}

inline bool simplicial_regular_mesh_element::operator<(const simplicial_regular_mesh_element& e) const
{
  if (corner < e.corner) return true;
  else if (corner == e.corner) return type < e.type;
  else return false;
}

inline bool simplicial_regular_mesh_element::operator==(const simplicial_regular_mesh_element& e) const
{
  return dim == e.dim && type == e.type && corner == e.corner;
}

inline bool simplicial_regular_mesh_element::valid(const simplicial_regular_mesh& m) const {
  if (type < 0 || type >= m.ntypes(dim)) return false;
  else {
    auto my_vertices = vertices(m);
    for (int i = 0; i < my_vertices.size(); i ++)
      for (int j = 0; j < my_vertices[i].size(); j ++)
        if (my_vertices[i][j] < m.lb(j) || my_vertices[i][j] > m.ub(j))
          return false;
    return true;
#if 0
    for (int i = 0; i < m.nd(); i ++) 
      if (corner[i] < m.lb(i) || corner[i] > m.ub(i))
        return false;
#endif
  }
}

inline std::vector<std::vector<int> > simplicial_regular_mesh_element::vertices(const simplicial_regular_mesh& m) const
{
  std::vector<std::vector<int>> vertices;

  const auto unit_vertices = m.unit_simplex(dim, type);
  vertices.resize(unit_vertices.size());
  for (size_t i = 0; i < unit_vertices.size(); i ++) {
    for (int j = 0; j < m.nd(); j ++) {
      vertices[i].push_back(corner[j] + unit_vertices[i][j]);
    }
  }

  return vertices;
}
 
inline void simplicial_regular_mesh_element::print(std::ostream& os, const simplicial_regular_mesh& m) const
{
  const auto &e = *this;

  os << "dim=" << e.dim << ",cornor={";
  for (size_t i = 0; i < e.corner.size(); i ++)
    if (i < e.corner.size()-1) os << e.corner[i] << ",";
    else os << e.corner[i];
  os << "},type=" << e.type << "/{";

  const auto s_vertices = m.unit_simplex(e.dim, e.type);
  for (size_t i = 0; i < s_vertices.size(); i ++) {
    for (size_t j = 0; j < s_vertices[i].size(); j ++)
      os << s_vertices[i][j];
    if (i < s_vertices.size()-1) os << ",";
  }
  os << "},";

  const auto vertices = e.vertices(m);
  os << "vertices={";
  for (size_t i = 0; i < vertices.size(); i ++) {
    os << "{";
    for (size_t j = 0; j < vertices[i].size(); j ++)
      if (j < vertices[i].size()-1) os << vertices[i][j] << ",";
      else os << vertices[i][j] << "}";;
    if (i < vertices.size()-1) os << ",";
    else os << "},";
  }

  // os << "int=" << e.to_integer(m) << ",";
  os << "valid=" << e.valid(m);
}

inline std::ostream& operator<<(std::ostream& os, const simplicial_regular_mesh_element& e)
{
  os << "dim=" << e.dim << ",cornor={";
  for (size_t i = 0; i < e.corner.size(); i ++)
    if (i < e.corner.size()-1) os << e.corner[i] << ",";
    else os << e.corner[i];
  os << "},type=" << e.type; // << "/{";

#if 0
  const auto s_vertices = m.unit_simplex(e.dim, e.type);
  for (size_t i = 0; i < s_vertices.size(); i ++) {
    for (size_t j = 0; j < s_vertices[i].size(); j ++)
      os << s_vertices[i][j];
    if (i < s_vertices.size()-1) os << ",";
  }
  os << "},";

  const auto vertices = e.vertices(m);
  os << "vertices={";
  for (size_t i = 0; i < vertices.size(); i ++) {
    os << "{";
    for (size_t j = 0; j < vertices[i].size(); j ++)
      if (j < vertices[i].size()-1) os << vertices[i][j] << ",";
      else os << vertices[i][j] << "}";;
    if (i < vertices.size()-1) os << ",";
    else os << "},";
  }

  os << "int=" << e.to_integer() << ",";
  os << "valid=" << e.valid();
#endif
 
#if 0
  os << ",";
  const auto sides = e.sides();
  os << "#sides=" << sides.size() << std::endl;

  for (auto side : sides)
    os << side;
#endif

#if 0
  os << ",";
  const auto side_of = e.side_of();
  os << "#side_of=" << side_of.size() << std::endl;

  for (auto e : side_of)
    os << e;
#endif
  
  return os;
}

inline size_t simplicial_regular_mesh_element::to_work_index(const simplicial_regular_mesh& m, const lattice& l, int scope) const
{
  size_t idx = l.to_integer(corner);
  return idx * m.ntypes(dim, scope);
}

inline void simplicial_regular_mesh_element::from_work_index(const simplicial_regular_mesh& m, size_t i, const lattice& l, int scope)
{
  const auto itype = i % m.ntypes(dim, scope);
  auto ii = i / m.ntypes(dim, scope);
 
  if (scope == ELEMENT_SCOPE_ORDINAL) type = m.unit_ordinal_simplex_types[dim][itype];
  else if (scope == ELEMENT_SCOPE_INTERVAL) type = m.unit_interval_simplex_types[dim][itype];
  else type = itype;

  corner = l.from_integer(ii);
  // fprintf(stderr, "i=%lu, itype=%lu, ii=%lu, ntypes=%d, scope=%d, corner=%d, %d, %d\n", 
  //     i, itype, ii, m.ntypes(dim, scope), scope, 
  //     corner[0], corner[1], corner[2]);
}

template <typename uint>
uint simplicial_regular_mesh_element::to_integer(const simplicial_regular_mesh& m) const
{
  uint corner_index = 0;
  for (size_t i = 0; i < m.nd(); i ++)
    corner_index += (corner[i] - m.lb(i)) * m.dimprod_[i];
  return corner_index * m.ntypes(dim) + type;
}

template <typename uint>
void simplicial_regular_mesh_element::from_integer(const simplicial_regular_mesh& m, uint index)
{
  type = index % m.ntypes(dim);
  uint corner_index = index / m.ntypes(dim); // m.dimprod_[m.nd()];

  for (int i = m.nd() - 1; i >= 0; i --) {
    corner[i] = corner_index / m.dimprod_[i]; 
    corner_index -= corner[i] * m.dimprod_[i];
  }
  for (int i = 0; i < m.nd(); i ++) 
    corner[i] += m.lb(i);

  // fprintf(stderr, "dimprod=%d, %d, %d, %d\n", m.dimprod_[0], m.dimprod_[1], m.dimprod_[2], m.dimprod_[3]);
  // fprintf(stderr, "lb=%d, %d, %d\n", m.lb(0), m.lb(1), m.lb(2));
  // std::cerr << "idx=" << index << "," << *this << std::endl;
}

#if 0
inline std::string simplicial_regular_mesh_element::to_string() const
{
  std::stringstream res;

  std::copy(corner.begin(), corner.end(), std::ostream_iterator<int>(res, ","));
  res<<type; 

  std::string res_str(res.str()); 

  return res_str; 
}

inline void simplicial_regular_mesh_element::from_string(const std::string& index)
{
  std::stringstream index_stream(index);
  std::vector<int> arr; 
  while( index_stream.good() )
  {
      std::string substr;
      getline( index_stream, substr, ',' );
      arr.push_back( std::stoi(substr) );
  }

  std::copy(arr.begin(), arr.end() - 1, corner.begin()); 
  // corner = std::vector<int>(arr.begin(), arr.end() - 1);
  
  type = *(arr.end()-1); 
}
#endif

#if 0
template <int nd> 
std::tuple<std::array<int, nd>, int> simplicial_regular_mesh_element::to_index() const
{
  std::array<int, nd> arr;
  std::copy_n(corner.begin(), nd, corner.end());
  return std::tuple<std::array<int, nd>, int>(arr, type);
}

template <int nd> 
void simplicial_regular_mesh_element::from_index(const std::tuple<std::array<int, nd>, int>& tuple)
{
  const std::array<int, nd> &arr = std::get<0>(tuple);
  corner = std::vector<int>(arr.begin(), arr.end());
  type = std::get<1>(tuple);
}
#endif

inline std::vector<simplicial_regular_mesh_element> simplicial_regular_mesh_element::sides(const simplicial_regular_mesh& m) const
{
  const auto &unit_simplex_sides = m.unit_simplex_sides[dim][type];
  std::vector<simplicial_regular_mesh_element> sides;

  for (auto s : unit_simplex_sides) {
    simplicial_regular_mesh_element side(corner, dim-1, std::get<0>(s));
    const auto &offset = std::get<1>(s);
    for (int i = 0; i < m.nd(); i ++) 
      side.corner[i] += offset[i];
    sides.push_back(side);
  }

  return sides;
}

inline std::vector<simplicial_regular_mesh_element> simplicial_regular_mesh_element::side_of(const simplicial_regular_mesh& m) const
{
  const auto &unit_simplex_side_of = m.unit_simplex_side_of[dim][type];
  std::vector<simplicial_regular_mesh_element> side_of;

  for (auto s : unit_simplex_side_of) {
    simplicial_regular_mesh_element e(corner, dim+1, std::get<0>(s));
    const auto &offset = std::get<1>(s);
    for (int i = 0; i < m.nd(); i ++)
      e.corner[i] += offset[i];
    side_of.push_back(e);
  }

  return side_of;
}

inline bool simplicial_regular_mesh_element::is_ordinal(const simplicial_regular_mesh& m) const
{
  return m.is_unit_simplex_type_ordinal[dim][type];
}

/////

inline int simplicial_regular_mesh::ntypes(int d, int scope) const 
{
  switch (scope) {
  case ELEMENT_SCOPE_ALL: return ntypes_.at(d);
  case ELEMENT_SCOPE_ORDINAL: return ntypes_ordinal_.at(d);
  case ELEMENT_SCOPE_INTERVAL: return ntypes_interval_.at(d);
  default: return 0;
  }
}

inline std::vector<std::vector<std::vector<int>>> simplicial_regular_mesh::subdivide_unit_cube(int n)
{
  std::vector<std::vector<std::vector<int>>> results;
  if (n == 1) {
    std::vector<std::vector<int>> simplex; 
    std::vector<int> v0, v1;
    v0.push_back(0);
    v1.push_back(1);
    simplex.push_back(v0);
    simplex.push_back(v1);
    results.push_back(simplex);
  } else {
    const auto results0 = subdivide_unit_cube(n-1);
    for (int i = 0; i < n; i ++) {
      for (int j = 0; j < results0.size(); j ++) {
        std::vector<std::vector<int>> simplex;
        for (int k = 0; k < n; k ++) {
          // std::string node = results0[j][k];
          std::vector<int> vertex = results0[j][k];
          // node.insert(i, "0");
          vertex.insert(vertex.begin()+i, 0);
          simplex.push_back(vertex);
        }
        std::vector<int> all_one_vertex(n);
        for (int k = 0; k < n; k ++)
          all_one_vertex[k] = 1;
        // simplex.push_back(std::string(n, '1'));
        simplex.push_back(all_one_vertex);
        results.push_back(simplex);
      }
    }
  }
  return results;
}

inline std::tuple<std::vector<std::vector<int>>, std::vector<int>>
simplicial_regular_mesh::reduce_unit_simplex(const std::vector<std::vector<int>> &simplex_) const
{
  auto simplex = simplex_;
  const size_t nvertices = simplex.size();
  std::vector<int> offset(nd(), 0);
  // std::string offset(nd(), '0');

  for (size_t i = 0; i < nd(); i ++) {
    // check if every coordinate on the i-th dimension is '1'
    bool all_one = true; 
    for (size_t j = 0; j < nvertices; j ++) {
      if (simplex[j][i] == 0) {
        all_one = false;
        break;
      }
    }

    // if every coordinate on the i-th dimension is '1', add the i-th
    // offset by 1 and set the i-th coordinate to '0'
    if (all_one) {
      offset[i] = 1;
      for (size_t j = 0; j < nvertices; j ++) 
        simplex[j][i] = 0;
    }
  }

  return std::make_tuple(simplex, offset);
}

inline std::vector<std::vector<std::vector<int>>> simplicial_regular_mesh::enumerate_unit_simplices(int n, int k)
{
  std::set<std::vector<std::vector<int>>> results;
  if (n == k) {
    auto subdivision = subdivide_unit_cube(n);
    for (auto simplex : subdivision) {
      std::sort(simplex.begin(), simplex.end());
      results.insert(simplex);
    }
  } else if (k < n) {
    auto hyper_simplices = enumerate_unit_simplices(n, k+1);
    for (auto hyper_simplex : hyper_simplices) {
      do {
        std::vector<std::vector<int>> simplex = hyper_simplex;
        simplex.resize(hyper_simplex.size() - 1);

        // auto [reduced_simplex, offset] = reduce_unit_simplex(simplex); // TODO: take care of offset
        auto tuple = reduce_unit_simplex(simplex); // TODO: take care of offset
        auto reduced_simplex = std::get<0>(tuple);
        auto offset = std::get<1>(tuple);
        std::sort(reduced_simplex.begin(), reduced_simplex.end());
        results.insert(reduced_simplex);
      } while (std::next_permutation(hyper_simplex.begin(), hyper_simplex.end()));
    }
  }
  
  std::vector<std::vector<std::vector<int>>> results1;
  for (const auto r : results)
    results1.push_back(r);
  return results1;
}

inline std::vector<std::tuple<int, std::vector<int>>> simplicial_regular_mesh::enumerate_unit_simplex_side_of(int k, int type)
{
  std::vector<std::tuple<int, std::vector<int>>> side_of;
  if (k == nd()) return side_of;

  std::set<std::tuple<int, std::vector<int>>> my_side_of;

  std::vector<int> zero_corner(nd(), 0);
  const auto k_simplex = simplicial_regular_mesh_element(zero_corner, k, type);
  const auto k_simplex_vertices = k_simplex.vertices(*this);

  int k_plus_one_simplex_type = 0;
  std::vector<int> corner(nd(), -1);
  while (1) {
    const simplicial_regular_mesh_element k_plus_one_simplex(corner, k+1, k_plus_one_simplex_type);
    const auto k_plus_one_simplex_vertices = k_plus_one_simplex.vertices(*this);

    // check if (k+1)-simplex contains the given k-simplex.
    if (std::includes(k_plus_one_simplex_vertices.begin(), k_plus_one_simplex_vertices.end(), 
          k_simplex_vertices.begin(), k_simplex_vertices.end())) 
      my_side_of.insert(std::make_tuple(k_plus_one_simplex_type, corner));

    if (k_plus_one_simplex_type == ntypes(k+1) - 1
        && std::all_of(corner.begin(), corner.end(), [](int i){return i == 1;})) break;
    else {
      if (k_plus_one_simplex_type == ntypes(k+1) - 1) {
        int i = 0;
        corner[i] ++;
        while (corner[i] > 1) {
          corner[i] = -1;
          corner[++i] ++;
        }
        k_plus_one_simplex_type = 0;
      } else k_plus_one_simplex_type ++;
    }
  };

  for (const auto i : my_side_of)
    side_of.push_back(i);

  return side_of;
}

inline std::vector<std::tuple<int, std::vector<int>>> simplicial_regular_mesh::enumerate_unit_simplex_sides(int k, int type)
{
  std::vector<std::tuple<int, std::vector<int>>> sides;
  if (k == 0) return sides;

  std::set<std::tuple<int, std::vector<int>>> mysides;

  auto simplex = unit_simplices[k][type]; 
  const auto k_minums_one_simplices = unit_simplices[k-1]; 
    
  do {
    std::vector<std::vector<int>> side = simplex;
    side.resize(simplex.size() - 1);

    // auto [reduced_side, offset] = reduce_unit_simplex(side);
    auto tuple = reduce_unit_simplex(side);
    auto reduced_side = std::get<0>(tuple);
    auto offset = std::get<1>(tuple);
    std::sort(reduced_side.begin(), reduced_side.end());

    int k_minums_one_type; // TODO: replace the linear search
    for (k_minums_one_type = 0; k_minums_one_type < k_minums_one_simplices.size(); k_minums_one_type ++) {
      if (reduced_side == k_minums_one_simplices[k_minums_one_type]) {
        // fprintf(stderr, "adding side, type=%d\n", k_minums_one_type);
        // sides.push_back(std::make_tuple(k_minums_one_type, offset));
        mysides.insert(std::make_tuple(k_minums_one_type, offset));
      }
    }

  } while (std::next_permutation(simplex.begin(), simplex.end())); // TODO: replace the permutation with combination

  for (const auto s : mysides)
    sides.push_back(s);

  // fprintf(stderr, "k=%d, type=%d, #sides=%lu\n", k, type, sides.size());
  // assert(sides.size() > 0);
  return sides;
}

inline void simplicial_regular_mesh::derive_ordinal_and_interval_simplices()
{
  for (int d = 0; d < nd()+1; d ++) {
    std::vector<int> ordinal_simplex_types, interval_simplex_types;
    std::vector<bool> is_type_ordinal;
    
    if (d == 0) {
      ordinal_simplex_types.push_back(0);
      is_type_ordinal.push_back(true);
    } else {
      for (int t = 0; t < ntypes(d); t ++) {
        const auto simplex = unit_simplices[d][t];
        int time = 0;
        for (int i = 0; i < simplex.size(); i ++) {
          time = time + simplex[i][nd()-1];
        }
        if (time == 0) {
          ordinal_simplex_types.push_back(t);
          is_type_ordinal.push_back(true);
        } else {
          interval_simplex_types.push_back(t);
          is_type_ordinal.push_back(false);
        }
      }
    }

    unit_ordinal_simplex_types.push_back(ordinal_simplex_types);
    unit_interval_simplex_types.push_back(interval_simplex_types);
    ntypes_ordinal_.push_back(ordinal_simplex_types.size());
    ntypes_interval_.push_back(interval_simplex_types.size());
    is_unit_simplex_type_ordinal.push_back(is_type_ordinal);
  }
}

#if 0
inline void simplicial_regular_mesh::partition(int np, std::vector<std::tuple<simplicial_regular_mesh, simplicial_regular_mesh>>& partitions) {
  std::vector<size_t> vector_zero = {0}; 
  partition(np, vector_zero, vector_zero, vector_zero, partitions); 
}

inline void simplicial_regular_mesh::partition(int np, const std::vector<size_t> &given, std::vector<std::tuple<simplicial_regular_mesh, simplicial_regular_mesh>>& partitions) {
  std::vector<size_t> vector_zero = {0}; 
  partition(np, given, vector_zero, vector_zero, partitions); 
}

inline void simplicial_regular_mesh::partition(int np, const std::vector<size_t> &given, const std::vector<size_t> &ghost, std::vector<std::tuple<simplicial_regular_mesh, simplicial_regular_mesh>>& partitions) {
  partition(np, given, ghost, ghost, partitions); 
}

inline void simplicial_regular_mesh::partition(int np, const std::vector<size_t> &given, const std::vector<size_t> &ghost_low, 
    const std::vector<size_t> &ghost_high, std::vector<std::tuple<simplicial_regular_mesh, simplicial_regular_mesh>>& partitions) {
  std::vector<std::tuple<lattice, lattice>> lattice_partitions;
  this->lattice_.partition(np, given, ghost_low, ghost_high, lattice_partitions); 

  // Lattice partitions 2 regular simplex mesh partitions
  for(auto& lattice_pair : lattice_partitions) {
    lattice& lattice_p = std::get<0>(lattice_pair); 
    lattice& lattice_ghost_p = std::get<1>(lattice_pair); 

    simplicial_regular_mesh p(this->nd(), lattice_p); 
    simplicial_regular_mesh ghost_p(this->nd(), lattice_ghost_p); 

    partitions.push_back(std::make_tuple(p, ghost_p)); 
  }
}
#endif

inline void simplicial_regular_mesh::print_unit_simplices(int nd, int d) const
{

}

inline void simplicial_regular_mesh::initialize_subdivision()
{
  ntypes_.resize(nd() + 1);
  unit_simplices.resize(nd() + 1);

  for (int k = 0; k <= nd(); k ++) {
    unit_simplices[k] = enumerate_unit_simplices(nd(), k);
    ntypes_[k] = unit_simplices[k].size();
#if 0
    for (const auto s : unit_simplices[k]) {
      for (const auto v : s)
        std::cerr << v << " ";
      std::cerr << std::endl;
    }
#endif
  }

  unit_simplex_sides.resize(nd()+1);
  for (int dim = 0; dim <= nd(); dim ++) {
    unit_simplex_sides[dim].resize(ntypes(dim));
    for (int type = 0; type < ntypes(dim); type ++) 
      unit_simplex_sides[dim][type] = enumerate_unit_simplex_sides(dim, type);
  }

  unit_simplex_side_of.resize(nd()+1);
  for (int dim = 0; dim <= nd(); dim ++) {
    unit_simplex_side_of[dim].resize(ntypes(dim));
    for (int type = 0; type < ntypes(dim); type ++) 
      unit_simplex_side_of[dim][type] = enumerate_unit_simplex_side_of(dim, type);
  }
  // enumerate_unit_simplex_side_of(0, 0);

  derive_ordinal_and_interval_simplices();
}

inline void simplicial_regular_mesh::set_lb_ub(const std::vector<int>& l, const std::vector<int>& u)
{
  lb_.resize(nd());
  ub_.resize(nd());
  dimprod_.resize(nd()+1);

  for (int i = 0; i < nd(); i ++) {
    lb_[i] = l[i];
    ub_[i] = u[i];
  }

  lattice_.reshape(lb_, sizes()); 

  for (int i = 0; i < nd()+1; i ++) {
    if (i == 0) dimprod_[i] = 1;
    else dimprod_[i] = (u[i-1] - l[i-1] + 1) * dimprod_[i-1];
  }
}

inline void simplicial_regular_mesh::set_lb_ub(const lattice& _lattice) {
  std::vector<size_t> _lattice_lb = _lattice.lower_bounds(); 
  std::vector<int> _lb(_lattice_lb.begin(), _lattice_lb.end());

  std::vector<size_t> _lattice_ub = _lattice.upper_bounds(); 
  std::vector<int> _ub(_lattice_ub.begin(), _lattice_ub.end()); 

  this->set_lb_ub(_lb, _ub);
}

inline std::vector<int> simplicial_regular_mesh::sizes() {
  std::vector<int> sizes; 
  for(int i = 0; i < nd(); ++i) {
    sizes.push_back(ub_[i] - lb_[i] + 1); 
  }

  return sizes; 
}

#if 0
inline simplicial_regular_mesh::iterator simplicial_regular_mesh::element_begin(int d, int scope)
{
  simplicial_regular_mesh_element e(*this, d);
  e.corner = lb();
  return e;
}

inline simplicial_regular_mesh::iterator simplicial_regular_mesh::element_end(int d, int scope)
{
  simplicial_regular_mesh_element e(*this, d);
  e.corner = ub();
  // e.type = ntypes(d) - 1;
  e.type = ntypes(d, scope); // the invalid type id combines with the ub corner encodes the end.
  return e;
}
#endif

inline size_t simplicial_regular_mesh::n(int d, int scope) const
{
  return (size_t)ntypes(d, scope) * dimprod_[nd()];
}
  
inline void simplicial_regular_mesh::element_for(int d, std::function<void(simplicial_regular_mesh_element)> f, int nthreads)
{
  element_for(d, lattice_, ELEMENT_SCOPE_ALL, f, nthreads);
}
  
inline void simplicial_regular_mesh::element_for_ordinal(int d, int t, std::function<void(simplicial_regular_mesh_element)> f, int nthreads)
{
  auto st = lattice_.starts(), sz = lattice_.sizes();
  st[nd()-1] = t;
  sz[nd()-1] = 1;

  lattice my_lattice(st, sz);
  // my_lattice.print(std::cerr);

  element_for(d, my_lattice, ELEMENT_SCOPE_ORDINAL, f, nthreads);
}
  
inline void simplicial_regular_mesh::element_for_interval(int d, int t0, int t1, std::function<void(simplicial_regular_mesh_element)> f, int nthreads)
{
  auto st = lattice_.starts(), sz = lattice_.sizes();
  st[nd()-1] = t0;
  sz[nd()-1] = 1; // TODO: t1

  lattice my_lattice(st, sz);

  element_for(d, my_lattice, ELEMENT_SCOPE_INTERVAL, f, nthreads);
}

inline void simplicial_regular_mesh::element_for(
    int d, const lattice& l, int scope, 
    std::function<void(simplicial_regular_mesh_element)> f,
    int nthreads)
{
  auto lambda = [=](size_t j) {
    simplicial_regular_mesh_element e(*this, d, j, l, scope);
    f(e);
  };

  const auto ntasks = l.n() * ntypes(d, scope);
  // fprintf(stderr,  "ntasks=%lu\n", ntasks);
#if FTK_HAVE_KOKKOS
  Kokkos::parallel_for("element_for", ntasks, KOKKOS_LAMBDA(const int& j) {f(j);});
#elif FTK_HAVE_TBB
  tbb::parallel_for(tbb::blocked_range<size_t>(0, ntasks),
      [=](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++ i) 
          lambda(i);
      });
#else
  std::vector<std::thread> workers;
  for (size_t i = 1; i < nthreads; i ++) {
    workers.push_back(std::thread([=]() {
      for (size_t j = i; j < ntasks; j += nthreads)
        lambda(j);
    }));
  }

  for (size_t j = 0; j < ntasks; j += nthreads) // the main thread
    lambda(j);

  std::for_each(workers.begin(), workers.end(), [](std::thread &t) {t.join();});
#endif
}

}


namespace diy {
  template <> struct Serialization<ftk::simplicial_regular_mesh_element> {
    static void save(diy::BinaryBuffer& bb, const ftk::simplicial_regular_mesh_element &e) {
      diy::save(bb, e.corner);
      diy::save(bb, e.dim);
      diy::save(bb, e.type);
    }

    static void load(diy::BinaryBuffer& bb, ftk::simplicial_regular_mesh_element& e) {
      diy::load(bb, e.corner); 
      diy::load(bb, e.dim);
      diy::load(bb, e.type);
    }
  };
}

#endif
