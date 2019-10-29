#ifndef _HYPERMESH_REGULAR_MESH_HH
#define _HYPERMESH_REGULAR_MESH_HH

#include <cmath>
#include <array>
#include <iostream>
#include <vector>

namespace ftk {

template <int N>
struct regular_mesh;

template <int N>
struct regular_mesh_element {
 public:
  friend class regular_mesh<N>;

  regular_mesh_element(const regular_mesh<N> &m_) : m(m_) { subspace.set(); }
  regular_mesh_element(const regular_mesh<N> &m_, const std::array<int, N>& c) :
      m(m_), coords(c) { subspace.set(); }
  regular_mesh_element(const regular_mesh<N> &m_, const std::bitset<N>& s) :
      m(m_), subspace(s) {}
  regular_mesh_element(
      const regular_mesh<N> &m_, const std::array<int, N>& c,
      const std::bitset<N>& s) : m(m_), coords(c), subspace(s) {}
  regular_mesh_element(const regular_mesh_element& element) :
      m(element.m) {
    coords = element.coords;
    subspace = element.subspace;
  }

  regular_mesh_element& operator++() {
    int d = m.dim(*this);

    auto it = std::find(m.element_types[d].begin(),
                        m.element_types[d].end(),
                        this->subspace);

    int element_types_i = std::distance(m.element_types[d].begin(), it);

    if (element_types_i < m.n_element_types(d) - 1) {
      this->subspace = m.element_types[d][element_types_i + 1];
    } else {
      int rez_i = 0;
      while (this->coords[rez_i] >= m.rez[rez_i] - 1) {
        this->coords[rez_i] = 0;
        rez_i ++;
      }
      this->coords[rez_i] ++;

      this->subspace = m.element_types[d][0];
    }

    return *this;
  }

  int dim() const {return subspace.count();}
  int n_sides() const {return 2 * dim();}
  // int n_nodes() const {return std::pow(2, dim());}

  std::vector<regular_mesh_element<N> > sides() const {
    std::vector<regular_mesh_element<N> > results;
    for (int i = 0; i < N; ++ i) {
      if (subspace[i]) {
        regular_mesh_element h0(*this);
        h0.subspace[i] = 0;
        results.push_back(h0);

        regular_mesh_element h1(*this);
        h1.subspace[i] = 0;
        h1.coords[i] ++;
        results.push_back(h1);
      }
    }
    return results;
  }

  std::vector<regular_mesh_element<N> > side_of() const {
    std::vector<regular_mesh_element<N> > results;
    for (int i = 0; i < N; ++ i) {
      if (!subspace[i]) {
        regular_mesh_element h0(*this);
        h0.subspace[i] = 1;
        results.push_back(h0);

        regular_mesh_element h1(*this);
        h1.subspace[i] = 1;
        h1.coords[i] --;
        results.push_back(h1);
      }
    }
    return results;
  }

  void print() const {
    std::cout << "coords: ";
    for (int i = N-1; i >= 0; -- i)
      std::cout << coords[i] << ", ";
    std::cout << "subspace: " << subspace << std::endl;
  }

 public:
  const regular_mesh<N> &m;

  std::array<int, N> coords = {};
  std::bitset<N> subspace;
};


template <int N>
struct regular_mesh {
 public:
  friend class regular_mesh_element<N>;
  typedef regular_mesh_element<N> iterator;

  regular_mesh() {
    initialize_element_types();
  }

  regular_mesh(const std::array<int, N>& r) : rez(r) {
    initialize_element_types();
  }

  iterator begin(int dim) {
    return iterator(*this, *(element_types[dim].begin()));
  }

  iterator end(int dim) {
    return iterator(*this, rezm1(), *(element_types[dim].end() - 1));
  }

  int n_element_types(int dim) const {
    if (dim < 0 || dim > N) return 0;
    else return element_types[dim].size();
  };

  std::array<int, N> rezm1() {
    std::array<int, N> result;
    for (int i = 0; i < N; ++i) {
      result[i] = rez[i] - 1;
    }
    return result;
  }

  std::vector<regular_mesh_element<N> > sides(const iterator& it) const {
    return it.sides();
  }

  std::vector<regular_mesh_element<N> > side_of(const iterator& it) const {
    return it.side_of();
  }

  int dim(const iterator& it) const {
    return it.subspace.count();
  }

 private:
  void initialize_element_types() {
    for (int k = 0; k < N+1; ++ k) {
      std::vector<bool> bitmask(k, 1); // K leading 1's
      bitmask.resize(N, 0); // N-K trailing 0's
      do {
        std::bitset<N> bitset;
        for (int i = 0; i < N; i ++)
          bitset[i] = bitmask[i];
        element_types[k].push_back(bitset);
      } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    }
  }

 public:
  std::array<int, N> rez;

  std::array<std::vector<std::bitset<N> >, N+1> element_types;
};

} // namespace hypermesh

#endif // _FTK_HYPERMESH_REGULAR_HH
