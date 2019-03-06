#ifndef _HYPERMESH_REGULAR_MESH_HH
#define _HYPERMESH_REGULAR_MESH_HH

#include <cmath>
#include <array>
#include <iostream>
#include <vector>

namespace hypermesh {

template <int N>
struct regular_element {
 public:
  regular_element() {subspace.set();}
  explicit regular_element(const std::array<int, N>& c) :
      coords(c) {subspace.set();}
  explicit regular_element(const std::bitset<N>& s) :
      subspace(s) {}
  regular_element(const std::array<int, N>& c, const std::bitset<N>& s) :
      coords(c), subspace(s) {}
  regular_element(const regular_element& element) {
    coords = element.coords;
    subspace = element.subspace;
  }

  int dim() const {return subspace.count();}
  int n_sides() const {return 2 * dim();}
  // int n_nodes() const {return std::pow(2, dim());}

  std::vector<regular_element<N> > sides() const {
    std::vector<regular_element<N> > results;
    for (int i = 0; i < N; ++ i) {
      if (subspace[i]) {
        regular_element h0(*this);
        h0.subspace[i] = 0;
        results.push_back(h0);

        regular_element h1(*this);
        h1.subspace[i] = 0;
        h1.coords[i] ++;
        results.push_back(h1);
      }
    }
    return results;
  }

  std::vector<regular_element<N> > side_of() const {
    std::vector<regular_element<N> > results;
    for (int i = 0; i < N; ++ i) {
      if (!subspace[i]) {
        regular_element h0(*this);
        h0.subspace[i] = 1;
        results.push_back(h0);

        regular_element h1(*this);
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
  std::array<int, N> coords;
  std::bitset<N> subspace;
};


template <int N>
struct regular {
 public:
  regular() {
    initialize_element_types();
  }

  regular(const std::array<int, N>& r) : rez(r) {
    initialize_element_types();
  }

  // iterator
  struct iterator {
    iterator(regular<N>& r) : regular(r), element() {}
    explicit iterator(regular<N>& r, const std::array<int, N>& c) :
        regular(r), element(c) {}
    explicit iterator(regular<N>& r, const std::bitset<N>& s) :
        regular(r), element(s) {}
    iterator(regular<N>& r, const std::array<int, N>& c,
             const std::bitset<N>& s) :
        regular(r), element(c, s) {}

    iterator& operator++() {
      int d = regular.dim(*this);

      auto it = std::find(regular.element_types[d].begin(),
                          regular.element_types[d].end(),
                          this->element.subspace);

      int element_types_i = std::distance(regular.element_types[d].begin(), it);

      if (element_types_i < regular.n_element_types(d) - 1) {
        this->element.subspace = regular.element_types[d][element_types_i + 1];
      } else {
        int rez_i = 0;
        while (this->element.coords[rez_i] >= regular.rez[rez_i] - 1) {
          this->element.coords[rez_i] = 0;
          rez_i ++;
        }
        this->element.coords[rez_i] ++;

        this->element.subspace = regular.element_types[d][0];
      }

      return *this;
    }

    regular<N>& regular;
    regular_element<N> element;
  };

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

  std::vector<regular_element<N> > sides(const iterator& it) const {
    return it.element.sides();
  }

  std::vector<regular_element<N> > side_of(const iterator& it) const {
    return it.element.side_of();
  }

  int dim(const iterator& it) const {
    return it.element.subspace.count();
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
