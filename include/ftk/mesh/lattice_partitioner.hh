#ifndef _FTK_LATTICE_PARTITIONER_HH
#define _FTK_LATTICE_PARTITIONER_HH

#include <ftk/mesh/lattice.hh>

namespace ftk {

struct lattice_partitioner {
  lattice_partitioner(const lattice& l_) : l(l_) {}
  
  friend std::ostream& operator<<(std::ostream& os, const lattice_partitioner&);

  size_t nd() const {return l.nd();}
  size_t np() const {return cores.size();}

public: // partitioning
  void partition(size_t np); // regular partition w/o ghosts
  void partition(size_t np, 
      const std::vector<size_t> &given); // regular partition with given number of cuts per dimension
  void partition(size_t np, 
      const std::vector<size_t> &given, 
      const std::vector<size_t> &ghost); // regular partition w/ given number of cuts and ghost size per dimension
  void partition(size_t np, 
      const std::vector<size_t> &given, 
      const std::vector<size_t> &ghost_low, 
      const std::vector<size_t> &ghost_high); // regular partition w/ given number of cuts and ghost sizes per dimension

  const lattice& get_core(size_t i) const {return cores[i];}
  const lattice& get_ext(size_t i) const {return extents[i];}

  size_t partition_id(const std::vector<size_t> &p) const;

private:
  void partition(std::vector<std::vector<size_t>> prime_factors_dims);

  template<typename uint = size_t>
  static std::vector<uint> prime_factorization(uint n);

  lattice add_ghost(const lattice& b, 
      const std::vector<size_t>&, 
      const std::vector<size_t>&) const;

protected:
  const lattice &l;

  std::vector<lattice> cores, extents; // core region and core+ghost regions
};

// factors of a given number n in decreasing order
template<typename uint>
inline std::vector<uint> lattice_partitioner::prime_factorization(uint n) {
  std::vector<uint> factors;

  while (n % 2 == 0) {  
    factors.push_back(2); 
    n = n/2; 
  }  

  for (int i = 3; i <= sqrt(n); i += 2) {
    while (n % i == 0) {
      factors.push_back(i); 
      n = n / i; 
    }
  }

  if(n > 1) {
    factors.push_back(n); 
  }

  std::reverse(factors.begin(), factors.end());
  return factors;
}

inline bool is_vector_zero(const std::vector<size_t> vec) {
  if(vec.size() == 0 || (vec.size() == 1 && vec[0] == 0)) {
    return true; 
  }

  return false; 
}

inline void lattice_partitioner::partition(size_t np) {
  std::vector<size_t> vector_zero(0, nd());
  partition(np, vector_zero, vector_zero, vector_zero);
}

inline void lattice_partitioner::partition(size_t np, const std::vector<size_t> &given) {
  std::vector<size_t> vector_zero(0, nd());
  partition(np, given, vector_zero, vector_zero); 
}

inline void lattice_partitioner::partition(size_t np, const std::vector<size_t> &given, const std::vector<size_t> &ghost)
{
  partition(np, given, ghost, ghost); 
}

inline void lattice_partitioner::partition(size_t np, 
    const std::vector<size_t> &given, 
    const std::vector<size_t> &ghost_low, 
    const std::vector<size_t> &ghost_high)
{
  if(np <= 0) {
    return ;
  }

  int ndim = l.nd_cuttable(); // # of cuttable dimensions

  // // Compute number of grid points by product
  // int ngridpoints = std::accumulate(sizes_.begin(), sizes_.end(), 1, std::multiplies<int>());  
  // if(ngridpoints <= np) {
  //   return empty; // return an empty vector 
  // }

  std::vector<std::vector<size_t>> prime_factors_dims(ndim, std::vector<size_t>()); // Prime factors of each dimension

  if(!is_vector_zero(given)) {
    bool has_zero = false; // Whither exist non-given dimensions
    // Compute np for non-given dimensions
    for(size_t nslice: given) {
      if(nslice != 0) {
        if(np % nslice != 0) { 
          return ;
        }

        np /= nslice; 
      } else {
        has_zero = true; 
      }
    }

    if(!has_zero && np > 1) { // If partitions of all dims are given, but np for non-given dim is not 1
      return ; 
    }

    for(auto i = 0; i < given.size(); ++i) {
      size_t nslice = given[i]; 
      if(nslice > 1) {
        std::vector<size_t> prime_factors_dim = prime_factorization(nslice);
        prime_factors_dims[i] = prime_factors_dim; 
      }
    }
  }
  
  // Get prime factors of non-given dimensions
  if(np > 1) {
    std::vector<size_t> prime_factors_all = prime_factorization(np);

    int curr = 0;
    for(size_t& prime : prime_factors_all) {
      while(!is_vector_zero(given) && given[curr] != 0) { // if current dim is given, then skip it
        curr = (curr + 1) % ndim;   
      }
      prime_factors_dims[curr].push_back(prime); 
      curr = (curr + 1) % ndim; 
    }
  }

  partition(prime_factors_dims); // get partitions given the prime factors of each dimension

  if (cores.size() == 0) return;

  // apply ghosts
  for(const auto& core : cores) {
    auto starts = core.starts(), 
         sizes = core.sizes();

    for(int d = 0; d < ndim; ++d) {
      // ghost_low layer
      {
        size_t offset = starts[d] - l.start(d); 
        if(ghost_low[d] < offset) {
          offset = ghost_low[d]; 
        }

        starts[d] -= offset; 
        sizes[d] += offset;
      }

      // ghost_high layer
      {
        size_t offset = (l.start(d) + l.size(d)) - (starts[d] + sizes[d]); 
        if(ghost_high[d] < offset) {
          offset = ghost_high[d]; 
        }

        sizes[d] += offset;
      }
    }

    lattice ext(starts, sizes);
    extents.push_back(ext);
  }
}


// Partition dimensions by given prime factors of each dimension
inline void lattice_partitioner::partition(std::vector<std::vector<size_t>> prime_factors_dims)
{
  int ndim = l.nd_cuttable(); // # of cuttable dimensions

  std::queue<lattice> lattice_queue; 
  lattice_queue.push(l);

  for(int curr = 0; curr < ndim; ++curr) { // current dim for cutting
    for(size_t& nslice : prime_factors_dims[curr]) {
      int n = lattice_queue.size();

      for (int j = 0; j < n; ++j) {
        auto p = lattice_queue.front();

        if(p.size(curr) < nslice) { // what if p.sizes(curr) < nslice? 
          return ; 
        }
        size_t ns = p.size(curr) / nslice; 

        std::vector<size_t> starts(p.starts()); 
        std::vector<size_t> sizes(p.sizes()); 
        sizes[curr] = ns;
        for(int k = 0; k < nslice - 1; ++k) {
          lattice_queue.push(lattice(starts, sizes));
          starts[curr] += ns;
        }

        sizes[curr] = p.size(curr) - ns * (nslice - 1);
        lattice_queue.push(lattice(starts, sizes));

        lattice_queue.pop();
      }
    }
  }

  // std::vector<std::tuple<lattice, lattice>> partitions;
  while (!lattice_queue.empty()) {

    lattice core = lattice_queue.front(); 
    // core.set_unlimited_time(l.unlimited_); //TODO

    lattice ghost(core.starts(), core.sizes()); 
    // ghost.set_unlimited_time(l.unlimited_); // TODO

    cores.push_back(core);
    // extents.push_back(ghost);
    
    lattice_queue.pop();
  }
}

// inline size_t lattice::global_index(const std::vector<size_t> &coords)
// {

// }

// // Each grid point is a regular lattice
// inline std::vector<lattice> lattice::partition_all() {
//   int ngridpoints = std::accumulate(sizes_.begin(), sizes_.end(), 1, std::multiplies<int>()); 

//   int ndim = unlimited_ ? nd() - 1 : nd(); // # of cuttable dimensions
  
//   std::vector<size_t> bounds(starts_);
//   for(int i = 0; i < ndim; ++i) {
//     bounds[i] += sizes_[i]; 
//   }

//   // Each time, the first dim plus one

//   std::vector<size_t> sizes(sizes_);
//   for(int i = 1; i < ndim; ++i) {
//     sizes[i] = 0; 
//   }
//   sizes[0] = 1;

//   std::vector<lattice> partitions;
//   std::vector<size_t> starts(starts_); 
//   for(int i = 0; i < ngridpoints; ++i) {
  
//     for(int j = 0; j < ndim; ++j) {
//       if(starts[j] < bounds[j]) {
//         starts[j] += 1;

//         break ;
//       } else {
//         starts[j] = 0; 
//       }
//     }
  
//     partitions.push_back(lattice(starts, sizes)); 
//   }

//   return partitions; 
// }
  
inline std::ostream& operator<<(std::ostream& os, const lattice_partitioner& partitioner)
{
  os << "lattice.nd=" << partitioner.nd() << ", np=" << partitioner.np() << std::endl;
  for (auto i = 0; i < partitioner.np(); i ++) {
    os << "--" << "id=" << i <<", "
               << "core: " << partitioner.cores[i] << "; "
               << "ext: " << partitioner.extents[i];
    if  (i < partitioner.np() - 1) os << std::endl;
  }
  
  return os;
}

};

#endif
