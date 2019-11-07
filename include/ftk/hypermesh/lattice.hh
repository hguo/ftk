#ifndef _FTK_LATTICE_HH
#define _FTK_LATTICE_HH

#include <cmath>
#include <vector>
#include <tuple>
#include <queue>
#include <ostream>

namespace ftk {

struct lattice {
  lattice() {}
  lattice(int n);
  lattice(const std::vector<size_t> &starts, const std::vector<size_t> &sizes) {reshape(starts, sizes);}
  lattice(const std::vector<size_t> &sizes) {reshape(sizes);}

  void print(std::ostream& os);

  size_t nd() const {return sizes_.size();}
  size_t nd_cuttable() const {return unlimited_ ? nd() - 1 : nd();}
  size_t start(size_t i) const {return starts_[i];}
  size_t size(size_t i) const {return sizes_[i];}
  size_t upper_bound(size_t i) const {return starts_[i] + sizes_[i] - 1;}

  const std::vector<size_t>& starts() const {return starts_;}
  const std::vector<size_t>& sizes() const {return sizes_;}

  std::vector<size_t> lower_bounds() const;
  std::vector<size_t> upper_bounds() const; 

  void reshape(const std::vector<size_t> &sizes);
  void reshape(const std::vector<int> &sizes);
  void reshape(const std::vector<size_t> &starts, const std::vector<size_t> &sizes);
  void reshape(const std::vector<int> &starts, const std::vector<int> &sizes);

  template <typename uint=uint64_t> uint to_integer(const std::vector<int> &coords) const;
  template <typename uint=uint64_t> std::vector<int> from_integer(uint i) const;

  // size_t global_index(const std::vector<size_t> &coords) const;
  size_t local_index(int p, const std::vector<size_t> &coords) const;

public: // the last dimension, aka time.  these functions are mainly for I/O and streaming purposes
  bool unlimited_time() const {return unlimited_;}
  void set_unlimited_time(bool u) {unlimited_ = u;}

  void advance_time(int nt = 1) {starts_[nd()-1] += nt;}
  void recess_time(int nt = 1) {starts_[nd()-1] -= nt;}


// partitioning
public: 
  void partition(int np, std::vector<std::tuple<lattice, lattice>>& partitions) const;
  void partition(int np, const std::vector<size_t> &given, std::vector<std::tuple<lattice, lattice>>& partitions) const;
  void partition(int np, const std::vector<size_t> &given, const std::vector<size_t> &ghost, std::vector<std::tuple<lattice, lattice>>& partitions) const;
  void partition(int np, const std::vector<size_t> &given, const std::vector<size_t> &ghost_low, const std::vector<size_t> &ghost_high, std::vector<std::tuple<lattice, lattice>>& partitions) const;
private:
  void partition(std::vector<std::vector<size_t>> prime_factors_dims, std::vector<std::tuple<lattice, lattice>>& partitions) const; 

  // std::vector<lattice> partition_all();
  size_t partition_id(const std::vector<size_t> &coords) const;

private:
  bool unlimited_ = false; 
  std::vector<size_t> starts_, sizes_; // the last dimension can be unlimited
  std::vector<size_t> prod_; 
};

/////

inline lattice::lattice(int n)
{
  starts_.resize(n);
  sizes_.resize(n);
}

inline void lattice::print(std::ostream& os)
{
  os << "starts: {";
  for (int i = 0; i < nd(); i ++)
    if (i < nd()-1) os << starts_[i] << ",";
    else os << starts_[i] << "}, sizes: {";

  for (int i = 0; i < nd(); i ++)
    if (i < nd()-1) os << sizes_[i] << ",";
    else os << sizes_[i] << "}" << std::endl;
}


inline std::vector<size_t> lattice::lower_bounds() const {
  std::vector<size_t> lb; 
  for(int i = 0; i < nd(); ++i) {
    lb.push_back(starts_[i]); 
  }

  return lb; 
}

inline std::vector<size_t> lattice::upper_bounds() const {
  std::vector<size_t> ub; 
  for(int i = 0; i < nd(); ++i) {
    ub.push_back(starts_[i] + sizes_[i] - 1); 
  }

  return ub; 
}

inline void lattice::reshape(const std::vector<size_t> &starts, const std::vector<size_t> &sizes)
{
  starts_ = starts;
  sizes_ = sizes;
  prod_.resize(sizes.size());

  for (int i = 0; i < nd(); i ++)
    if (i == 0)
      prod_[i] = 1;
    else
      prod_[i] = prod_[i-1] * sizes[i];
}

inline void lattice::reshape(const std::vector<int> &starts, const std::vector<int> &sizes) 
{
  // Convert int to size_t

  std::vector<size_t> _starts(starts.begin(), starts.end()); 
  std::vector<size_t> _sizes(sizes.begin(), sizes.end()); 

  reshape(_starts, _sizes); 
}

inline void lattice::reshape(const std::vector<size_t> &sizes)
{
  starts_.resize(sizes.size());
  reshape(starts_, sizes);
}

inline void lattice::reshape(const std::vector<int> &sizes)
{
  // Convert int to size_t
  std::vector<size_t> _sizes(sizes.size()); 
  for(int i = 0; i < sizes.size(); ++i) {
    _sizes[i] = sizes[i]; 
  } 

  reshape(_sizes); 
}

template <typename uint> 
inline uint lattice::to_integer(const std::vector<int> &idx1) const
{
  std::vector<int> idx(idx1);
  for (auto j = 0; j < nd(); j ++)
    idx[j] -= starts_[j];

  uint i(idx[0]);
  for (auto j = 1; j < nd(); j ++)
    i += idx[j] * prod_[j];
  return i;
}

template <typename uint>
inline std::vector<int> lattice::from_integer(uint i) const
{
  std::vector<int> idx(nd());
  for (auto j = nd()-1; j > 0; j --) {
    idx[j] = i / prod_[j];
    i -= idx[j] * prod_[j];
  }
  idx[0] = i;

  for (auto j = 0; j < nd(); j ++)
    idx[j] += starts_[j];

  return idx;
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

// factors of a given number n in decreasing order
inline void prime_factorization(int n, std::vector<size_t>& factors) {  
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
}

inline bool is_vector_zero(const std::vector<size_t> vec) {
  if(vec.size() == 0 || (vec.size() == 1 && vec[0] == 0)) {
    return true; 
  }

  return false; 
}

inline void lattice::partition(int np, std::vector<std::tuple<lattice, lattice>>& partitions) const {  
  std::vector<size_t> vector_zero = {0}; 
  partition(np, vector_zero, vector_zero, vector_zero, partitions); 
}

inline void lattice::partition(int np, const std::vector<size_t> &given, std::vector<std::tuple<lattice, lattice>>& partitions) const {
  std::vector<size_t> vector_zero = {0}; 
  partition(np, given, vector_zero, vector_zero, partitions); 
}

inline void lattice::partition(int np, const std::vector<size_t> &given, const std::vector<size_t> &ghost, std::vector<std::tuple<lattice, lattice>>& partitions) const {
  partition(np, given, ghost, ghost, partitions); 
}

inline void lattice::partition(int np, const std::vector<size_t> &given, const std::vector<size_t> &ghost_low, const std::vector<size_t> &ghost_high, std::vector<std::tuple<lattice, lattice>>& partitions) const 
{
  if(np <= 0) {
    return ;
  }

  int ndim = this->nd_cuttable(); // # of cuttable dimensions

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

    for(int i = 0; i < given.size(); ++i) {
      int nslice = given[i]; 
      if(nslice > 1) {
        std::vector<size_t> prime_factors_dim;
        prime_factorization(nslice, prime_factors_dim); 

        prime_factors_dims[i] = prime_factors_dim; 
      }
    }
  }
  
  // Get prime factors of non-given dimensions
  if(np > 1) {
    std::vector<size_t> prime_factors_all;
    prime_factorization(np, prime_factors_all); 

    int curr = 0;
    for(size_t& prime : prime_factors_all) {
      while(!is_vector_zero(given) && given[curr] != 0) { // if current dim is given, then skip it
        curr = (curr + 1) % ndim;   
      }
      prime_factors_dims[curr].push_back(prime); 
      curr = (curr + 1) % ndim; 
    }
  }

  partition(prime_factors_dims, partitions); // get partitions given the prime factors of each dimension

  if(partitions.size() == 0) {
    return ;
  }

  // Add ghost_low layer
  if(!is_vector_zero(ghost_low)) {
    for(auto& p : partitions) {
      for(int d = 0; d < ndim; ++d) {
        lattice& ghost_p = std::get<1>(p); 

        std::vector<size_t> starts(ghost_p.starts());  
        std::vector<size_t> sizes(ghost_p.sizes()); 

        size_t offset = starts[d] - this->start(d); 
        if(ghost_low[d] < offset) {
          offset = ghost_low[d]; 
        }

        starts[d] -= offset; 
        sizes[d] += offset; 

        ghost_p.reshape(starts, sizes); 
      }
    }
  }

  // Add ghost_high layer
  if(!is_vector_zero(ghost_high)) {
    for(auto& p : partitions) {
      for(int d = 0; d < ndim; ++d) {
        lattice& ghost_p = std::get<1>(p); 

        std::vector<size_t> sizes(ghost_p.sizes()); 

        size_t offset = (this->start(d) + this->size(d)) - (ghost_p.start(d) + sizes[d]); 
        if(ghost_high[d] < offset) {
          offset = ghost_high[d]; 
        }

        sizes[d] += offset;      

        ghost_p.reshape(ghost_p.starts(), sizes); 
      }
    }
  }
}


// Partition dimensions by given prime factors of each dimension
inline void lattice::partition(std::vector<std::vector<size_t>> prime_factors_dims, std::vector<std::tuple<lattice, lattice>>& partitions) const {

  int ndim = this->nd_cuttable(); // # of cuttable dimensions

  std::queue<lattice> lattice_queue; 
  lattice_queue.push(*this);

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
    core.set_unlimited_time(this->unlimited_); 

    lattice ghost(core.starts(), core.sizes()); 
    ghost.set_unlimited_time(this->unlimited_);     

    std::tuple<lattice, lattice> p = std::make_tuple(core, ghost); 
    
    partitions.push_back(p); 
    lattice_queue.pop();
  }
}

}

#endif
