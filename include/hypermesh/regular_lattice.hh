#ifndef _FTK_REGULAR_LATTICE_HH
#define _FTK_REGULAR_LATTICE_HH

#include <cmath>

#include <vector>
#include <queue>
#include <ostream>

namespace hypermesh {

struct regular_lattice {
  regular_lattice() {}
  regular_lattice(int n);
  regular_lattice(const std::vector<size_t> &starts, const std::vector<size_t> &sizes) {reshape(starts, sizes);}
  regular_lattice(const std::vector<size_t> &sizes) {reshape(sizes);}

  void print(std::ostream& os);

  size_t nd() const {return sizes_.size();}
  size_t start(size_t i) {return starts_[i];}
  size_t size(size_t i) {return sizes_[i];}

  std::vector<size_t>& starts() {return starts_;}
  std::vector<size_t>& sizes() {return sizes_;}

  const std::vector<size_t>& starts() const {return starts_;}
  const std::vector<size_t>& sizes() const {return sizes_;}

  void reshape(const std::vector<size_t> &sizes);
  void reshape(const std::vector<size_t> &starts, const std::vector<size_t> &sizes);

  size_t global_index(const std::vector<size_t> &coords) const;
  size_t local_index(int p, const std::vector<size_t> &coords) const;

public: // the last dimension, aka time.  these functions are mainly for I/O and streaming purposes
  bool unlimited_time() const {return unlimited_;}
  void set_unlimited_time(bool u) {unlimited_ = u;}

  void advance_time(int nt = 1) {starts_[nd()-1] += nt;}
  void recess_time(int nt = 1) {starts_[nd()-1] -= nt;}

public: // partitioning
  std::vector<regular_lattice> partition(int np);
  std::vector<regular_lattice> partition_all();
  size_t partition_id(const std::vector<size_t> &coords) const;

private:
  bool unlimited_ = false;
  std::vector<size_t> starts_, sizes_; // the last dimension can be unlimited
  std::vector<size_t> local_prod_, global_prod_;
};

/////

regular_lattice::regular_lattice(int n)
{
  starts_.resize(n);
  sizes_.resize(n);
}

void regular_lattice::print(std::ostream& os)
{
  os << "starts: {";
  for (int i = 0; i < nd(); i ++)
    if (i < nd()-1) os << starts_[i] << ",";
    else os << starts_[i] << "}, sizes: {";

  for (int i = 0; i < nd(); i ++)
    if (i < nd()-1) os << sizes_[i] << ",";
    else os << sizes_[i] << "}" << std::endl;
}

void regular_lattice::reshape(const std::vector<size_t> &starts, const std::vector<size_t> &sizes)
{
  starts_ = starts;
  sizes_ = sizes;

  local_prod_.resize(sizes.size());
  global_prod_.resize(sizes.size());

  for (int i = 0; i < nd(); i ++) {
    if (i == 0) {
      local_prod_[i] = 1;
      global_prod_[i] = 1;
    } else {
      local_prod_[i] = local_prod_[i-1] * sizes[i];
      global_prod_[i] = global_prod_[i-1] * (starts[i] + sizes[i]);
    }
  }
}

void regular_lattice::reshape(const std::vector<size_t> &sizes)
{
  starts_.resize(sizes.size());
  reshape(starts_, sizes);
}

// size_t regular_lattice::global_index(const std::vector<size_t> &coords)
// {

// }


// factors of a given number n in decreasing order
void prime_factorization(int n, std::vector<int>& factors) {  
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

  std::reverse(factors.begin(), factors.end()); 
}  

std::vector<regular_lattice> regular_lattice::partition_all() {
  int ngridpoints = std::accumulate(sizes_.begin(), sizes_.end(), 1, std::multiplies<int>()); 

  int ndim = unlimited_ ? nd() - 1 : nd(); // # of cuttable dimensions
  
  std::vector<size_t> bounds(starts_);
  for(int i = 0; i < ndim; ++i) {
    bounds[i] += sizes_[i]; 
  }

  // Each time, the first dim plus one

  std::vector<size_t> sizes(sizes_);
  for(int i = 1; i < ndim; ++i) {
    sizes[i] = 0; 
  }
  sizes[0] = 1;

  std::vector<regular_lattice> partitions;
  std::vector<size_t> starts(starts_); 
  for(int i = 0; i < ngridpoints; ++i) {
  
    for(int j = 0; j < ndim; ++j) {
      if(starts[j] < bounds[j]) {
        starts[j] += 1;

        break ;
      } else {
        starts[j] = 0; 
      }
    }
  
    partitions.push_back(regular_lattice(starts, sizes)); 
  }

  return partitions; 
}


std::vector<regular_lattice> regular_lattice::partition(int np) {
  // Compute number of grid points by product
  int ngridpoints = std::accumulate(sizes_.begin(), sizes_.end(), 1, std::multiplies<int>());  
  if(ngridpoints <= np) {
    return partition_all(); 
  }

  std::vector<int> prime_factors;
  prime_factorization(np, prime_factors); 

  // int ncut = std::log2(np);
  // int ncut = prime_factors.size(); 
  int ndim = unlimited_ ? nd() - 1 : nd(); // # of cuttable dimensions

  std::queue<regular_lattice> lattice_queue; 
  lattice_queue.push(*this);
  int curr = 0; // current dim for cutting

  // for (int i = 0; i < ncut; ++i) {

  for(int& nslice : prime_factors): 
    int n = lattice_queue.size();

    for (int j = 0; j < n; ++j) {
      auto p = lattice_queue.front();

      size_t ns = p.size(curr) / nslice; 

      std::vector<size_t> starts(p.starts()); 
      std::vector<size_t> sizes(p.sizes()); 
      sizes[curr] = ns;
      for(int k = 0; k < nslice - 1; ++k) {
        lattice_queue.push(regular_lattice(starts, sizes));
        starts[curr] += ns;
      }

      sizes[curr] = p.size(curr) - ns * nslice;
      lattice_queue.push(regular_lattice(starts, sizes));

      lattice_queue.pop();
    }

    curr = (curr + 1) % ndim; 
  }

  std::vector<regular_lattice> partitions;
  while (!lattice_queue.empty()) {
    partitions.emplace_back(std::move(lattice_queue.front()));
    lattice_queue.pop();
  }

  return partitions;
}

}

#endif
