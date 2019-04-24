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

std::vector<regular_lattice> regular_lattice::partition(int np) {
  std::vector<regular_lattice> partitions;

  int ncut = std::log2(np);
  int cuttable = unlimited_ ? nd() - 1 : nd();

  std::queue<regular_lattice> tmp;
  tmp.push(*this);
  int curr = 0;

  for (int i = 0; i < ncut; ++i) {
    int n = tmp.size();
    for (int j = 0; j < n; ++j) {
      auto p = tmp.front();
      size_t ns = p.size(curr) >> 1;

      std::vector<size_t> starts(p.starts()), sizes(p.sizes());
      sizes[curr] = ns;

      tmp.push(regular_lattice(starts, sizes));

      starts[curr] += ns;
      sizes[curr] = p.size(curr) - ns;

      tmp.push(regular_lattice(starts, sizes));

      tmp.pop();
    }

    curr = curr + 1 == cuttable ? 0 : curr + 1;
  }

  while (!tmp.empty()) {
    partitions.push_back(tmp.front());
    tmp.pop();
  }

  return partitions;
}

}

#endif
