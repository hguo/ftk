#ifndef _FTK_REGULAR_LATTICE_HH
#define _FTK_REGULAR_LATTICE_HH

#include <vector>
#include <ostream>

namespace ftk {

struct regular_lattice {
  regular_lattice() {}
  regular_lattice(int n);
  regular_lattice(const std::vector<size_t> &starts, const std::vector<size_t> &sizes) {reshape(starts, sizes);}
  regular_lattice(const std::vector<size_t> &sizes) {reshape(sizes);}

  void print(std::ostream& os);

  size_t nd() const {return sizes_.size();}
  size_t start(size_t i) {return starts_[i];}
  size_t size(size_t i) {return sizes_[i];}

  void reshape(const std::vector<size_t> &starts);
  void reshape(const std::vector<size_t> &starts, const std::vector<size_t> &sizes);

  size_t global_index(const std::vector<size_t> &coords) const;
  size_t local_index(int p, const std::vector<size_t> &coords) const;
  
public: // the last dimension, aka time.  these functions are mainly for I/O and streaming purposes
  bool unlimited_time() const {return unlimited_;}
  void set_unlimited_time(bool u) {unlimited_ = u;}
  
  void advance_time(int nt = 1) {starts_[nd()-1] += nt;}
  void recess_time(int nt = 1) {starts_[nd()-1] -= nt;}

public: // partitioning
  void partition(int np);
  size_t partition_id(const std::vector<size_t> &coords) const;

private:
  bool unlimited_ = false;
  std::vector<size_t> starts_, sizes_; // the last dimension can be unlimited
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
}

void regular_lattice::reshape(const std::vector<size_t> &sizes)
{
  sizes_ = sizes;
  starts_.resize(sizes.size(), 0);
}

}

#endif
