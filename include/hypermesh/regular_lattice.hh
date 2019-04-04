#ifndef _FTK_REGULAR_LATTICE_HH
#define _FTK_REGULAR_LATTICE_HH

#include <vector>

namespace ftk {

struct regular_lattice {
  regular_lattice() {}
  regular_lattice(int n);
  regular_lattice(const std::vector<size_t> &starts, const std::vector<size_t> &sizes) {reshape(starts, sizes);}
  regular_lattice(const std::vector<size_t> &sizes) {reshape(sizes);}

  size_t nd() const {return sizes.size();}
  size_t start(size_t i) {return starts[i];}
  size_t size(size_t i) {return sizes[i];}

  void reshape(const std::vector<size_t> &starts);
  void reshape(const std::vector<size_t> &starts, const std::vector<size_t> &sizes);

  size_t global_index(const std::vector<size_t> &coords) const;
  size_t local_index(int p, const std::vector<size_t> &coords) const;

public: // partitioning
  void partition(int np);

private:
  std::vector<size_t> starts, sizes; // the last dimension can be unlimited
  std::vector<regular_lattice> partitions
};

/////

void regular_lattice::reshape(const std::vector<size_t> &starts_, const std::vector<size_t> &sizes_)
{
  starts = starts_;
  sizes = sizes_;
}

void regular_lattice::reshape(const std::vector<size_t> &sizes_)
{
  sizes = sizes_;
  starts.resize(sizes.size(), 0);
}

void regular_lattice::partition(int np)
{
  int rem = np;

  while (1) {

  }
}

}

#endif
