#ifndef _FTK_LATTICE_HH
#define _FTK_LATTICE_HH

#include <cmath>
#include <vector>
#include <tuple>
#include <queue>
#include <limits>
#include <ostream>
#include <limits>

namespace ftk {

struct lattice {
  friend class lattice_partitioner;

  lattice() {}
  lattice(int n);
  lattice(const std::vector<size_t> &starts, const std::vector<size_t> &sizes) {reshape(starts, sizes);}
  lattice(const std::vector<size_t> &sizes) {reshape(sizes);}

  friend std::ostream& operator<<(std::ostream& os, const lattice&);

  size_t nd() const {return sizes_.size();}
  size_t nd_cuttable() const {return unlimited_time() ? nd() - 1 : nd();}
  size_t start(size_t i) const {return starts_[i];}
  size_t size(size_t i) const {return sizes_[i];}
  size_t upper_bound(size_t i) const {return starts_[i] + sizes_[i] - 1;}

  size_t n() const {return prod_[nd()-1] * sizes_[nd()-1];}
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
  bool unlimited_time() const {return size(nd()-1) == std::numeric_limits<int>::max();}// return unlimited_;}
  // void set_unlimited_time(bool u) {unlimited_ = u;}

  // void advance_time(int nt = 1) {starts_[nd()-1] += nt;} // deprecated
  // void recess_time(int nt = 1) {starts_[nd()-1] -= nt;}

public:
  // bool unlimited_ = false; 
  std::vector<size_t> starts_, sizes_; // the last dimension can be unlimited
  std::vector<size_t> prod_; 
};

/////

inline lattice::lattice(int n)
{
  starts_.resize(n);
  sizes_.resize(n);
}

inline std::ostream& operator<<(std::ostream& os, const lattice& l)
{
  os << "starts={";
  for (int i = 0; i < l.nd(); i ++)
    if (i < l.nd()-1) os << l.starts_[i] << ",";
    else os << l.starts_[i] << "}, sizes={";

  for (int i = 0; i < l.nd(); i ++)
    if (i < l.nd()-1) os << l.sizes_[i] << ",";
    else os << l.sizes_[i] << "}, prod={";
  
  for (int i = 0; i < l.nd(); i ++)
    if (i < l.nd()-1) os << l.prod_[i] << ",";
    else os << l.prod_[i] << "}"; // << std::endl;

  return os;
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
      prod_[i] = prod_[i-1] * sizes[i-1];
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

}

#endif
