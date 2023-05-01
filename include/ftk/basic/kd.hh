#ifndef _FTK_KD_HH
#define _FTK_KD_HH

#include <ftk/config.hh>
#include <vector>
#include <list>

namespace ftk {

template <typename F, size_t n>
struct kd_node_t {
  size_t id; // median point id;
  size_t offset, size; // for building kd tree
  std::shared_ptr<kd_node_t> left, right;

  bool leaf() const { return left == nullptr && right == nullptr; }
};

template <typename F /*coordinate type, e.g., double*/, size_t n>
struct kd {
  void set_inputs(size_t npts, const F* p, size_t stride = n);
  void build(size_t max_level = std::numeric_limits<size_t>::max()); // no limit on max level
  std::shared_ptr<kd_node_t<F, n>> build_recursive(
      size_t max_level,
      size_t level,
      size_t offset, // 0 for the root
      size_t size,
      size_t ids[]);

  void print() const;
  void print_recursive(
      std::shared_ptr<kd_node_t<F, n>> node, 
      size_t depth) const;

  // size_t find_nearest(const F*) const;
  size_t find_nearest(const std::array<F, n> &x) const;
  size_t find_nearest_naive(const std::array<F, n> &x) const; // brute force
  std::shared_ptr<kd_node_t<F, n>> find_nearest_recursive(
      std::shared_ptr<kd_node_t<F, n>> node, 
      const std::array<F, n> &x, 
      size_t depth) const;

  std::shared_ptr<kd_node_t<F, n>> find_leaf_recursive(
      std::shared_ptr<kd_node_t<F, n>> node, 
      const std::array<F, n>& x, 
      size_t depth) const;
  std::shared_ptr<kd_node_t<F, n>> find_leaf(const std::array<F, n>& x) const { 
    return find_leaf_recursive(root, x, 0);
  }

  std::shared_ptr<kd_node_t<F, n>> closest(const std::array<F, n>& x, 
      std::shared_ptr<kd_node_t<F, n>> node1, 
      std::shared_ptr<kd_node_t<F, n>> node2) const;

  static F dist2(const std::array<F, n>& x, const std::array<F, n>& y);

  std::vector<std::array<F, n>> pts; // input points
  std::vector<size_t> ids;
  std::shared_ptr<kd_node_t<F, n>> root;
  
  // size_t parent(size_t i) const { return i / 2; }
  // size_t child(size_t i, unsigned char side) const { return i * 2 + side; }
};

////
template <typename F, size_t n>
void kd<F, n>::build(size_t max_level)
{
  ids.resize(pts.size());
  for (size_t i = 0; i < ids.size(); i ++) 
    ids[i] = i;

  root = build_recursive(max_level, 0, 0, pts.size(), &ids[0]);
}
  
template <typename F, size_t n>
F kd<F, n>::dist2(const std::array<F, n>& x, const std::array<F, n>& y)
{
  F dist(0);
  for (size_t k = 0; k < n; k ++) {
    const F d = x[k] - y[k];
    dist += d * d;
  }
  return dist;
}

template <typename F, size_t n>
void kd<F, n>::set_inputs(size_t npts, const F* p, size_t stride)
{
  pts.resize(npts);
  for (size_t i = 0; i < npts; i ++)
    for (size_t k = 0; k < n; k ++)
      pts[i][k] = p[i * stride + k];
}

template <typename F, size_t n>
std::shared_ptr<kd_node_t<F, n>> kd<F, n>::build_recursive(
    size_t max_level,
    size_t level,
    size_t offset, // 0 for the root
    size_t size,
    size_t ids[])
{
  std::shared_ptr<kd_node_t<F, n>> node(new kd_node_t<F, n>);
  // node->axis = axis;
  node->offset = offset;
  node->size = size;
  // node->derive_bounds(pts, ids);
    
  const size_t axis = level % n;

  // std::stable_sort(ids + offset, ids + offset + size, [axis, this](size_t i, size_t j) {
  //     return pts[i][axis] < pts[j][axis];
  // });
  std::nth_element(ids + offset, ids + offset + size / 2, ids + offset + size, 
      [axis, this](size_t i, size_t j) {
        return pts[i][axis] < pts[j][axis];
      });

  const size_t i = size / 2; // TODO: the median value may not be the split value
  node->id = ids[offset + i];
  // node->median = pts[ids[offset + i]][axis];

  if (level < max_level && size > 1) {
    node->left  = build_recursive(max_level, level + 1, offset, i, ids);
    node->right = build_recursive(max_level, level + 1, offset + i, size - i, ids);
  }

  return node;
}

template <typename F, size_t n>
void kd<F, n>::print() const
{
  print_recursive(root, 0);
}

template <typename F, size_t n>
void kd<F, n>::print_recursive(std::shared_ptr<kd_node_t<F, n>> node, size_t depth) const
{
  for (size_t i = 0; i < depth; i ++)
    printf("-");

  printf("id=%zu, pt=%f, %f, offset=%zu, size=%zu, leaf=%d, depth=%zu, axis=%zu\n", 
      node->id, pts[node->id][0], pts[node->id][1], 
      node->offset, node->size, node->leaf(), depth, depth % n);

  if (node->left) print_recursive(node->left, depth + 1);
  if (node->right) print_recursive(node->right, depth + 1);
}

template <typename F, size_t n>
size_t kd<F, n>::find_nearest(
    const std::array<F, n> &x) const 
{
  auto node = find_nearest_recursive(root, x, 0);
  return node->id;
}

template <typename F, size_t n>
size_t kd<F, n>::find_nearest_naive(
    const std::array<F, n> &x) const 
{
  size_t mi = 0;
  F mindist = std::numeric_limits<F>::max();
  for (size_t i = 0; i < pts.size(); i ++) {
    const F d2 = dist2(x, pts[i]);
    if (mindist > d2) {
      mindist = d2; 
      mi = i;
    }
  }
  return mi;
}
  
template <typename F, size_t n>
std::shared_ptr<kd_node_t<F, n>> kd<F, n>::closest(const std::array<F, n>& x, 
    std::shared_ptr<kd_node_t<F, n>> node1,
    std::shared_ptr<kd_node_t<F, n>> node2) const
{
  if (node1 == nullptr) return node2;
  else if (node2 == nullptr) return node1;
  else { 
    const F d1 = dist2(x, pts[node1->id]),
            d2 = dist2(x, pts[node2->id]);

    if (d1 < d2) return node1;
    else return node2;
  }
}
  
template <typename F, size_t n>
std::shared_ptr<kd_node_t<F, n>> kd<F, n>::find_nearest_recursive(
    std::shared_ptr<kd_node_t<F, n>> node, 
    const std::array<F, n> &x, 
    size_t depth) const
{
  if (!node) return nullptr;

  const size_t axis = depth % n;

  std::shared_ptr<kd_node_t<F, n>> next, other;
  if (x[axis] < pts[node->id][axis]) {
    next = node->left;
    other = node->right;
  } else {
    next = node->right;
    other = node->left;
  }
  
  auto temp = find_nearest_recursive(next, x, depth + 1);
  auto best = closest(x, temp, node);
  
  // fprintf(stderr, "depth=%zu, node=%p, temp=%p, best=%p\n",
  //     depth, node.get(), temp.get(), best.get());

  const F d2 = dist2(x, pts[best->id]);

  // check the other side
  const F dp = x[axis] - pts[node->id][axis];
  const F dp2 = dp * dp;

  if (d2 >= dp2) {
    temp = find_nearest_recursive(other, x, depth + 1);
    best = closest(x, temp, best);
  }

  return best;
}

template <typename F, size_t n>
std::shared_ptr<kd_node_t<F, n>> kd<F, n>::find_leaf_recursive(
    std::shared_ptr<kd_node_t<F, n>> node,
    const std::array<F, n> &x, 
    size_t depth) const
{
  if (node->leaf) 
    return node;
  else {
    if (x[depth % n] < pts[node->id][depth % n])
      return find_leaf_recursive(node->left, x, depth + 1);
    else 
      return find_leaf_recursive(node->right, x, depth + 1);
  }
}

} // namespace ftk

#endif
