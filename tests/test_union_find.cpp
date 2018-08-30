#include <ftk/mesh/mesh.h>
#include <ftk/mesh/access.h>
#include <ftk/basic/union_find.h>
#include <functional>
#include <sstream>
#include <vector>
#include <queue>
#include <set>
#include <map>

template <class IdType>
struct Tree {
  void add_node(IdType i) { nodes.insert(i); }
  bool has_node(IdType i) const { return nodes.find(i) != nodes.end(); }

  void add_arc(IdType lo, IdType hi) {
    upper_links[lo].insert(hi);
    lower_links[hi].insert(lo);
  }
  
  bool reduce_node(IdType i) { // remove a node whose up-degree and down-degree are both less or equal to 1
    size_t ud = upper_degree(i), ld = lower_degree(i);
    if (ud > 1 || ld > 1) return false; // cannot do the node reduction
    else if (ud == 1 && ld == 0) {
      lower_links[upper_node(i)].erase(i);
    } else if (ud == 0 && ld == 1) {
      upper_links[lower_node(i)].erase(i);
    } else { // ud == 1 && ld == 1
      auto hi = upper_node(i), lo = lower_node(i);
      lower_links[hi].erase(i);
      upper_links[lo].erase(i);
      lower_links[hi].insert(lo);
      upper_links[lo].insert(hi);
    }
    nodes.erase(i);
    upper_links.erase(i);
    lower_links.erase(i);
    return true;
  }
 
  bool is_leaf(IdType i) const {
    return upper_degree(i) + lower_degree(i) == 1;
  }

  size_t upper_degree(IdType lo) const {
    auto it = upper_links.find(lo);
    if (it != upper_links.end()) return it->second.size();
    else return 0;
  }

  size_t lower_degree(IdType hi) const {
    auto it = lower_links.find(hi);
    if (it != lower_links.end()) return it->second.size();
    else return 0;
  }

  IdType upper_node(IdType lo) const { // return the first upper node
    auto it = upper_links.find(lo);
    if (it != upper_links.end()) return *it->second.begin();
    else return IdType(-1);
  }

  IdType lower_node(IdType hi) const { // return the first lower node
    auto it = lower_links.find(hi);
    if (it != lower_links.end()) return *it->second.begin();
    else return IdType(-1);
  }

  const std::set<IdType>& upper_nodes(IdType i) const {
    auto it = upper_links.find(i);
    if (it != upper_links.end()) return it.second;
    else return std::set<IdType>();
  }

  const std::set<IdType>& lower_nodes(IdType i) const {
    auto it = lower_links.find(i);
    if (it != lower_links.end()) return it.second;
    else return std::set<IdType>();
  }

  template <class ValueType>
  void print_with_values(const std::function<ValueType(IdType)>& value) const {
    for (const auto &kv : upper_links) {
      IdType i = kv.first;
      for (const auto j : kv.second) {
        fprintf(stderr, "%zu(%f) <--> %zu(%f)\n", 
            i, value(i), j, value(j));
      }
    }
  }

protected:
  std::set<IdType> nodes;
  std::map<IdType, std::set<IdType> > upper_links, lower_links;
};

template <class IdType>
Tree<IdType> merge_trees(IdType nn, Tree<IdType>& jt, Tree<IdType>& st)
{
  Tree<IdType> ct;
  std::queue<IdType> Q;

  for (IdType i = 0; i < nn; i ++) {
    ct.add_node(i);
    if (jt.upper_degree(i) + st.lower_degree(i) == 1) {
      Q.push(i);
    }
  }

  while (!Q.empty()) {
    auto i = Q.front();
    Q.pop();

    size_t j;
    if (jt.upper_degree(i) == 0 && st.lower_degree(i) == 1) { // upper leaf in split tree
      j = jt.lower_node(i);
      ct.add_arc(j, i);
    } else if (st.lower_degree(i) == 0 && jt.upper_degree(i) == 1) { // lower leaf in join tree
      j = st.upper_node(i);
      ct.add_arc(i, j);
    } else {
      continue; // i is no longer a leaf.
    }
      
    jt.reduce_node(i);
    st.reduce_node(i);

    if (jt.upper_degree(j) + st.lower_degree(j) == 1) { // new leaf node
      Q.push(j);
    }
  }

  return ct;
}

template <class IdType>
Tree<IdType> build_join_tree(IdType nn, 
    const std::vector<IdType> &order,
    const std::vector<IdType> &inverseOrder,
    const std::function<std::set<size_t>(size_t)> &neighbors)
{
  ftk::quick_union<IdType> uf(nn);
  Tree<IdType> jt;

  for (int i=nn-1; i>=0; i--) {
    jt.add_node(order[i]);
    for (auto j : neighbors(order[i])) {
      j = inverseOrder[j];
      if (j > i) {
        size_t ri = uf.root(i), rj = uf.root(j);
        if (ri != rj) {
          jt.add_arc(order[ri], order[rj]);
          uf.unite(rj, ri);
        }
      }
    }
  }

  return jt;
}

template <class IdType>
Tree<IdType> build_split_tree(IdType nn, 
    const std::vector<IdType> &order,
    const std::vector<IdType> &inverseOrder,
    const std::function<std::set<size_t>(size_t)> &neighbors)
{
  ftk::quick_union<IdType> uf(nn);
  Tree<IdType> st;

  for (int i=0; i<nn; i++) {
    st.add_node(order[i]);
    for (auto j : neighbors(order[i])) {
      j = inverseOrder[j];
      if (j < i) {
        size_t ri = uf.root(i), rj = uf.root(j);
        if (ri != rj) {
          st.add_arc(order[rj], order[ri]);
          uf.unite(rj, ri);
        }
      }
    }
  }

  return st;
}

template <class IdType>
Tree<IdType> sweep_and_merge(IdType nn,
    const std::vector<IdType> &order,
    const std::vector<IdType> &inverseOrder,
    const std::function<std::set<size_t>(size_t)> &neighbors)
{
  Tree<IdType> st = build_split_tree(nn, order, inverseOrder, neighbors), 
               jt = build_join_tree(nn, order, inverseOrder, neighbors);
  
  return merge_trees<IdType>(nn, jt, st);
}

int main(int argc, char **argv)
{
#if 1
  const int W = 5, H = 3;
  std::vector<double> values = {
    0, 15, 4, 7, 2, 
    17, 20, 5, 10, 8, 
    13, 16, 1, 6, 3};
#endif
#if 0
  const int W = 3, H = 3;
  std::vector<double> values = {
    8, 7, 4, 
    6, 5, 2, 
    3, 1, 0};
#endif

  std::vector<size_t> order(W*H), inverseOrder(W*H);
  for (size_t i=0; i<W*H; i++) 
    order[i] = i;

  std::stable_sort(order.begin(), order.end(), 
      [&values](size_t p, size_t q) {
        return values[p] < values[q];
      });

  for (size_t i=0; i<W*H; i++)
    inverseOrder[order[i]] = i;

  auto ct = sweep_and_merge<size_t>(W*H, order, inverseOrder,
      std::bind(ftk::Get6Neighbors2DRegular<size_t>, W, H, std::placeholders::_1));

  ct.print_with_values<double>([&values](size_t i) {return values[i];});

  return 0;
}
