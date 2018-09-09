#ifndef _FTK_SWEEP_AND_MERGE_H
#define _FTK_SWEEP_AND_MERGE_H

#include <queue>
#include "ftk/basic/union_find.h"
#include "ftk/basic/contour_tree.h"

namespace ftk {

template <class IdType, class ValueType>
contour_tree<IdType> build_contour_tree(IdType nn,
    const std::vector<ValueType> &vector,
    const std::function<std::set<IdType>(IdType)> &neighbors);

template <class IdType, class ValueType>
contour_tree<IdType> build_contour_tree(IdType nn,
    const std::function<ValueType(IdType)> &value,
    const std::function<std::set<IdType>(IdType)> &neighbors);

// internal functions
template <class IdType>
contour_tree<IdType> merge_join_and_split_trees(IdType nn, contour_tree<IdType>& jt, contour_tree<IdType>& st);

template <class IdType>
contour_tree<IdType> build_join_tree(IdType nn, 
    const std::vector<IdType> &order,
    const std::vector<IdType> &inverse_order,
    const std::function<std::set<IdType>(IdType)> &neighbors);

template <class IdType>
contour_tree<IdType> build_split_tree(IdType nn, 
    const std::vector<IdType> &order,
    const std::vector<IdType> &inverse_order,
    const std::function<std::set<IdType>(IdType)> &neighbors);
}


////////////////////////////
// implementations
namespace ftk {

template <class IdType>
contour_tree<IdType> merge_join_and_split_trees(IdType nn, contour_tree<IdType>& jt, contour_tree<IdType>& st)
{
  contour_tree<IdType> ct;
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

    IdType j;
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
contour_tree<IdType> build_join_tree(IdType nn, 
    const std::vector<IdType> &order,
    const std::vector<IdType> &inverse_order,
    const std::function<std::set<IdType>(IdType)> &neighbors)
{
  quick_union<IdType> uf(nn);
  contour_tree<IdType> jt;

  for (IdType i=nn-1; i<nn; i--) {
    jt.add_node(order[i]);
    bool is_minima = true;
    for (auto j : neighbors(order[i])) {
      j = inverse_order[j];
      // fprintf(stderr, "i=%zu(%zu), j=%zu(%zu)\n", order[i], i, order[j], j);
      if (j > i) {
        IdType ri = uf.root(i), rj = uf.root(j);
        if (ri != rj) {
          jt.add_arc(order[ri], order[rj]);
          uf.unite(rj, ri);
        } 
      } else
        is_minima = false;
    }
    if (is_minima)
      fprintf(stderr, "minima: %zu\n", order[i]);
  }

  return jt;
}

template <class IdType>
contour_tree<IdType> build_split_tree(IdType nn, 
    const std::vector<IdType> &order,
    const std::vector<IdType> &inverse_order,
    const std::function<std::set<IdType>(IdType)> &neighbors)
{
  quick_union<IdType> uf(nn);
  contour_tree<IdType> st;

  for (IdType i=0; i<nn; i++) {
    st.add_node(order[i]);
    bool is_maxima = true;
    for (auto j : neighbors(order[i])) {
      j = inverse_order[j];
      if (j < i) {
        IdType ri = uf.root(i), rj = uf.root(j);
        if (ri != rj) {
          st.add_arc(order[rj], order[ri]);
          uf.unite(rj, ri);
        }
      } else 
        is_maxima = false;
    }
    if (is_maxima) 
      fprintf(stderr, "maxima: %zu\n", order[i]);
  }

  return st;
}

template <class IdType>
contour_tree<IdType> sweep_and_merge(IdType nn,
    const std::vector<IdType> &order,
    const std::vector<IdType> &inverse_order,
    const std::function<std::set<IdType>(IdType)> &neighbors)
{
  contour_tree<IdType> 
    st = build_split_tree(nn, order, inverse_order, neighbors), 
    jt = build_join_tree(nn, order, inverse_order, neighbors);
  
  return merge_join_and_split_trees<IdType>(nn, jt, st);
}

template <class IdType, class ValueType>
contour_tree<IdType> build_contour_tree(IdType nn,
    const std::function<ValueType(IdType)> &value,
    const std::function<std::set<IdType>(IdType)> &neighbors)
{
  std::vector<IdType> order(nn), inverse_order(nn);
  for (IdType i=0; i<nn; i++) 
    order[i] = i;

  std::stable_sort(order.begin(), order.end(), 
      [&](IdType p, IdType q) {
        return value(p) < value(q);
      });

  for (IdType i=0; i<nn; i++)
    inverse_order[order[i]] = i;

  return sweep_and_merge<IdType>(nn, order, inverse_order, neighbors);
}
  
template <class IdType, class ValueType>
contour_tree<IdType> build_contour_tree(IdType nn,
    const std::vector<ValueType> &vector,
    const std::function<std::set<IdType>(IdType)> &neighbors)
{
  // auto value = std::bind(&std::vector<double>::at, vector, std::placeholders::_1);
  return build_contour_tree<IdType, ValueType>(nn, 
      [&vector](IdType i) {return vector[i];},
      neighbors);
}

} // namespace ftk

#endif
