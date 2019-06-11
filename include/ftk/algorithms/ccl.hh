#ifndef _CCL_H
#define _CCL_H

#include <set>
#include <vector>
#include "ftk/algorithms/bfs.h"

namespace ftk {

template <class IdType=size_t, class ArrayType=std::vector<size_t>, class ContainerType=std::set<size_t>, class LabelType=size_t>
LabelType labelConnectedComponents(
    IdType nNodes,
    ArrayType &labels,
    const std::function<ContainerType(IdType) > &neighbors,
    const std::function<bool(IdType)> &criterion)
{
  LabelType componentCount = 0;

  for (IdType seed = 0; seed < nNodes; seed ++) {
    if (labels[seed] == 0 && criterion(seed)) {
      componentCount ++;
      bfs<IdType, ContainerType>(
          seed,
          neighbors, 
          [&labels, componentCount](IdType i) {
            labels[i] = componentCount;
          }, criterion);
    }
  }

  return componentCount;
}

template <class LabelType, class ContainerType>
std::set<std::pair<LabelType, LabelType> > trackConnectedComponentsByLabels(
    const ContainerType& labels0, 
    const ContainerType& labels1) 
{
  std::set<std::pair<LabelType, LabelType> > results;

  for (size_t i=0; i<labels0.size(); i++) {
    auto p = labels0[i], q = labels1[i];
    if (p != q && p != 0 && q != 0)
      results.insert(std::make_pair(p, q));
  }

  return results;
}

}

#endif
