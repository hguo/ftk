#ifndef _CCA_H
#define _CCA_H

#include "ftk/transition/transitionMatrix.h"

namespace ftk {

template <class IdType>
TransitionMatrix associateConnectedComponents(
    const std::vector<std::set<IdType> > &components0,
    const std::vector<std::set<IdType> > &components1) 
{
  auto checkIntersection = [](const std::set<IdType>& s0, const std::set<IdType>& s1) {
    for (const auto &e0 : s0) {
      if (s1.find(e0) != s1.end()) return true;
    }
    return false;
  };

  TransitionMatrix mat;

  for (size_t i = 0; i < components0.size(); i ++) {
    for (size_t j = 0; j < components1.size(); j ++) {
      if (checkIntersection(components0[i], components1[j])) {
        mat.set(i, j);
      }
    }
  }
}

template <class IdType>
std::vector<std::set<IdType> > extractConnectedComponents(
    const std::function<std::set<IdType>(size_t) >& neighbors,
    std::set<IdType> qualified)
{
  // extract connected components
  std::vector<std::set<IdType> > components;

  std::set<IdType> Q;

  while (!qualified.empty()) {
    Q.insert(*qualified.begin());

    std::set<IdType> visited;
    while (!Q.empty()) {
      IdType current = *Q.begin();
      Q.erase(current);
      visited.insert(current);

      for (auto neighbor : neighbors(current)) {
        if (qualified.find(neighbor) != qualified.end()
            && visited.find(neighbor) == visited.end()
            && Q.find(neighbor) == Q.end()) {
          Q.insert(neighbor);
        }
      }
    }

    for (auto v : visited)
      qualified.erase(v);

    components.push_back(visited);
  }

  return components;
}

template <class IdType>
std::vector<std::set<IdType> > extractConnectedComponents(
    IdType nNodes,
    const std::function<std::set<IdType>(size_t) >& neighbors,
    const std::function<bool(IdType)>& criterion)
{
  fprintf(stderr, "finding qualified...\n");
  // find qualified
  std::set<IdType> qualified;
  for (IdType i=0; i<nNodes; i++) 
    if (criterion(i)) 
      qualified.insert(i);
  
  fprintf(stderr, "#qualified=%zu\n", qualified.size());

  return extractConnectedComponents(neighbors, qualified);
}

}

#endif
