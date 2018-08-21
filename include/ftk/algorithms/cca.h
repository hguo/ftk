#ifndef _CCA_H
#define _CCA_H

// #include "ftk/transition/transitionMatrix.h"
#include "ftk/algorithms/bfs.h"

namespace ftk {

template <class IdType>
std::vector<std::pair<IdType, IdType> > trackConnectedComponents(
    const std::vector<std::set<IdType> > &components0,
    const std::vector<std::set<IdType> > &components1) 
{
  auto overlaps = [](const std::set<IdType>& s0, const std::set<IdType>& s1) {
    for (const auto &e0 : s0) 
      if (s1.find(e0) != s1.end()) return true;
    return false;
  };

  std::vector<std::pair<IdType, IdType> > mat;

  for (size_t i = 0; i < components0.size(); i ++) {
    for (size_t j = 0; j < components1.size(); j ++) {
      if (overlaps(components0[i], components1[j])) {
        mat.push_back(std::make_pair(i, j));
      }
    }
  }

  return mat;
}


template <class IdType>
std::vector<std::set<IdType> > extractConnectedComponents(
    const std::function<std::set<IdType>(IdType) >& neighbors,
    const std::set<IdType> &qualified_)
{
#if 0
  // extract connected components
  // std::set<IdType> qualified(qualified_);
  std::vector<std::set<IdType> > components;
  std::set<IdType> visited;

  while (!qualified.empty()) {
    IdType seed = *qualified.begin();
    std::set<IdType> component;

    bfs<IdType>(seed, neighbors, 
        [&component, &visited](IdType i) {component.insert(i); visited.insert(i);}, 
        [&qualified, &visited](IdType i) {return qualified.find(i) != qualified.end() && visited.find(i) == visited.end();});
    
    components.emplace_back(component);
  }

  return components;
#endif

  std::set<IdType> qualified(qualified_);
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
    const std::function<std::set<IdType>(IdType) >& neighbors,
    const std::function<bool(IdType)>& criterion)
{
  // fprintf(stderr, "finding qualified...\n");
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
