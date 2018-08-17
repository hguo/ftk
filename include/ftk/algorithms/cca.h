#ifndef _CCA_H
#define _CCA_H

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

#endif
