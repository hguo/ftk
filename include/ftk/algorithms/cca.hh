#ifndef _FTK_CCA_H
#define _FTK_CCA_H

#include <set>
#include <map>
#include <string>

#include <ftk/basic/union_find.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>

namespace ftk {

// ================================================================
// Extract connected components by given mesh elements

#if 0
typedef simplicial_regular_mesh_element EleType;
template <class ContainerType>
std::vector<std::set<EleType> > extract_connected_components_element(
    const std::function<ContainerType(EleType) >& neighbors,
    const std::set<EleType> &qualified_)
{
  std::set<EleType> qualified(qualified_);

  union_find<std::string> UF; 
  std::map<std::string, EleType> id2ele; 
  // for(auto ite = qualified.begin(); ite != qualified.end(); ++ite) {
  for(auto& ele : qualified) {
    std::string id = ele.to_string(); 

    UF.add(id); 
    id2ele.insert(std::make_pair (id, ele)); 
  }

  for(auto ite = qualified.begin(); ite != qualified.end(); ++ite) {
    std::string current = ite->to_string(); 

    for (auto nei_ite : neighbors(*ite)) {
      std::string neighbor = nei_ite.to_string(); 

      if (UF.has(neighbor)) {
        UF.unite(current, neighbor); 
      }
    }
  }

  std::vector<std::set<std::string> > components_str = UF.get_sets();

  // Convert element ids to elements
  std::vector<std::set<EleType> >  components; 
  for(auto comp_str = components_str.begin(); comp_str != components_str.end(); ++comp_str) {
    std::set<EleType> comp; 
    for(auto ele_id = comp_str->begin(); ele_id != comp_str->end(); ++ele_id) {
      comp.insert(id2ele.find(*ele_id)->second); 
    }

    components.push_back(comp); 
  }

  return components;
}
#endif

// ================================================================
// Extract connected components by given IDs

template <class IdType>
std::set<std::pair<IdType, IdType> > track_connected_components(
    const std::vector<std::set<IdType> > &components0,
    const std::vector<std::set<IdType> > &components1) 
{
  auto overlaps = [](const std::set<IdType>& s0, const std::set<IdType>& s1) {
    for (const auto &e0 : s0) 
      if (s1.find(e0) != s1.end()) return true;
    return false;
  };

  std::set<std::pair<IdType, IdType> > results;

  for (size_t i = 0; i < components0.size(); i ++) {
    for (size_t j = 0; j < components1.size(); j ++) {
      if (overlaps(components0[i], components1[j])) {
        results.insert(std::make_pair(i, j));
      }
    }
  }

  return results;
}

template <class IdType, class ContainerType>
std::vector<std::set<IdType> > extract_connected_components(
    const std::function<ContainerType(IdType) >& neighbors,
    const std::set<IdType> &qualified_)
{
  std::set<IdType> qualified(qualified_);

  union_find<IdType> UF; 
  for(auto ite = qualified.begin(); ite != qualified.end(); ++ite) {
    UF.add(*ite); 
  }

  for(auto ite = qualified.begin(); ite != qualified.end(); ++ite) {
    IdType current = *ite; 

    for (auto neighbor : neighbors(current)) {
      if (UF.has(neighbor)) {
        UF.unite(current, neighbor); 
      }
    }
  }  

  std::vector<std::set<IdType> > components = UF.get_sets();

  return components;
}

template <class IdType, class ContainerType>
std::vector<std::set<IdType> > extract_connected_components(
    IdType nNodes,
    const std::function<ContainerType(IdType) >& neighbors,
    const std::function<bool(IdType)>& criterion)
{
  // fprintf(stderr, "finding qualified...\n");
  // find qualified
  std::set<IdType> qualified;
  for (IdType i=0; i<nNodes; i++) 
    if (criterion(i)) 
      qualified.insert(i);
  
  fprintf(stderr, "#qualified=%zu\n", qualified.size());

  return extract_connected_components(neighbors, qualified);
}

template <class IdType, class ContainerType>
std::vector<std::set<IdType> > extract_connected_components(
    const ContainerType& nodes,
    const std::function<ContainerType(IdType) >& neighbors,
    const std::function<bool(IdType)>& criterion)
{
  std::set<IdType> qualified;
  for (const auto &n : nodes)
    if (criterion(n))
      qualified.insert(n);

  return extract_connected_components(neighbors, qualified);
}

}

#endif
