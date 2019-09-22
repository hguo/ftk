#ifndef _FTK_DISTRIBUTED_UNION_FIND_H
#define _FTK_DISTRIBUTED_UNION_FIND_H

#include <vector>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <iostream>

#include <ftk/external/diy/mpi.hpp>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/assigner.hpp>
#include <ftk/external/diy/serialization.hpp>

// Union-find algorithm with distributed-memory parallelism
  // All elements need to be added by invoking the add function at the initialization stage. 

// Reference 
  // Paper: "Evaluation of connected-component labeling algorithms for distributed-memory systems"

// Add the sparse representation by using Hash map/table 

#define ISDEBUG   0
#define OUTPUT_TIME_EACH_ROUND 0

namespace ftk {
template <class IdType=std::string>
struct distributed_union_find
{
  distributed_union_find() : eles(), id2parent(), parent2localchildren(), parent2nonlocalchildren() {

  }

  // Initialization
  // Add and initialize elements
  void add(IdType i) {
    eles.insert(i); 
    id2parent.insert(std::make_pair(i, i)); 
    parent2localchildren.insert(std::make_pair(i, std::vector<IdType>())); 
    parent2nonlocalchildren.insert(std::make_pair(i, std::vector<IdType>())); 
  }

  // Operations

  void set_parent(IdType i, IdType par) {
    assert(has(i));

    // Ensure the id of new parent is smaller than the id of old parent
      // Otherwise, may occur override and erase the updated parent
    if(par < id2parent[i]) {
      id2parent[i] = par; 
    }
    // else {
    //   std::cout<< "Parent is override! " <<std::endl;
    // }
  }

  // void add_child(IdType par, IdType ele) {
  //   assert(has(par));

  //   if(this->has(ele)) {
  //     parent2localchildren[par].push_back(ele);   
  //   } else {
  //     parent2nonlocalchildren[par].push_back(ele);   
  //   }
  // }

  void add_local_child(IdType par, IdType ele) {
    assert(has(par));

    parent2localchildren[par].push_back(ele);   
  }

  void add_nonlocal_child(IdType par, IdType ele) {
    assert(has(par));

    parent2nonlocalchildren[par].push_back(ele);   
  }

  void add_nonlocal_intermediate_root(IdType intermediate_root) {
    nonlocal_intermediate_roots.insert(intermediate_root); 
  }

  void erase_nonlocal_intermediate_root(IdType intermediate_root) {
    nonlocal_intermediate_roots.erase(intermediate_root); 
  }

  // Queries

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
  }

  IdType parent(IdType i) {
    assert(has(i));

    return id2parent[i]; 
  }

  std::vector<IdType>& local_children(IdType i) {
    assert(has(i));

    return parent2localchildren[i]; 
  }

  std::vector<IdType>& nonlocal_children(IdType i) {
    assert(has(i));

    return parent2nonlocalchildren[i]; 
  }

  typename std::vector<IdType>::iterator erase_local_child(IdType par, typename std::vector<IdType>::iterator ite) {
    return parent2localchildren[par].erase(ite); 
  }

  typename std::vector<IdType>::iterator erase_nonlocal_child(IdType par, typename std::vector<IdType>::iterator ite) {
    return parent2nonlocalchildren[par].erase(ite); 
  }

  void clear_local_children(std::string i) {
    assert(has(i));

    parent2localchildren[i].clear(); 
  }

  void clear_nonlocal_children(std::string i) {
    assert(has(i));

    parent2nonlocalchildren[i].clear(); 
  }

  bool is_nonlocal_intermediate_root(IdType i) {
    if(this->nonlocal_intermediate_roots.find(i) != this->nonlocal_intermediate_roots.end()) {
      return true ;
    }

    return false ;
  }

  bool is_root(IdType i) {
    assert(has(i));

    return i == id2parent[i]; 
  }

public:
  std::set<IdType> eles; 

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id2parent;

  // std::map<IdType, std::vector<IdType>> parent2children;

  std::map<IdType, std::vector<IdType>> parent2localchildren;
  std::map<IdType, std::vector<IdType>> parent2nonlocalchildren;

  std::set<IdType> nonlocal_intermediate_roots;
  
};

}

// DIY Block for distributed union-find
struct Block_Union_Find : public ftk::distributed_union_find<std::string> {
  Block_Union_Find(): nchanges(0), related_elements(), all_related_elements(), ele2gid(), distributed_union_find() { 
    
  }

  // add an element
  void add(std::string ele) {
    this->nchanges += 1;

    if(this->has(ele)) {
      std::cout<<"This ele has been added. "<<ele<<std::endl; 
    } else {
      distributed_union_find::add(ele); 
      this->related_elements.insert(std::make_pair(ele, std::set<std::string>())); 
    }
  }

  void erase_element(std::string ele) {
    this->nchanges += 1;

    this->eles.erase(ele);
  }

  bool has_related_element(std::string ele, std::string related_ele) {
    assert(this->has(ele)); 

    if(this->related_elements[ele].find(related_ele) == this->related_elements[ele].end()) {
      return false; 
    }

    return true; 
  }

  void add_related_element(std::string ele, std::string related_ele) {
    this->nchanges += 1;

    // std::cout<<"Add union: "<<ele<<"-"<<related_ele<<std::endl; 

    assert(this->has(ele)); 

    if(related_ele < ele) {
      
      if(this->all_related_elements[ele].find(related_ele) == this->all_related_elements[ele].end()) {
        this->related_elements[ele].insert(related_ele); 
        this->all_related_elements[ele].insert(related_ele); 
      }

    } else {
      std::cout<<"Related element is larger! "<<ele<<" "<<related_ele<<std::endl; 
      exit(0);   
    }

  }

  const std::set<std::string>& get_related_elements(std::string ele) {
    assert(this->has(ele)); 
    
    return this->related_elements[ele]; 
  }

  void clear_related_elements(std::string ele) {
    this->nchanges += 1;

    assert(this->has(ele)); 

    this->related_elements[ele].clear(); 
  }

  void erase_related_element(std::string ele, std::string related_element) {
    this->nchanges += 1;

    assert(this->has(ele)); 
      
    this->related_elements[ele].erase(related_element); 
  }

  bool has_gid(std::string ele) {
    if(this->ele2gid.find(ele) == this->ele2gid.end()) {
      return false; 
    }

    return true; 
  }

  void set_gid(std::string ele, int gid) {
    this->nchanges += 1;

    if(gid >= 0) {
      this->ele2gid[ele] = gid; 
    }
  }

  int get_gid(std::string ele) {
    assert(this->has_gid(ele)); 

    return this->ele2gid[ele]; 
  }

  void set_parent(std::string i, std::string par) {
    this->nchanges += 1;

    distributed_union_find::set_parent(i, par); 
  }

  void add_local_child(std::string par, std::string ele) {
    this->nchanges += 1;

    distributed_union_find::add_local_child(par, ele); 
  }

  void add_nonlocal_child(std::string par, std::string ele) {
    this->nchanges += 1;

    distributed_union_find::add_nonlocal_child(par, ele); 
  }


  std::vector<std::string>::iterator erase_local_child(std::string par, std::vector<std::string>::iterator ite) {
    this->nchanges += 1;

    return distributed_union_find::erase_local_child(par, ite); 
  }

  std::vector<std::string>::iterator erase_nonlocal_child(std::string par, std::vector<std::string>::iterator ite) {
    this->nchanges += 1;

    return distributed_union_find::erase_nonlocal_child(par, ite); 
  }

  void clear_local_children(std::string par) {
    this->nchanges += 1;

    distributed_union_find::clear_local_children(par); 
  }

  void clear_nonlocal_children(std::string par) {
    this->nchanges += 1;

    distributed_union_find::clear_nonlocal_children(par); 
  }

  void add_nonlocal_intermediate_root(std::string intermediate_root) {
    this->nchanges += 1;

    distributed_union_find::add_nonlocal_intermediate_root(intermediate_root); 
  }

  void erase_nonlocal_intermediate_root(std::string intermediate_root) {
    this->nchanges += 1;

    distributed_union_find::erase_nonlocal_intermediate_root(intermediate_root); 
  }

  bool is_intermediate_root(std::string i) {
    if(this->has(i)) {
      if(this->is_root(i) && this->get_related_elements(i).size() == 0) {
        return true ;
      } else {
        return false ;  
      }
    } else {
      return distributed_union_find::is_nonlocal_intermediate_root(i); 
    }

    return false ;
  }


  void get_sets(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results);

public: 

  std::map<std::string, int> ele2gid; 

  int nchanges = 0; // # of processed unions per round = valid unions (united unions) + passing unions

  #if OUTPUT_TIME_EACH_ROUND
    double time = 0; // time for measure duration of each round;
  #endif
  
  // std::map<int, std::set<std::string>> cache_send_intermediate_roots ;

private:
  // map element id to ids of its related elements
  // Can be optimized by ordered the related elements, put related elements on process first and ordered decreingly by ids
  
  std::map<std::string, std::set<std::string>> related_elements; 
  std::map<std::string, std::set<std::string>> all_related_elements;   
};

// ==========================================================

struct Message_Union_Find {

  Message_Union_Find() : tag(), strs() {

  }

// Send message for query
  void send_gid_query(std::string& ele) {
    tag = "gid_query"; 
    strs.push_back(ele);  
  }

  void send_gid_response(std::string& ele, int& gid) {
    tag = "gid_response"; 

    strs.push_back(ele); 
    strs.push_back(std::to_string(gid)); 
  }

  void send_gparent(std::string& ele, std::string& gparent, int& gid_gparent, bool& is_known) {
    tag = "gparent";

    strs.push_back(ele); 
    strs.push_back(gparent); 
    strs.push_back(std::to_string(gid_gparent)); 
    strs.push_back(std::to_string(is_known)); 
  }

  void send_child(const std::string& parent, const std::string& child, const int& gid_child, const bool& is_known) {
    tag = "child";

    strs.push_back(parent); 
    strs.push_back(child); 
    strs.push_back(std::to_string(gid_child)); 
    strs.push_back(std::to_string(is_known)); 
  }

  void send_add_intermediate_root(const std::string& parent) {
    tag = "add_intermediate_root";

    strs.push_back(parent); 
  }

  void send_erase_intermediate_root(const std::string& parent) {
    tag = "erase_intermediate_root";

    strs.push_back(parent); 
  }

  // Send the union (ele, related_ele) to ele
  void send_union(const std::string& ele, const std::string& related_ele, const int& gid_related_ele) {
    tag = "union"; 

    strs.push_back(ele); 
    strs.push_back(related_ele); 
    strs.push_back(std::to_string(gid_related_ele)); 
  }

  void send_ele_parent_pair(std::pair<std::string, std::string>& pair) {
    tag = "ele_parent_pair"; 

    strs.push_back(pair.first); 
    strs.push_back(pair.second); 
  }

  void send_ele_parent_pairs(std::vector<std::pair<std::string, std::string>>& pairs) {
    tag = "ele_parent_pairs"; 

    for(auto& pair : pairs) {
      strs.push_back(pair.first); 
      strs.push_back(pair.second); 
    }
  }

// Receive message

  void rec_gid_query(std::string& ele) {
    ele = strs[0]; 
  }

  void rec_gid_response(std::string& ele, int& gid) {
    ele = strs[0]; 
    gid = std::stoi(strs[1]);  
  }

  void rec_gparent_query(std::string& ele, std::string& parent) {
    ele = strs[0]; 
    parent = strs[1]; 
  }

  void rec_gparent(std::string& ele, std::string& gparent, int& gid_gparent, bool& is_known) {
    ele = strs[0]; 
    gparent = strs[1];

    gid_gparent = std::stoi(strs[2]); 
    is_known = std::stoi(strs[3]); 
  }

  void rec_child(std::string& parent, std::string& child, int& gid_child, bool& is_known) {
    parent = strs[0]; 
    child = strs[1];

    gid_child = std::stoi(strs[2]); 
    is_known = std::stoi(strs[3]); 
  }

  void rec_add_intermediate_root(std::string& parent) {
    parent = strs[0];
  }

  void rec_erase_intermediate_root(std::string& parent) {
    parent = strs[0];
  }

  void rec_union(std::string& ele, std::string& related_ele, int& rgid) {
    ele = strs[0]; 
    related_ele = strs[1]; 

    rgid = std::stoi(strs[2]); 
  }

  void receive_ele_parent_pair(std::pair<std::string, std::string>& pair) {
    pair = std::make_pair(strs[0], strs[1]); 
  }

  void receive_ele_parent_pairs(std::vector<std::pair<std::string, std::string>>& pairs) {
    for(int i = 0; i < strs.size(); i += 2) {
      pairs.push_back(std::make_pair(strs[i], strs[i+1]));
    }
  }

public:
  std::string tag; 

  std::vector<std::string> strs; 
}; 

// ==========================================================

namespace diy
{
    template<>
        struct Serialization<Message_Union_Find>
    {
        static void save(BinaryBuffer& bb, const Message_Union_Find& msg)
        {
            diy::save(bb, msg.tag);
            // diy::save(bb, msg.str);
            diy::save(bb, msg.strs);
        }

        static void load(BinaryBuffer& bb, Message_Union_Find& msg)
        {
            diy::load(bb, msg.tag);
            // diy::load(bb, msg.str);
            diy::load(bb, msg.strs);
        }
    };
}

// ==========================================================

void import_data(std::vector<Block_Union_Find*>& blocks, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<int>& gids) {
  // for the local blocks in this processor
  for (unsigned i = 0; i < gids.size(); ++i) {
    int gid = gids[i];

    diy::Link*    link = new diy::Link; 
    for(unsigned j = 0; j < assigner.nblocks(); ++j) {
      // int neighbor_gid = gids[j]; 
      int neighbor_gid = j; 
      if(gid != neighbor_gid) {
        diy::BlockID  neighbor;
        neighbor.gid  = neighbor_gid;  // gid of the neighbor block
        neighbor.proc = assigner.rank(neighbor.gid);  // process of the neighbor block
        link->add_neighbor(neighbor); 
      }
    }

    Block_Union_Find* b = blocks[i]; 
    master.add(gid, b, link); 
    // master.replace_link(0, link);
  }
}

void query_gid(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  // std::cout<<"Query: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 
  
  // Can be optimized more, link by a queue. 
  for(auto it = b->ele2gid.begin(); it != b->ele2gid.end(); ++it) {
    if(it->second == -1) { // If the gid of this element is unknown

      // std::cout<<"A gid is missing: " << it->first<<std::endl; 

      std::string ele = it->first; 
      if(b->has(ele)) {
        it->second = gid; 
      } else {
        for (int i = 0; i < l->size(); ++i) {
          Message_Union_Find msg = Message_Union_Find(); 
          msg.send_gid_query(ele); 

          cp.enqueue(l->target(i), msg);
        }
      }

    }
  }

  // std::cout<<std::endl; 
}


void answer_gid(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg, int from_gid) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();
  // diy::Master* master = cp.master(); 

  // std::cout<<"Answer: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 

  std::string ele; 
  msg.rec_gid_query(ele); 

  if(b->has(ele)) {
    Message_Union_Find send_msg; 
    send_msg.send_gid_response(ele, gid); 

    // std::cout<<ele<<"-";
    // std::cout<<gid;
    // std::cout<<std::endl; 

    cp.enqueue(l->target(l->find(from_gid)), send_msg);
    // cp.enqueue(l->target(l->find(in[i])), 1);
  }

  // std::cout<<std::endl; 
}

void save_gid(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  // std::cout<<"Save: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 

  std::string ele; 
  int gid_ele; 
  msg.rec_gid_response(ele, gid_ele); 

  // std::cout<<_pair.first<<" - "<<_pair.second<<std::endl; 

  b->set_gid(ele, gid_ele); 

  // std::cout<<"i-"<<i<<'-'<<in[i]<<std::endl;
  // int t; 
  // cp.dequeue(in[i], t); 
  // std::cout<<t<<std::endl; 


  // std::cout<<std::endl; 
}


void unite_once(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  for(auto& ele : b->eles) {

    if(b->is_root(ele)) {
      auto& related_elements = b->get_related_elements(ele); 
      // cannot use auto&, since we will remove some elements from the original set; auto& will cause segmentation erro

      std::string found_related_ele = ele; // Find a smallest related element has a smaller id than ele

      // Find from local related elements
      for(auto& related_ele : related_elements) {
        if(b->has(related_ele)) {
          if (related_ele < found_related_ele) {
            found_related_ele = related_ele; 
          }
        }
      }

      if(found_related_ele == ele) { // Local related elements do not have the satisfied one
        // Find from non-local related elements
        for(auto& related_ele : related_elements) {
          if(!b->has(related_ele) && b->has_gid(related_ele)) {
            if (related_ele < found_related_ele) {
              found_related_ele = related_ele; 
            }
          }
        }
      }

      // auto& related_elements = b->get_related_elements(ele); 
      // std::string found_related_ele = ele; // Find a smallest related element has a smaller id than ele
      // for(auto& related_ele : related_elements) {
        
      //   if(!b->has_gid(related_ele)) {
      //     continue ; 
      //   }

      //   // Unite a root element with larger id or smaller id is better? 
      //     // Here is a smaller id

      //   if (related_ele < found_related_ele) {
      //     found_related_ele = related_ele; 
      //   }
      // }

      if(found_related_ele < ele) {
        int rgid = b->get_gid(found_related_ele); 
        // std::cout<<gid<<" "<<rgid<<std::endl; 
        
        b->set_parent(ele, found_related_ele); 

        if(b->has(found_related_ele)) {
          b->add_local_child(found_related_ele, ele); 
        } else {
          Message_Union_Find send_msg; 

          bool is_known = b->is_intermediate_root(found_related_ele); 
          send_msg.send_child(found_related_ele, ele, gid, is_known); 

          cp.enqueue(l->target(l->find(rgid)), send_msg);
        }

        #if ISDEBUG
          std::cout<<ele<<" -1> "<<found_related_ele<<std::endl; 
        #endif

        b->erase_related_element(ele, found_related_ele); 
      } 

    }
  }

}


void compress_path(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  diy::Link* l = cp.link();

  for(auto& parent : b->eles) {
    if(!b->is_root(parent)) {
      std::string grandparent = b->parent(parent); 
      bool is_local_grandparent = b->has(grandparent); // if the grandparent is not in this block
      if(!is_local_grandparent && !b->has_gid(grandparent)) { // if currently the gid of grandparent is not available
        continue ;
      }
      int gid_grandparent = b->get_gid(grandparent); 

      {
        auto& nonlocal_children = b->nonlocal_children(parent); 
        if(nonlocal_children.size() > 0) {
          std::vector<std::string> cache; 
          for(auto& child : nonlocal_children) {

            if(!b->has_gid(child)) { // if currently the gid of child is not available
              cache.push_back(child); 
              continue ;
            }

            int gid_child = b->get_gid(child); 

            // Set child's parent to grandparent
            Message_Union_Find send_msg; 
            bool is_known = b->is_intermediate_root(grandparent); 
            send_msg.send_gparent(child, grandparent, gid_grandparent, is_known); 
            // std::cout<<*ele_ptr<<" - "<<*parent_ptr<<" - "<<grandparent<<" - "<<b->ele2gid[grandparent]<<std::endl; 
            cp.enqueue(l->target(l->find(gid_child)), send_msg); 

            // Add child to grandparent's child list
            if(is_local_grandparent) {
              b->add_nonlocal_child(grandparent, child); 
            } else {
              Message_Union_Find send_msg; 

              bool is_known = b->is_intermediate_root(grandparent); 
              send_msg.send_child(grandparent, child, gid_child, is_known); 
              cp.enqueue(l->target(l->find(gid_grandparent)), send_msg); 
            }
          }
          b->clear_nonlocal_children(parent); 

          if(cache.size() > 0) {
            for(auto& cache_ele : cache) {
              b->add_nonlocal_child(parent, cache_ele); 
            }
          }
        }
      }


      if(is_local_grandparent || b->is_intermediate_root(grandparent)) {
        auto& local_children = b->local_children(parent); 
        if(local_children.size() > 0) {
          for(auto& child : local_children) {

            // Set child's parent to grandparent
            b->set_parent(child, grandparent); 

            // Add child to grandparent's child list
            if(is_local_grandparent) {
              b->add_local_child(grandparent, child); 
            } else {
              int gid_child = b->get_gid(child); 

              Message_Union_Find send_msg; 
              bool is_known = b->is_intermediate_root(grandparent); 
              send_msg.send_child(grandparent, child, gid_child, is_known); 

              cp.enqueue(l->target(l->find(gid_grandparent)), send_msg); 
            }
          }
          b->clear_local_children(parent); 

        }
      }

    }
  }
}


// Distributed path compression - Step Two
void distributed_save_gparent(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  std::string ele; 
  std::string grandpar; 
  int gid_grandparent; 
  bool is_known; 

  msg.rec_gparent(ele, grandpar, gid_grandparent, is_known); 
  
  // std::cout<<*ele_ptr << " - " << b->parent(*ele_ptr) <<" - "<< *grandpar_ptr<<" - "<<*grandpar_gid_ptr<<std::endl; 

  // if(*ele_ptr == "4") {
  //   std::cout<<*ele_ptr<<" ~ "<<*grandpar_ptr<<std::endl; 
  // }


  b->set_parent(ele, grandpar); 
  b->set_gid(grandpar, gid_grandparent); 
  if(is_known && !b->has(grandpar)) {
    b->add_nonlocal_intermediate_root(grandpar);
  }

  #if ISDEBUG
    std::cout<<ele<<" -2> "<<grandpar<<std::endl; 
  #endif
}

// Distributed path compression - Step Three
void distributed_save_child(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  diy::Link* l = cp.link();

  std::string par; 
  std::string child; 
  int gid_child; 
  
  // Indicate whether child knows the new parent is an intermediate root or not: True - child knows; False - child does not know
    // If child does not know, by default, child thinks the new parent is NOT an intermediate root. 
  bool is_known; 

  msg.rec_child(par, child, gid_child, is_known); 

  b->set_gid(child, gid_child); 

  // Possible that the child has the same process as the new parent

  if(b->has(child)) {
    b->add_local_child(par, child); 
  } else { // if child is not in this block
    b->add_nonlocal_child(par, child); 

    if(is_known) {  // if child knows this parent is an intermediate root
      if(!b->is_intermediate_root(par))  { // however, currently it is not
        // Tell the child that it is not the intermediate root anymore

        // b->cache_send_intermediate_roots[gid_child].insert(par); 

        Message_Union_Find send_msg; 
        send_msg.send_erase_intermediate_root(par); // Add an intermediate root to the process of child

        cp.enqueue(l->target(l->find(gid_child)), send_msg); 
      }
    } else { // if child does not know
      if(b->is_intermediate_root(par)) { // however, currently it is
        // Tell the child that it is the intermediate root

        Message_Union_Find send_msg; 
        send_msg.send_add_intermediate_root(par); // Add an intermediate root to the process of child

        cp.enqueue(l->target(l->find(gid_child)), send_msg); 
      }
    }
  }  

}


// Receive adding an intermediate root
void distributed_add_intermediate_root(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  std::string par;
  msg.rec_add_intermediate_root(par); // Add an intermediate root

  b->add_nonlocal_intermediate_root(par); 
}

void distributed_erase_intermediate_root (Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  std::string parent;
  msg.rec_erase_intermediate_root(parent); // Erase an intermediate root

  b->erase_nonlocal_intermediate_root(parent); 
}

void pass_unions(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  // Local computation
  // Pass unions of elements in this block to their parents, save these unions
  for(auto& ele : b->eles) {

    if(!b->is_root(ele)) {

      auto& src = b->get_related_elements(ele);
      if(src.size() > 0) {
        std::string par = b->parent(ele); 
        if(!b->is_intermediate_root(par)) { // Previously test best
          continue ;
        }
        // if(!is_local_parent && !b->is_intermediate_root(par)) { 
        //   continue ;
        // }

        bool is_local_parent = true; 
        int p_gid = gid; 
        if(!b->has(par)) {
          if(!b->has_gid(par)) {
            continue ;
          }

          is_local_parent = false;
          p_gid = b->get_gid(par); 
        }

        std::vector<std::string> cache;
        for(auto& related_ele : src) {

          if(related_ele < par) {

            if(is_local_parent) {
              b->add_related_element(par, related_ele); 
            } else {
              // Communicate with other processes

              if(!b->has_gid(related_ele)) {
                cache.push_back(related_ele); 
                continue ;
              }

              int r_gid = b->get_gid(related_ele); 

              Message_Union_Find send_msg; 
              send_msg.send_union(par, related_ele, r_gid); 

              cp.enqueue(l->target(l->find(p_gid)), send_msg); 

            }

          } else if(par < related_ele) {

            if(b->has(related_ele)) {
              b->add_related_element(related_ele, par); 
            } else {
              if(!b->has_gid(related_ele)) {
                cache.push_back(related_ele); 
                continue ;
              }

              int r_gid = b->get_gid(related_ele); 

              Message_Union_Find send_msg; 
              send_msg.send_union(related_ele, par, p_gid); 

              cp.enqueue(l->target(l->find(r_gid)), send_msg); 
            }

          }

        }

        b->clear_related_elements(ele); 

        if(cache.size() > 0) {
          for(auto& cache_ele : cache) {
            b->add_related_element(ele, cache_ele); 
          }
        }

      }

    }

  }
}


// Update unions of related elements
void distributed_save_union(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  // std::cout<<"Save Unions: "<<std::endl; 
  // std::cout<<"Block ID: "<<cp.gid()<<std::endl; 

  std::string ele; 
  std::string related_ele; 
  int rgid; 
  
  msg.rec_union(ele, related_ele, rgid); 

  // std::cout<<*ele_ptr<<" - "<<*par_ptr << " - " << *related_ele_ptr <<" - "<< *gid_ptr<<std::endl; 

  b->add_related_element(ele, related_ele); 
  b->set_gid(related_ele, rgid); 

  // std::cout<<std::endl; 
}

void local_computation(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {

  #if ISDEBUG
    int gid = cp.gid();
    std::cout<<"Start Local Computation. "<<"gid: "<<gid<<std::endl; 
  #endif

  unite_once(b, cp); 

  #if ISDEBUG
    int gid = cp.gid();
    std::cout<<"Finish unite once. "<<"gid: "<<gid<<std::endl; 
  #endif

  compress_path(b, cp); 

  #if ISDEBUG
    int gid = cp.gid();
    std::cout<<"Finish compress_path. "<<"gid: "<<gid<<std::endl; 
  #endif

  pass_unions(b, cp); 

  #if ISDEBUG
    int gid = cp.gid();
    std::cout<<"Finish pass_unions. "<<"gid: "<<gid<<std::endl; 
  #endif


  // diy::Link* l = cp.link();
  // for(auto& pair : b->cache_send_intermediate_roots) {
  //   if(pair.second.size() > 0) {
  //     auto& target = l->target(l->find(pair.first)); 
  //     for(auto& par : pair.second) {
  //       Message_Union_Find send_msg; 
  //       // send_msg.send_add_intermediate_root(par); // Add an intermediate root to the process of child
  //       send_msg.send_erase_intermediate_root(par); // Erase an intermediate root to the process of child

  //       cp.enqueue(target, send_msg);   
  //     }
  //     pair.second.clear(); 
  //   } 
  // }
  // // b->cache_send_intermediate_roots.clear(); 

}


void received_msg(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg, int from_gid) {
  // std::cout<<msg.tag<<std::endl; 

  // if(msg.tag == "gid_query") {
  //   answer_gid(b, cp, msg, from_gid); 
  // } else if (msg.tag == "gid_response") {
  //   save_gid(b, cp, msg); 
  // } 

  #if ISDEBUG
    std::cout<<"Start Receiving Msg: "<<msg.tag<<std::endl; 
  #endif

  if(msg.tag == "gparent") {
    distributed_save_gparent(b, cp, msg); 
  } else if(msg.tag == "child") {
    distributed_save_child(b, cp, msg); 
  } else if(msg.tag == "add_intermediate_root") {
    distributed_add_intermediate_root(b, cp, msg);
  } else if(msg.tag == "erase_intermediate_root") {
    distributed_erase_intermediate_root(b, cp, msg);
  } else if(msg.tag == "union") {
    distributed_save_union(b, cp, msg); 
  } else if(msg.tag == "gid_query") {
    answer_gid(b, cp, msg, from_gid); 
  } else if(msg.tag == "gid_response") {
    save_gid(b, cp, msg); 
  } else {
    std::cout<<"Error! "<<std::endl; 
  }
}

void receive_msg(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();
  
  while(!cp.empty_incoming_queues()) {
    std::vector<int> in;
    cp.incoming(in);

    for (unsigned i = 0; i < in.size(); ++i) {
      if(cp.incoming(in[i])) {
        Message_Union_Find msg; 
        cp.dequeue(in[i], msg);

        received_msg(b, cp, msg, in[i]); 
      }
    }
  }
}

bool union_find_iexchange(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  b->nchanges = 0; 
  
  receive_msg(b, cp); 
  local_computation(b, cp); 
  
  #if ISDEBUG
    int gid = cp.gid(); 
    std::cout<<"gid: "<<gid<<"================================="<<std::endl; 
  #endif

  #if OUTPUT_TIME_EACH_ROUND
    b->time = MPI_Wtime();
  #endif

  return b->nchanges == 0; 
}


void iexchange_process(diy::Master& master) {
  master.iexchange(&union_find_iexchange); 
}

void total_changes(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  cp.collectives()->clear();

  cp.all_reduce(b->nchanges, std::plus<int>()); 
  b->nchanges = 0;

  #if OUTPUT_TIME_EACH_ROUND
    b->time = MPI_Wtime();
  #endif
}

void exchange_process(diy::Master& master) {
  master.foreach(&receive_msg);
  master.foreach(&local_computation);
  master.foreach(&total_changes);
  master.exchange();

  #if ISDEBUG
    std::cout<<"================================="<<std::endl; 
  #endif
}

void exec_distributed_union_find(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<Block_Union_Find*>& blocks, bool is_iexchange = true, std::string filename_time_uf_w="") {

  std::vector<int> gids;                     // global ids of local blocks
  assigner.local_gids(world.rank(), gids);   // get the gids of local blocks for a given process rank 

  import_data(blocks, master, assigner, gids); 

  #if OUTPUT_TIME_EACH_ROUND
    #ifdef FTK_HAVE_MPI
      std::stringstream ss;
      ss << gids[0]; 

      MPI_Barrier(world); 
      double start = MPI_Wtime();
    #endif
  #endif

  if(is_iexchange) { // for iexchange
    master.foreach(&query_gid);
    iexchange_process(master);  

    // =========================================
    // Debug and Print

    #if OUTPUT_TIME_EACH_ROUND
      #ifdef FTK_HAVE_MPI
        if(!filename_time_uf_w.empty()) {
          double duration = blocks[0]->time - start; 
          ss << " " << duration ;

          MPI_Status status;
          MPI_File fh;

          ss<<std::endl;
          const std::string buf = ss.str();

          MPI_File_open(world, filename_time_uf_w.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
          
          MPI_File_write_ordered(fh, buf.c_str(), buf.length(), MPI_CHAR, &status);

          MPI_File_close(&fh);
        }
      #endif
    #endif

  } else { // for exchange
    bool all_done = false;

    master.foreach(&query_gid);
    master.exchange();

    #if OUTPUT_TIME_EACH_ROUND
      int round_cnt = 0; 
    #endif

    while(!all_done) {
      exchange_process(master); 
      int total_changes = master.proxy(master.loaded_block()).read<int>();
      // int total_changes = master.proxy(master.loaded_block()).get<int>();
      all_done = total_changes == 0;

      // =========================================
      // Debug and Print

      #if ISDEBUG
        std::cout<<total_changes<<std::endl; 
      #endif

      #if OUTPUT_TIME_EACH_ROUND
        #ifdef FTK_HAVE_MPI
          double duration = blocks[0]->time - start; 

          // ss << "Round "<<round_cnt << ":"; 
          ss << " " << duration ;

          round_cnt ++; 

          MPI_Barrier(world); 
          start = MPI_Wtime(); 
        #endif
      #endif
    }


    #if OUTPUT_TIME_EACH_ROUND
      #ifdef FTK_HAVE_MPI
        if(!filename_time_uf_w.empty()) {
          // std::string filename = filename_time_uf_w + std::to_string(world.rank()) + ".time"; 

          MPI_Status status;
          MPI_File fh;

          ss<<std::endl;
          const std::string buf = ss.str();

          MPI_File_open(world, filename_time_uf_w.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
          
          MPI_File_write_ordered(fh, buf.c_str(), buf.length(), MPI_CHAR, &status);

          MPI_File_close(&fh);
        }
      #endif
    #endif
  }

  // master.iexchange(&union_find_iexchange, 16, 1000);
}


// =================================================

int hash_string(std::string& ele, size_t nprocess) {
  // Algorithms come from: https://cp-algorithms.com/string/string-hashing.html
  
  const int p = 11; // since the characters are 0~9 plus ',' for separation
  
  int val = 0; 
  int p_pow = 1;
  for(char c : ele) {
    int c_val = c - '0'; // ASCII of comma ',' is 44, of '0' is 48
    if(c == ',') {
      c_val = 10; 
    }

    val = (val + c_val * p_pow) % nprocess; 
    p_pow = (p_pow * p) % nprocess; 
  }

  return val; 
}


// Gather all element-root information to the process of the root
  // Step one: send to target processes
  // Step two: gather on target processes

void send_2_target_processes(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();
  auto master = cp.master(); 
  auto& world = master->communicator(); 

  int nblocks = world.size();

  // auto eles = b->eles; 

  auto i = b->eles.begin();
  while(i != b->eles.end()) {
    auto& ele = (*i);

    std::string root = b->parent(ele); 
  
    int gid_root = hash_string(root, nblocks); 

    if(gid_root != gid) {
      std::pair<std::string, std::string> local_pair(ele, root); 

      auto& target = l->target(l->find(gid_root)); 

      Message_Union_Find send_msg; 
      send_msg.send_ele_parent_pair(local_pair); 

      cp.enqueue(target, send_msg); 

      i = b->eles.erase(i);
    } else {
      ++i ;
    }
    
  }
}

bool gather_on_target_processes(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {

  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  while(!cp.empty_incoming_queues()) {
    // Save unions from other blocks
    std::vector<int> in; // gids of incoming neighbors in the link
    cp.incoming(in);

    // for all neighbor blocks
    // dequeue data received from this neighbor block in the last exchange
    for (unsigned i = 0; i < in.size(); ++i) {
      if(cp.incoming(in[i])) {
        Message_Union_Find msg; 
        cp.dequeue(in[i], msg);

        if(msg.tag == "ele_parent_pair") {
            std::pair<std::string, std::string> pair; 
            msg.receive_ele_parent_pair(pair); 

            b->add(pair.first); 
            b->set_parent(pair.first, pair.second); 
        } else {
          std::cout<<"Wrong! Tag is not correct: "<<msg.tag<<std::endl; 
        }

      }
    }
  }

  return true;
}

// Get sets of elements by redistributing data
inline void Block_Union_Find::get_sets(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {
  #ifdef FTK_HAVE_MPI
    double start = MPI_Wtime();
  #endif

  master.foreach(&send_2_target_processes); 
  master.iexchange(&gather_on_target_processes); 

  Block_Union_Find* b = this;

  std::map<std::string, std::set<std::string>> root2set; 

  for(auto& ele : b->eles) {
    // if(!b->is_root(b->parent(ele))) {
    //   std::cout<<"Wrong! The parent is not root! get_sets_redistributed()"<< std::endl; 
    //   std::cout<<ele<<" "<<b->parent(ele)<<" "<<b->parent(b->parent(ele))<<std::endl; 
    // }
    assert(b->is_root(b->parent(ele))); 

    root2set[b->parent(ele)].insert(ele);  
  }

  for(auto ite = root2set.begin(); ite != root2set.end(); ++ite) {
    results.push_back(ite->second);   
  }

  #if ISDEBUG
    for(int i = 0; i < results.size(); ++i) {
      std::cout<<"Set "<<i<<":"<<std::endl; 
      auto& ele_set = results[i]; 
      for(auto& ele : ele_set) {
        std::cout<<ele<<" "; 
      }
      std::cout<<std::endl; 
    }
  #endif
}

#endif
