#ifndef _FTK_DISTRIBUTED_UNION_FIND_H
#define _FTK_DISTRIBUTED_UNION_FIND_H

#include <vector>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <iostream>

#include <ftk/ftk_config.hh>
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
#define OUTPUT_TIME_EACH_ROUND false

namespace ftk {
template <class IdType=std::string>
struct distributed_union_find
{
  distributed_union_find() : eles(), id2parent(), parent2oldchildren(), nonlocal_temporary_roots(), not_nonlocal_temporary_roots() {

  }

  // Initialization
  // Add and initialize elements
  void add(IdType i) {
    eles.insert(i); 
    id2parent.insert(std::make_pair(i, i)); 
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

  void add_nonlocal_temporary_root(IdType temporary_root) {
    if(not_nonlocal_temporary_roots.find(temporary_root) == not_nonlocal_temporary_roots.end()) {
      nonlocal_temporary_roots.insert(temporary_root); 
    }
  }

  void erase_nonlocal_temporary_root(IdType temporary_root) {
    nonlocal_temporary_roots.erase(temporary_root); 
    not_nonlocal_temporary_roots.insert(temporary_root); 
  }

  // Queries

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
  }

  IdType parent(IdType i) {
    assert(has(i));

    return id2parent[i]; 
  }

  const std::vector<IdType>& old_children(IdType i) {
    assert(has(i));

    return parent2oldchildren[i]; 
  }

  bool is_nonlocal_temporary_root(IdType i) {
    if(this->nonlocal_temporary_roots.find(i) != this->nonlocal_temporary_roots.end()) {
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
  std::map<IdType, std::vector<IdType>> parent2oldchildren; 

  std::set<IdType> nonlocal_temporary_roots;

  // Previous is nonlocal_temporary_root, but now they are not anymore
  std::set<IdType> not_nonlocal_temporary_roots;
};

}

// DIY Block for distributed union-find
struct Block_Union_Find : public ftk::distributed_union_find<std::string> {
  Block_Union_Find(): nchanges(0), related_elements(), all_related_elements(), temporary_root_2_gids(), nonlocal_temporary_roots_2_grandparents(), ele2gid(), distributed_union_find() { 
    
  }

  // add an element
  void add(std::string ele) {
    this->nchanges += 1;

    if(this->has(ele)) {
      std::cout<<"This ele has been added. "<<ele<<std::endl; 
    } else {
      distributed_union_find::add(ele); 
      this->related_elements.insert(std::make_pair(ele, std::set<std::string>())); 
      this->has_sent_gparent_query[ele] = false;
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

  void add_nonlocal_temporary_root(std::string temporary_root) {
    this->nchanges += 1;

    distributed_union_find::add_nonlocal_temporary_root(temporary_root); 
  }

  void erase_nonlocal_temporary_root(std::string temporary_root) {
    this->nchanges += 1;

    distributed_union_find::erase_nonlocal_temporary_root(temporary_root); 
  }

  bool is_temporary_root(std::string i) {
    if(this->has(i)) {
      // if(this->is_root(i) && this->get_related_elements(i).size() == 0) {
      if(this->is_root(i)) {
        return true ;
      } else {
        return false ;
      }
    } else {
      return distributed_union_find::is_nonlocal_temporary_root(i); 
    }

    return false ;
  }


  void get_sets(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results);

public: 

  std::map<std::string, int> ele2gid; 

  int nchanges = 0; // # of processed unions per round = valid unions (united unions) + passing unions

  #if OUTPUT_TIME_EACH_ROUND
    double time = 0; // time for measure duration of each round;
    double time_final_update_start;
  #endif
  
  // std::map<int, std::set<std::string>> cache_send_temporary_roots ;

  std::map<std::string, std::set<int>> temporary_root_2_gids; 

  std::map<std::string, std::string> nonlocal_temporary_roots_2_grandparents; 

  std::map<std::string, bool> has_sent_gparent_query; 

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

  void send_gparent_query(const std::string& child, const std::string& parent, const int& gid_child) {
    tag = "gparent_query";

    strs.push_back(child); 
    strs.push_back(parent); 
    strs.push_back(std::to_string(gid_child)); 
  }

  void send_gparent(std::string& ele, std::string& gparent, int& gid_gparent, bool& is_known) {
    tag = "gparent";

    strs.push_back(ele); 
    strs.push_back(gparent); 
    strs.push_back(std::to_string(gid_gparent)); 
    strs.push_back(std::to_string(is_known)); 
  }

  void send_add_temporary_root(const std::string& parent) {
    tag = "add_temporary_root";

    strs.push_back(parent); 
  }

  void send_erase_temporary_root(const std::string& parent, const std::string& grandparent, const int& gid_grandparent, const bool& is_known ) {
    tag = "erase_temporary_root";

    strs.push_back(parent); 
    strs.push_back(grandparent); 
    strs.push_back(std::to_string(gid_grandparent)); 
    strs.push_back(std::to_string(is_known)); 
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

  void rec_gparent_query(std::string& child, std::string& parent, int& gid_child) {
    child = strs[0]; 
    parent = strs[1]; 
    gid_child = std::stoi(strs[2]);  
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

  void rec_add_temporary_root(std::string& parent) {
    parent = strs[0];
  }

  void rec_erase_temporary_root(std::string& parent, std::string& grandparent, int& gid_grandparent, bool& is_known) {
    parent = strs[0];

    grandparent = strs[1]; 
    gid_grandparent = std::stoi(strs[2]); 
    is_known = std::stoi(strs[3]); 
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


void send_erase_temporary_root_to_all_processes(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, const std::string& parent, const std::string& grandparent, const int& gid_grandparent, const bool& is_known) {
    
    if(b->temporary_root_2_gids.find(parent) != b->temporary_root_2_gids.end()) {
      diy::Link* l = cp.link();
      for(auto& gid : b->temporary_root_2_gids[parent]) {
        Message_Union_Find send_msg; 
        send_msg.send_erase_temporary_root(parent, grandparent, gid_grandparent, is_known); // Erase an temporary root to the process of child
        cp.enqueue(l->target(l->find(gid)), send_msg); 

        if(is_known) {
          b->temporary_root_2_gids[grandparent].insert(gid); 
        }

      }
      b->temporary_root_2_gids.erase(parent); 
    }

}


void unite_once(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  for(auto& ele : b->eles) {

    if(b->is_root(ele)) {

      // If having an assign_global_roots
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

      // ================================================================

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


      // ================================================================

      if(found_related_ele < ele) {

        int rgid = b->get_gid(found_related_ele); 
        // std::cout<<gid<<" "<<rgid<<std::endl; 
        
        b->set_parent(ele, found_related_ele); 

        #if ISDEBUG
          std::cout<<ele<<" -1> "<<found_related_ele<<std::endl; 
        #endif

        b->erase_related_element(ele, found_related_ele); 

        bool is_known = b->is_temporary_root(found_related_ele);
        send_erase_temporary_root_to_all_processes(b, cp, ele, found_related_ele, rgid, is_known);
      } 

    }
  }

}


void compress_path(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  if(b->nonlocal_temporary_roots_2_grandparents.size() > 0) {
    for(auto& ele : b->eles) {
      if(!b->is_root(ele)) {
        std::string parent = b->parent(ele);
        if(b->nonlocal_temporary_roots_2_grandparents.find(parent) != b->nonlocal_temporary_roots_2_grandparents.end()) {
          std::string grandparent = b->nonlocal_temporary_roots_2_grandparents[parent]; 
          b->set_parent(ele, grandparent); 
        }
      }
    }
    b->nonlocal_temporary_roots_2_grandparents.clear();
  }

  for(auto& ele : b->eles) {

    if(!b->is_root(ele)) {
      std::string parent = b->parent(ele);

      if(b->has(parent)) {
        
        if(!b->is_root(parent)) {
          std::string grandparent = b->parent(parent);

          if(b->has(grandparent)) { // If having an assign_global_roots
            b->set_parent(ele, grandparent);    
          }

          // if(b->has(grandparent) || b->is_temporary_root(grandparent)) { // Hollow tree structure
          //   b->set_parent(ele, grandparent);    
          // }

          // b->set_parent(ele, grandparent); // Freely

        }

      } else {
        
        bool is_known = b->is_temporary_root(parent); 
        if(!is_known) {
          if(b->has_sent_gparent_query[ele]) {
            continue ;
          }

          if(b->has_gid(parent)) {
            int gid_parent = b->get_gid(parent);

            Message_Union_Find send_msg; 
            
            send_msg.send_gparent_query(ele, parent, gid); 
            cp.enqueue(l->target(l->find(gid_parent)), send_msg); 

            b->has_sent_gparent_query[ele] = true;
          }
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

  b->has_sent_gparent_query[ele] = false;

  b->set_parent(ele, grandpar); 
  b->set_gid(grandpar, gid_grandparent); 
  if(is_known && !b->has(grandpar)) {
    b->add_nonlocal_temporary_root(grandpar);
  }

  #if ISDEBUG
    std::cout<<ele<<" -2> "<<grandpar<<std::endl; 
  #endif
}

// Distributed path compression - Step Three
void distributed_answer_gparent_query(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  diy::Link* l = cp.link();

  std::string child; 
  std::string parent; 
  int gid_child; 

  msg.rec_gparent_query(child, parent, gid_child); 

  b->set_gid(child, gid_child);

  std::string grandparent = b->parent(parent); 
  int gid_grandparent = b->get_gid(grandparent);

  // Set child's parent to grandparent
  Message_Union_Find send_msg; 
  bool is_known = b->is_temporary_root(grandparent); 

  if(is_known) {
    b->temporary_root_2_gids[grandparent].insert(gid_child);
  }

  send_msg.send_gparent(child, grandparent, gid_grandparent, is_known); 
  cp.enqueue(l->target(l->find(gid_child)), send_msg); 
}

// Receive adding an temporary root
void distributed_add_temporary_root(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  std::string par;
  msg.rec_add_temporary_root(par); // Add an temporary root

  b->add_nonlocal_temporary_root(par); 
}

void distributed_erase_temporary_root (Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  std::string parent;
  std::string grandparent;
  int gid_grandparent; 
  bool is_known;
  msg.rec_erase_temporary_root(parent, grandparent, gid_grandparent, is_known); // Erase an temporary root, and if some elements' parents are this temporary root before, we replace with the grandparent

  b->set_gid(grandparent, gid_grandparent); 
  b->nonlocal_temporary_roots_2_grandparents[parent] = grandparent;
  b->erase_nonlocal_temporary_root(parent); 

  if(is_known) {
    b->add_nonlocal_temporary_root(grandparent); 
  }

  is_known = is_known || b->is_temporary_root(grandparent); 
  send_erase_temporary_root_to_all_processes(b, cp, parent, grandparent, gid_grandparent, is_known);
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

        // if(!b->is_temporary_root(par)) { // Previously test best
        //   continue ;
        // }
        // if(!is_local_parent && !b->is_temporary_root(par)) { 
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

        if(!(b->is_temporary_root(par) || (is_local_parent && !b->has(b->parent(par))))) {  // If having an assign_global_roots
          continue ;
        }

        std::vector<std::string> cache;
        for(auto& related_ele : src) {

          if(related_ele < par) {

            if(is_local_parent) {
              // send_erase_temporary_root_to_all_processes(b, cp, par);

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
              // send_erase_temporary_root_to_all_processes(b, cp, related_ele);

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
  diy::Link* l = cp.link();

  // std::cout<<"Save Unions: "<<std::endl; 
  // std::cout<<"Block ID: "<<cp.gid()<<std::endl; 

  std::string ele; 
  std::string related_ele; 
  int rgid; 
  
  msg.rec_union(ele, related_ele, rgid); 

  // std::cout<<*ele_ptr<<" - "<<*par_ptr << " - " << *related_ele_ptr <<" - "<< *gid_ptr<<std::endl; 

  // send_erase_temporary_root_to_all_processes(b, cp, ele); 

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


  // if(b->cache_send_temporary_roots.size() > 0) {
  //   diy::Link* l = cp.link();
  //   for(auto& pair : b->cache_send_temporary_roots) {
  //     if(pair.second.size() > 0) {
  //       auto& target = l->target(l->find(pair.first)); 
  //       for(auto& par : pair.second) {
  //         Message_Union_Find send_msg; 
  //         send_msg.send_add_temporary_root(par); // Add an temporary root to the process of child

  //         cp.enqueue(target, send_msg);   
  //       }
  //       pair.second.clear(); 
  //     } 
  //   }
  //   b->cache_send_temporary_roots.clear(); 
  // }

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
  } else if(msg.tag == "gparent_query") {
    distributed_answer_gparent_query(b, cp, msg); 
  } else if(msg.tag == "add_temporary_root") {
    distributed_add_temporary_root(b, cp, msg);
  } else if(msg.tag == "erase_temporary_root") {
    distributed_erase_temporary_root(b, cp, msg);
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

void assign_global_roots(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  #if OUTPUT_TIME_EACH_ROUND
    b->time_final_update_start = MPI_Wtime();
  #endif
  for(auto& ele : b->eles) {

    if(!b->is_root(ele)) {
      std::string parent = b->parent(ele);

      if(b->has(parent)) {
        if(!b->is_root(parent)) {
          std::string grandparent = b->parent(parent);

          b->set_parent(ele, grandparent);
        }
      }
    }
  }
}


void iexchange_process(diy::Master& master) {
  master.iexchange(&union_find_iexchange); 
  master.foreach(&assign_global_roots); 
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
  
  // master.exchange();

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
          double time_final_update_end = MPI_Wtime(); 

          double duration = blocks[0]->time - start; 
          ss << " " << duration ;

          // double duration_final_update = time_final_update_end - blocks[0]->time_final_update_start; 
          // ss << " " << duration_final_update ;

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

      master.exchange();
      int total_changes = master.proxy(master.loaded_block()).read<int>();
      // int total_changes = master.proxy(master.loaded_block()).get<int>();
      all_done = total_changes == 0;

      // =========================================
      // Debug and Print

      #if ISDEBUG
        std::cout<<total_changes<<std::endl; 
      #endif

    }

    #if OUTPUT_TIME_EACH_ROUND
      #ifdef FTK_HAVE_MPI
        double duration = MPI_Wtime() - start; 

        // ss << "Round "<<round_cnt << ":"; 
        ss << " " << duration ;

        round_cnt ++; 

        MPI_Barrier(world); 
        start = MPI_Wtime(); 
      #endif
    #endif

    // master.foreach(&assign_global_roots); 
    // #if OUTPUT_TIME_EACH_ROUND
    //   #ifdef FTK_HAVE_MPI
    //     double time_final_update_end = MPI_Wtime(); 
    //     double duration_final_update = time_final_update_end - blocks[0]->time_final_update_start; 

    //     // ss << "Round "<<round_cnt << ":"; 
    //     ss << " " << duration_final_update ;

    //     round_cnt ++; 

    //     MPI_Barrier(world); 
    //     start = MPI_Wtime(); 
    //   #endif
    // #endif


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
