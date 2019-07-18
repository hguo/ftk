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

#define IEXCHANGE 1
#define ISDEBUG   0

namespace ftk {
template <class IdType=std::string>
struct distributed_union_find
{
  distributed_union_find() : eles(), id2parent(), parent2children() {

  }

  // Initialization
  // Add and initialize elements
  void add(IdType i) {
    eles.insert(i); 
    id2parent.insert(std::make_pair(i, i)); 
    parent2children.insert(std::make_pair(i, std::vector<IdType>())); 
  }

  // Operations

  void set_parent(IdType i, IdType par) {
    if(!has(i)) {
      std::cout<< "No element set_parent(). " <<std::endl;
      exit(0); 
    }

    // Ensure the id of new parent is smaller than the id of old parent
      // Otherwise, may occur override and erase the updated parent
    if(par < id2parent[i]) {
      id2parent[i] = par; 
    }
    // else {
    //   std::cout<< "Parent is override! " <<std::endl;
    // }
  }

  void add_child(IdType par, IdType ele) {
    if(!has(par)) {
      std::cout<< "No element parent add_child(). " <<std::endl;
      exit(0); 
    }

    parent2children[par].push_back(ele); 
  }

  // Queries

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
  }

  IdType parent(IdType i) {
    if(!has(i)) {
      std::cout<< "No element parent parent(). " <<std::endl;
      exit(0); 
    }

    return id2parent[i]; 
  }

  const std::vector<IdType>& children(IdType i) {
    if(!has(i)) {
      std::cout<< "No element child children(). " <<std::endl;
      exit(0); 
    }

    return parent2children[i]; 
  }

  void clear_children(std::string i) {
    if(!has(i)) {
      std::cout<< "No element child children(). " <<std::endl;
      exit(0); 
    }

    parent2children[i].clear(); 
  }

  bool is_root(IdType i) {
    if(!has(i)) {
      std::cout<< "No element is_root(). " <<i<<std::endl;
      exit(0); 
    }

    return i == id2parent[i]; 
  }

public:
  std::set<IdType> eles; 

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id2parent;
  std::map<IdType, std::vector<IdType>> parent2children;
};

}

// ==========================================================

struct intersection_t {
  float x[3]; // the spacetime coordinates of the trajectory
  float val; // scalar value at the intersection
  
  std::string eid; // the "size_t" type is not scalable for element id, for example, when number of elements is large

  template <class Archive> void serialize(Archive & ar) {
    ar(eid, x[0], x[1], x[2], val);
  }
};

// ==========================================================

struct Block : public ftk::distributed_union_find<std::string> {
  Block(): nchanges(0), related_elements(), ele2gid(), intersections(), distributed_union_find() { 
    
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
    if(this->has(ele)) {
      if(this->related_elements[ele].find(related_ele) == this->related_elements[ele].end()) {
        return false; 
      }

      return true; 
    } else {
      std::cout<<"Don't have this element "<<ele<<std::endl; 
      exit(0); 
    }
  }

  void add_related_element(std::string ele, std::string related_ele) {
    this->nchanges += 1;

    if(this->has(ele)) {
      if(ele != related_ele) {
        this->related_elements[ele].insert(related_ele);   
      }
    } else {
      std::cout<<"Don't have this element "<<ele<<std::endl; 
      exit(0); 
    }
  }

  const std::set<std::string>& get_related_elements(std::string ele) {
    if(this->has(ele)) {
      return this->related_elements[ele]; 
    } else {
      std::cout<<"Don't have this element "<<ele<<std::endl; 
      exit(0); 
    }
  }

  void clear_related_elements(std::string ele) {
    this->nchanges += 1;

    if(this->has(ele)) {
      this->related_elements[ele].clear(); 
    } else {
      std::cout<<"Don't have this element "<<ele<<std::endl; 
      exit(0); 
    }
  }

  void erase_related_element(std::string ele, std::string related_element) {
    this->nchanges += 1;

    if(this->has(ele)) {
      this->related_elements[ele].erase(related_element); 
    } else {
      std::cout<<"Don't have this element "<<ele<<std::endl; 
      exit(0); 
    }
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
    if(this->has_gid(ele)) {
      return this->ele2gid[ele]; 
    } else {
      std::cout<<"Don't have gid of "<<ele<<std::endl; 
      exit(0); 

      return -1; 
    }
  }

  void set_parent(std::string i, std::string par) {
    this->nchanges += 1;

    distributed_union_find::set_parent(i, par); 
  }

  void add_child(std::string par, std::string ele) {
    this->nchanges += 1;

    distributed_union_find::add_child(par, ele); 
  }

  void clear_children(std::string par) {
    this->nchanges += 1;

    distributed_union_find::clear_children(par); 
  }

public: 

  std::map<std::string, intersection_t> intersections; 
  std::map<std::string, int> ele2gid; 

  int nchanges; // # of processed unions per round = valid unions (united unions) + passing unions
  
private:
  // map element id to ids of its related elements
  // Can be optimized by ordered the related elements, put related elements on process first and ordered decreingly by ids
  std::map<std::string, std::set<std::string>> related_elements; 
};

// ==========================================================

struct Message {

  Message() : tag(), strs() {

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

  void send_gparent(std::string& ele, std::string& gparent, int& gid_gparent) {
    tag = "gparent";

    strs.push_back(ele); 
    strs.push_back(gparent); 
    strs.push_back(std::to_string(gid_gparent)); 
  }

  void send_child(const std::string& parent, const std::string& child, const int& gid_child) {
    tag = "child";

    strs.push_back(parent); 
    strs.push_back(child); 
    strs.push_back(std::to_string(gid_child)); 
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

  void send_intersection(const std::pair<std::string, intersection_t>& pair) {
    tag = "intersection"; 

    strs.push_back(pair.first); 

    auto& I = pair.second; 
    strs.push_back(std::to_string(I.x[0])); 
    strs.push_back(std::to_string(I.x[1])); 
    strs.push_back(std::to_string(I.x[2])); 
    strs.push_back(std::to_string(I.val)); 
    strs.push_back(I.eid); 
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

  void rec_gparent(std::string& ele, std::string& gparent, int& gid_gparent) {
    ele = strs[0]; 
    gparent = strs[1];

    gid_gparent = std::stoi(strs[2]); 
  }

  void rec_child(std::string& parent, std::string& child, int& gid_child) {
    parent = strs[0]; 
    child = strs[1];

    gid_child = std::stoi(strs[2]); 
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

  void receive_intersection(std::pair<std::string, intersection_t>& pair) {
    pair.first = strs[0]; 

    intersection_t& I = pair.second; 
    I.x[0] = std::stof(strs[1]); 
    I.x[1] = std::stof(strs[2]); 
    I.x[2] = std::stof(strs[3]); 
    I.val = std::stof(strs[4]); 
    I.eid = strs[5];
  }

public:
  std::string tag; 

  std::vector<std::string> strs; 
}; 

// ==========================================================

namespace diy
{
    template<>
        struct Serialization<Message>
    {
        static void save(BinaryBuffer& bb, const Message& msg)
        {
            diy::save(bb, msg.tag);
            // diy::save(bb, msg.str);
            diy::save(bb, msg.strs);
        }

        static void load(BinaryBuffer& bb, Message& msg)
        {
            diy::load(bb, msg.tag);
            // diy::load(bb, msg.str);
            diy::load(bb, msg.strs);
        }
    };
}

// ==========================================================

void query_gid(Block* b, const diy::Master::ProxyWithLink& cp) {
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
          Message msg = Message(); 
          msg.send_gid_query(ele); 

          cp.enqueue(l->target(i), msg);
        }
      }

    }
  }

  // std::cout<<std::endl; 
}


void answer_gid(Block* b, const diy::Master::ProxyWithLink& cp, Message msg, int from_gid) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();
  // diy::Master* master = cp.master(); 

  // std::cout<<"Answer: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 

  std::string ele; 
  msg.rec_gid_query(ele); 

  if(b->has(ele)) {
    Message send_msg; 
    send_msg.send_gid_response(ele, gid); 

    // std::cout<<ele<<"-";
    // std::cout<<gid;
    // std::cout<<std::endl; 

    cp.enqueue(l->target(l->find(from_gid)), send_msg);
    // cp.enqueue(l->target(l->find(in[i])), 1);
  }

  // std::cout<<std::endl; 
}

void save_gid(Block* b, const diy::Master::ProxyWithLink& cp, Message msg) {
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


void unite_once(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  for(auto& ele : b->eles) {

    // if(ele == "7") {
    //   std::cout<<ele<<": "<<b->related_elements[ele].size()<<std::endl; 
    //   for(std::vector<std::string>::iterator it_vec = b->related_elements[ele].begin(); it_vec != b->related_elements[ele].end(); ++it_vec) {
    //     std::string related_ele = *it_vec; 

    //     std::cout<<ele<<"~~~~"<<related_ele<<std::endl; 
    //   }
    // }

    if(b->is_root(ele)) {
      auto related_elements = b->get_related_elements(ele); // cannot use auto&, since we will remove some elements from the original set; auto& will cause segmentation erro

      for(auto& related_ele : related_elements) {
        
        if(!b->has_gid(related_ele)) {
          continue ; 
        }

        int rgid = b->get_gid(related_ele); 
        // std::cout<<gid<<" "<<rgid<<std::endl; 

        // Unite a root element with larger id or smaller id is better? 
          // Here is a smaller id

        if(related_ele < ele) {
          b->set_parent(ele, related_ele); 

          if(b->has(related_ele)) {
            b->add_child(related_ele, ele); 
          } else {
            Message send_msg; 
            send_msg.send_child(related_ele, ele, gid); 

            cp.enqueue(l->target(l->find(rgid)), send_msg);
          }

          if(ISDEBUG) {
            std::cout<<ele<<" -1> "<<related_ele<<std::endl; 
          }

          b->erase_related_element(ele, related_ele); 

          break ; 
        } else {

          if(b->has(related_ele)) {
            b->add_related_element(related_ele, ele); 
          } else {
            Message send_msg; 
            send_msg.send_union(related_ele, ele, gid); 

            cp.enqueue(l->target(l->find(rgid)), send_msg);
          }

          b->erase_related_element(ele, related_ele); 
        }

      }
    }
  }

}


// Compress path - Step One
void compress_path(Block* b, const diy::Master::ProxyWithLink& cp) {
  diy::Link* l = cp.link();

  for(auto& parent : b->eles) {
    if(!b->is_root(parent)) {
      auto& children = b->children(parent); 

      if(children.size() > 0) {
        std::string grandparent = b->parent(parent); 
        bool is_local_grandparent = true; 
        int gid_grandparent = b->get_gid(grandparent); 
        
        if(!b->has(grandparent)) { // if the grandparent is not in this block
          is_local_grandparent = false ;
          if(gid_grandparent == -1) { // if currently the gid of grandparent is not available
            continue ;
          }
        }

        std::vector<std::string> cache; 
        for(auto child : children) {
          
          // Set child's parent to grandparent
          if(b->has(child)) {
            b->set_parent(child, grandparent); 
          } else { // if the child is not in this block
            if(!b->has_gid(child)) { // if currently the gid of child is not available
              cache.push_back(child); 
              continue ;
            }

            int gid_child = b->get_gid(child); 

            Message send_msg; 
            send_msg.send_gparent(child, grandparent, gid_grandparent); 

            // std::cout<<*ele_ptr<<" - "<<*parent_ptr<<" - "<<grandparent<<" - "<<b->ele2gid[grandparent]<<std::endl; 

            cp.enqueue(l->target(l->find(gid_child)), send_msg); 
          }

          // Add child to grandparent's child list
          if(is_local_grandparent) {
            b->add_child(grandparent, child); 
          } else {
            int gid_child = b->get_gid(child); 

            Message send_msg; 
            send_msg.send_child(grandparent, child, gid_child); 

            cp.enqueue(l->target(l->find(gid_grandparent)), send_msg); 
          }
        }

        b->clear_children(parent); 

        if(cache.size() > 0) {
          for(auto& cache_ele : cache) {
            b->add_child(parent, cache_ele); 
          }
        }
      }
    }
  }
}


// Distributed path compression - Step Two
void distributed_save_gparent(Block* b, const diy::Master::ProxyWithLink& cp, Message msg) {
  std::string ele; 
  std::string grandpar; 
  int gid_grandparent; 

  msg.rec_gparent(ele, grandpar, gid_grandparent); 
  
  // std::cout<<*ele_ptr << " - " << b->parent(*ele_ptr) <<" - "<< *grandpar_ptr<<" - "<<*grandpar_gid_ptr<<std::endl; 

  // if(*ele_ptr == "4") {
  //   std::cout<<*ele_ptr<<" ~ "<<*grandpar_ptr<<std::endl; 
  // }


  b->set_parent(ele, grandpar); 
  b->set_gid(grandpar, gid_grandparent); 

  if(ISDEBUG) {
    std::cout<<ele<<" -2> "<<grandpar<<std::endl; 
  }
}

// Distributed path compression - Step Three
void distributed_save_child(Block* b, const diy::Master::ProxyWithLink& cp, Message msg) {
  std::string par; 
  std::string child; 
  int gid_child; 

  msg.rec_child(par, child, gid_child); 

  b->add_child(par, child); 
  b->set_gid(child, gid_child); 
}


void pass_unions(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  // Local computation
  // Pass unions of elements in this block to their parents, save these unions
  for(auto& ele : b->eles) {

    if(!b->is_root(ele)) {

      auto& src = b->get_related_elements(ele);
      std::vector<std::string> cache;

      if(src.size() > 0) {
        std::string par = b->parent(ele); 

        if(b->has(par)) {
          // Update directly, if the realted element is in the block

          for(auto& related_ele : src) {
            b->add_related_element(par, related_ele); 
          } 
        } else {
          // Communicate with other processes
          int pgid = b->get_gid(par); 

          for(auto& related_ele : src) {
            if(!b->has_gid(related_ele)) {
              cache.push_back(related_ele); 
              continue ;
            }
            int r_gid = b->get_gid(related_ele); 

            Message send_msg; 
            send_msg.send_union(par, related_ele, r_gid); 

            cp.enqueue(l->target(l->find(pgid)), send_msg); 
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
void distributed_save_union(Block* b, const diy::Master::ProxyWithLink& cp, Message msg) {
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

void local_computation(Block* b, const diy::Master::ProxyWithLink& cp) {

  if(ISDEBUG) {
    int gid = cp.gid();
    std::cout<<"Start Local Computation. "<<"gid: "<<gid<<std::endl; 
  }

  unite_once(b, cp); 

  if(ISDEBUG) {
    int gid = cp.gid();
    std::cout<<"Finish unite once. "<<"gid: "<<gid<<std::endl; 
  }

  compress_path(b, cp); 

  if(ISDEBUG) {
    int gid = cp.gid();
    std::cout<<"Finish compress_path. "<<"gid: "<<gid<<std::endl; 
  }

  pass_unions(b, cp); 

  if(ISDEBUG) {
    int gid = cp.gid();
    std::cout<<"Finish pass_unions. "<<"gid: "<<gid<<std::endl; 
  }
}


void received_msg(Block* b, const diy::Master::ProxyWithLink& cp, Message msg, int from_gid) {
  // std::cout<<msg.tag<<std::endl; 

  // if(msg.tag == "gid_query") {
  //   answer_gid(b, cp, msg, from_gid); 
  // } else if (msg.tag == "gid_response") {
  //   save_gid(b, cp, msg); 
  // } 

  if(ISDEBUG) {
    std::cout<<"Start Receiving Msg: "<<msg.tag<<std::endl; 
  }

  if(msg.tag == "gparent") {
    distributed_save_gparent(b, cp, msg); 
  } else if(msg.tag == "child") {
    distributed_save_child(b, cp, msg); 
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

void receive_msg(Block* b, const diy::Master::ProxyWithLink& cp) {
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
        Message msg; 
        cp.dequeue(in[i], msg);

        received_msg(b, cp, msg, in[i]); 
      }
    }
  }
}

void total_changes(Block* b, const diy::Master::ProxyWithLink& cp) {
  cp.collectives()->clear();

  cp.all_reduce(b->nchanges, std::plus<int>()); 
  b->nchanges = 0;
}

bool union_find_iexchange(Block* b, const diy::Master::ProxyWithLink& cp) {
  b->nchanges = 0; 

  int gid = cp.gid();

  // std::cout<<"gid: "<<gid<<"=====Phase 1============================"<<std::endl; 

  receive_msg(b, cp); 

  // std::cout<<"gid: "<<gid<<"=====Phase 2============================"<<std::endl; 

  local_computation(b, cp); 

  // std::cout<<"gid: "<<gid<<"=====Phase 3============================"<<std::endl; 
  
  if(ISDEBUG) {
    int gid = cp.gid(); 

    std::cout<<"gid: "<<gid<<"================================="<<std::endl; 
  }

  return b->nchanges == 0; 
}

void iexchange_process(diy::Master& master) {
  master.iexchange(&union_find_iexchange); 

  // Test whether the iexchange ends but the union find does not end
  // master.iexchange(&union_find_iexchange); 
}

void union_find_exchange(Block* b, const diy::Master::ProxyWithLink& cp) {
  compress_path(b, cp); 
  pass_unions(b, cp); 
  
  receive_msg(b, cp); 

  total_changes(b, cp); 
}

// Method 1
// gparent query and gparent answer should be synchronized. 
  // otherwise it is possible that one query message is transferring, but the program ends. 
  // Other parts are not required to synchronize. 
void exchange_process(diy::Master& master) {
  master.foreach(&unite_once);
  master.foreach(&compress_path);
  master.exchange();                 
  master.foreach(&receive_msg); // distributed_answer_gparent + additional responses

  master.exchange();
  master.foreach(&receive_msg); // distributed_save_gparent + additional responses

  // master.foreach(&union_find_exchange);
  // master.exchange();

  master.foreach(&compress_path);
  master.foreach(&pass_unions);
  master.foreach(&pass_unions);
  master.foreach(&total_changes);
  master.exchange();

  // While two consecutive exchange() without receiving will clear the buffer
    // Which means, in-between each pair of exchange should contain at least one receive
    // Hence, an additional receive is needed
  master.foreach(&receive_msg);

  if(ISDEBUG) {
    std::cout<<"================================="<<std::endl; 
  }
}

void import_data(std::vector<Block*>& blocks, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<int>& gids) {
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

    Block* b = blocks[i]; 
    master.add(gid, b, link); 
  }
}

void exec_distributed_union_find(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<Block*>& blocks) {

  std::vector<int> gids;                     // global ids of local blocks
  assigner.local_gids(world.rank(), gids);   // get the gids of local blocks for a given process rank 

  import_data(blocks, master, assigner, gids); 
  
  if(IEXCHANGE) { // for iexchange
    // #ifndef DIY_NO_MPI
    // #ifdef FTK_HAVE_MPI
    //   // MPI_Barrier(world); 
    //   double start = MPI_Wtime();
    // #endif

    master.foreach(&query_gid);
    iexchange_process(master);  

    // #ifdef FTK_HAVE_MPI
    //   double end = MPI_Wtime();
    //   if(world.rank() == 0) {
    //     std::cout << "Union-Find: " << end - start << " seconds. " << std::endl;
    //   }
    // #endif

  } else { // for exchange
    bool all_done = false;
    while(!all_done) {
      exchange_process(master); 

      int total_changes = master.proxy(master.loaded_block()).read<int>();
      // int total_changes = master.proxy(master.loaded_block()).get<int>();

      if(ISDEBUG) {
        std::cout<<total_changes<<std::endl; 
      }

      all_done = total_changes == 0;
    }
  }

  // master.iexchange(&union_find_iexchange, 16, 1000);
}




// Generate sets of elements

// Method 1:
  // Send to p0
void send_2_p0(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  if(gid != 0) {
    // std::vector<std::pair<std::string, std::string>> local_pairs; 

    auto& target = l->target(l->find(0));
    for(auto& ele : b->eles) {
      std::string parent = b->parent(ele); 

      // local_pairs.push_back(std::make_pair(ele, parent)); 

      std::pair<std::string, std::string> local_pair(ele, parent); 

      Message send_msg; 
      send_msg.send_ele_parent_pair(local_pair); 

      cp.enqueue(target, send_msg); 
    }

    for(auto& pair : b->intersections) {
      Message send_msg; 
      send_msg.send_intersection(pair); 

      cp.enqueue(target, send_msg); 
    }
  }
}


// Gather all element-parent information to the first block
bool gather_2_p0(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  if(gid == 0) {
    while(!cp.empty_incoming_queues()) {
      // Save unions from other blocks
      std::vector<int> in; // gids of incoming neighbors in the link
      cp.incoming(in);

      // for all neighbor blocks
      // dequeue data received from this neighbor block in the last exchange
      for (unsigned i = 0; i < in.size(); ++i) {
        if(cp.incoming(in[i])) {
          Message msg; 
          cp.dequeue(in[i], msg);

          if(msg.tag == "ele_parent_pair") {
            std::pair<std::string, std::string> pair; 
            msg.receive_ele_parent_pair(pair); 

            b->add(pair.first); 
            b->set_parent(pair.first, pair.second); 
          } else if(msg.tag == "intersection") {
            std::pair<std::string, intersection_t> pair; 
            msg.receive_intersection(pair);

            b->intersections.insert(pair); 
          } else {
            std::cout<<"Wrong! Tag is not correct: "<<msg.tag<<std::endl; 
          }

        }
      }
    }
  }

  return true;
}


// Get sets of elements
void get_sets_on_p0(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {
  master.foreach(&send_2_p0); 
  master.iexchange(&gather_2_p0); 

  if(world.rank() == 0) {
    Block* b = static_cast<Block*> (master.get(0)); 

    std::map<std::string, std::set<std::string>> root2set; 

    #ifdef FTK_HAVE_MPI
      std::cout<<"# of elements on proc. " << world.rank() <<" : "<< b->eles.size()<<std::endl; 
    #endif

    for(auto& ele : b->eles) {
      if(!b->is_root(b->parent(ele))) {
        std::cout<<"Wrong! The parent is not root! "<< std::endl; 
        std::cout<<ele<<" "<<b->parent(ele)<<" "<<b->parent(b->parent(ele))<<std::endl; 
      }

      root2set[b->parent(ele)].insert(ele);  
    }

    for(auto ite = root2set.begin(); ite != root2set.end(); ++ite) {
      results.push_back(ite->second);   
    }

    if(ISDEBUG) {
      for(int i = 0; i < results.size(); ++i) {
        std::cout<<"Set "<<i<<":"<<std::endl; 
        std::set<std::string>& ele_set = results[i]; 
        for(auto& ele : ele_set) {
          std::cout<<ele<<" "; 
        }
        std::cout<<std::endl; 
      }
    }
  }
}


// Method 2:
// Gather all element-root information to the process of the root
  // Since the roots have all children, can directly use the information to reconstruct the sets; the only lacking part is the intersection. 
  // Step one: send to root
  // Step two: gather on root
void send_2_roots(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  for(auto& ele : b->eles) {
    std::string root = b->parent(ele); 
    if(b->has(root)) continue ; 

    int gid_root = b->get_gid(root); 
    auto& target = l->target(l->find(gid_root)); 

    Message send_msg_intersection; 
    send_msg_intersection.send_intersection(*(b->intersections.find(ele))); 

    cp.enqueue(target, send_msg_intersection); 
  }
}

bool gather_on_roots(Block* b, const diy::Master::ProxyWithLink& cp) {

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
        Message msg; 
        cp.dequeue(in[i], msg);

        if(msg.tag == "intersection") {
          std::pair<std::string, intersection_t> pair; 
          msg.receive_intersection(pair);

          b->intersections.insert(pair); 
        } else {
          std::cout<<"Wrong! Tag is not correct: "<<msg.tag<<std::endl; 
        }

      }
    }
  }

  return true;
}

// Get sets of elements on processes of roots
void get_sets_on_roots(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {
  master.foreach(&send_2_roots); 
  master.iexchange(&gather_on_roots); 

  // std::vector<int> gids;                     // global ids of local blocks
  // assigner.local_gids(world.rank(), gids);

  Block* b = static_cast<Block*> (master.get(0)); // load block with local id 0

  std::map<std::string, std::set<std::string>> root2set; 

  for(auto& ele : b->eles) {
    if(b->is_root(ele)) {
      root2set[ele].insert(ele);  
      auto& children = b->children(ele); 
      for(auto& child : children) {
        root2set[ele].insert(child); 
      }
    }
  }

  for(auto ite = root2set.begin(); ite != root2set.end(); ++ite) {
    results.push_back(ite->second);   
  }

  if(ISDEBUG) {
    for(int i = 0; i < results.size(); ++i) {
      std::cout<<"Set "<<i<<":"<<std::endl; 
      std::set<std::string>& ele_set = results[i]; 
      for(auto& ele : ele_set) {
        std::cout<<ele<<" "; 
      }
      std::cout<<std::endl; 
    }
  }
}



// Method 3:
// Gather all element-root information to the process of the root
  // Since the roots have all children, can directly use the information to reconstruct the sets; the only lacking part is the intersection. 
  // Step one: send to root
  // Step two: gather on root

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


void send_2_redistributed_processes(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();
  auto master = cp.master(); 
  auto& world = master->communicator(); 

  int nblocks = world.size();

  auto eles = b->eles; 
  for(auto& ele : eles) {
    std::string root = b->parent(ele); 
  
    int gid_root = hash_string(root, nblocks); 

    if(gid_root != gid) {
      std::pair<std::string, std::string> local_pair(ele, root); 

      auto& target = l->target(l->find(gid_root)); 

      Message send_msg; 
      send_msg.send_ele_parent_pair(local_pair); 

      cp.enqueue(target, send_msg); 

      // ==============

      Message send_msg_intersection; 
      send_msg_intersection.send_intersection(*(b->intersections.find(ele))); 

      cp.enqueue(target, send_msg_intersection); 

      b->erase_element(ele);  
    }    
  }
}

bool gather_on_redistributed_processes(Block* b, const diy::Master::ProxyWithLink& cp) {

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
        Message msg; 
        cp.dequeue(in[i], msg);

        if(msg.tag == "ele_parent_pair") {
            std::pair<std::string, std::string> pair; 
            msg.receive_ele_parent_pair(pair); 

            b->add(pair.first); 
            b->set_parent(pair.first, pair.second); 
        } else if(msg.tag == "intersection") {
          std::pair<std::string, intersection_t> pair; 
          msg.receive_intersection(pair);

          b->intersections.insert(pair); 
        } else {
          std::cout<<"Wrong! Tag is not correct: "<<msg.tag<<std::endl; 
        }

      }
    }
  }

  return true;
}

// Get sets of elements by redistributing data
void get_sets_redistributed(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {

  #ifdef FTK_HAVE_MPI
    double start = MPI_Wtime();
  #endif

  master.foreach(&send_2_redistributed_processes); 
  master.iexchange(&gather_on_redistributed_processes); 

  // #ifdef FTK_HAVE_MPI
  //   MPI_Barrier(world); 
  //   double end = MPI_Wtime();
  //   if(world.rank() == 0) {
  //     std::cout << "CCL: Gather Connected Components - Communication: " << end - start << " seconds. " << std::endl;
  //   }
  //   start = end; 
  // #endif

  Block* b = static_cast<Block*> (master.get(0)); // load block with local id 0

  std::map<std::string, std::set<std::string>> root2set; 

  // #ifdef FTK_HAVE_MPI
    // std::cout<<"# of elements on proc. " << world.rank() <<" : "<< b->eles.size()<<std::endl; 
  // #endif
  for(auto& ele : b->eles) {
    if(!b->is_root(b->parent(ele))) {
      std::cout<<"Wrong! The parent is not root! get_sets_redistributed()"<< std::endl; 
      std::cout<<ele<<" "<<b->parent(ele)<<" "<<b->parent(b->parent(ele))<<std::endl; 
    }

    root2set[b->parent(ele)].insert(ele);  
  }

  for(auto ite = root2set.begin(); ite != root2set.end(); ++ite) {
    results.push_back(ite->second);   
  }

  if(ISDEBUG) {
    for(int i = 0; i < results.size(); ++i) {
      std::cout<<"Set "<<i<<":"<<std::endl; 
      auto& ele_set = results[i]; 
      for(auto& ele : ele_set) {
        std::cout<<ele<<" "; 
      }
      std::cout<<std::endl; 
    }
  }

  // #ifdef FTK_HAVE_MPI
  //   MPI_Barrier(world); 
  //   end = MPI_Wtime();
  //   if(world.rank() == 0) {
  //     std::cout << "CCL: Gather Connected Components - Computation: " << end - start << " seconds. " << std::endl;
  //   }
  //   start = end; 
  // #endif
}

// Get sets of elements
void get_sets(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {
  // get_sets_on_p0(world, master, assigner, results); 
  // get_sets_on_roots(world, master, assigner, results); 
  get_sets_redistributed(world, master, assigner, results); 
}


#endif
