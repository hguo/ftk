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
#define OUTPUT_TIME_EACH_ROUND true
#define TRACK_PEAK_MEMORY false

namespace ftk {
template <class IdType=std::string>
struct distributed_union_find
{
  distributed_union_find() : eles(), id2parent() {

  }

  // Initialization
  // Add and initialize elements
  void add(IdType i) {
    eles.insert(i); 
    id2parent.insert(std::make_pair(i, i)); 

    #if TRACK_PEAK_MEMORY
      int64_t _memory_residule = i.length() * sizeof(char) * 3  + sizeof(std::string) * 3; 

      this->cur_memory += _memory_residule; 
      this->accumulative_memory += _memory_residule; 
    #endif
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

  // Queries

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
  }

  IdType parent(IdType i) {
    assert(has(i));

    return id2parent[i]; 
  }

  bool is_root(IdType i) {
    assert(has(i));

    return i == id2parent[i]; 
  }

public:
  std::set<IdType> eles; 

  // Statistics
  int nrounds = 0; // # of rounds
  int64_t peak_memory = 0;  // # peak used memory
  int64_t accumulative_memory = 0; // # accumulative used memory
  #if TRACK_PEAK_MEMORY
    int64_t cur_memory = 0;  // # current used memory
  #endif

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id2parent;
};

}

// DIY Block for distributed union-find
struct Block_Union_Find : public ftk::distributed_union_find<std::string> {
  Block_Union_Find(): nchanges(0), related_elements(), all_related_elements(), ele2gid(), distributed_union_find() { 
    #if TRACK_PEAK_MEMORY
      this->cur_memory = sizeof(*this) - sizeof(int) - sizeof(int64_t) * 3; // reduce the size of nrounds, peak_memory, accumulative_memory, and cur_memory
      this->accumulative_memory = this->cur_memory; 
      this->peak_memory = this->cur_memory; 
    #endif
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

      #if TRACK_PEAK_MEMORY
        int64_t  _memory_residule = ele.length() * sizeof(char) * 2 + sizeof(std::string) * 2 + sizeof(std::set<std::string>) + sizeof(bool); 
        this->cur_memory += _memory_residule; 
        this->accumulative_memory += _memory_residule; 
      #endif
    }
  }

  void erase_element(std::string ele) {
    this->nchanges += 1;

    #if TRACK_PEAK_MEMORY
      if(this->eles.find(ele) != this->eles.end()) {
        this->cur_memory -= ele.length() * sizeof(char) + sizeof(std::string); 
      }
    #endif

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

        #if TRACK_PEAK_MEMORY
          int64_t _memory_residule = related_ele.length() * sizeof(char) * 2 + sizeof(std::string) * 2; 
          this->cur_memory += _memory_residule; 
          this->accumulative_memory += _memory_residule; 
        #endif
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

    #if TRACK_PEAK_MEMORY
      for(auto& related_ele : this->related_elements[ele]) {
        this->cur_memory -= related_ele.length() * sizeof(char) + sizeof(std::string);   
      }
    #endif

    this->related_elements[ele].clear(); 
  }

  void erase_related_element(std::string ele, std::string related_element) {
    this->nchanges += 1;

    assert(this->has(ele)); 

    #if TRACK_PEAK_MEMORY
      if(this->related_elements[ele].find(related_element) != this->related_elements[ele].end()) {
        this->cur_memory -= related_element.length() * sizeof(char) + sizeof(std::string); 
      }
    #endif
      
    this->related_elements[ele].erase(related_element); 
  }

  bool has_gid(std::string ele) {
    if(this->ele2gid.find(ele) == this->ele2gid.end()) {
      return false; 
    }

    return true; 
  }

  void set_gid(std::string ele, int gid) {
    if(!this->has_gid(ele)) {
      this->nchanges += 1;

      if(gid >= 0) {
        this->ele2gid[ele] = gid; 

        #if TRACK_PEAK_MEMORY
          int64_t _memory_residule = ele.length() * sizeof(char) + sizeof(std::string) + sizeof(int); 
          this->cur_memory += _memory_residule; 
          this->accumulative_memory += _memory_residule; 
        #endif
      }
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


  void get_sets(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results);
  
  #if TRACK_PEAK_MEMORY
    void update_peak_memory() {
      #if TRACK_PEAK_MEMORY
        if(this->cur_memory > this->peak_memory) {
          this->peak_memory = this->cur_memory;   
        }
      #endif
    }
  #endif

public: 

  std::map<std::string, int> ele2gid; 

  int nchanges = 0; // # of processed unions per round = valid unions (united unions) + passing unions

  #if OUTPUT_TIME_EACH_ROUND
    std::stringstream time_ss; // store time for measure duration of each round;
    double time_start; 
  #endif

  std::map<std::string, bool> has_sent_gparent_query; 

private:
  // map element id to ids of its related elements
  // Can be optimized by ordered the related elements, put related elements on process first and ordered descendingly by ids
  
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

  void send_gparent(std::string& ele, std::string& gparent, int& gid_gparent) {
    tag = "gparent";

    strs.push_back(ele); 
    strs.push_back(gparent); 
    strs.push_back(std::to_string(gid_gparent)); 
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
  #if TRACK_PEAK_MEMORY
    b->update_peak_memory();
  #endif

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
      // If having a root

      auto& related_elements = b->get_related_elements(ele); 
      
      std::string found_related_ele = ele; // Find a the first related element has a smaller id than ele // cannot use auto&, since we will remove some elements from the original set; auto& will cause segmentation error
      for(auto& related_ele : related_elements) {
        if(b->has_gid(related_ele)) {
          if (related_ele < found_related_ele) {
            found_related_ele = related_ele; 
            break ;
          }
        }
      }

      // ================================================================

      if(found_related_ele < ele) {

        int rgid = b->get_gid(found_related_ele); 
        // std::cout<<gid<<" "<<rgid<<std::endl; 
        
        b->set_parent(ele, found_related_ele); 

        #if ISDEBUG
          std::cout<<ele<<" -1> "<<found_related_ele<<std::endl; 
        #endif

        b->erase_related_element(ele, found_related_ele); 
      } 

    }
  }

}


void compress_path(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  for(auto& ele : b->eles) {

    if(!b->is_root(ele)) {
      std::string parent = b->parent(ele);

      if(b->has(parent)) {
        
        if(!b->is_root(parent)) {
          std::string grandparent = b->parent(parent);
          b->set_parent(ele, grandparent); // Freely
        }

      } else {
        
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


// Distributed path compression - Step Two
void distributed_save_gparent(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  std::string ele; 
  std::string grandpar; 
  int gid_grandparent; 

  msg.rec_gparent(ele, grandpar, gid_grandparent); 
  b->set_gid(grandpar, gid_grandparent); 
  
  std::string parent = b->parent(ele);
  if(grandpar < parent) {
    b->set_parent(ele, grandpar);   
  }

  #if ISDEBUG
    std::cout<<ele<<" -2> "<<grandpar<<std::endl; 
  #endif

  b->has_sent_gparent_query[ele] = false;
}

// Distributed path compression - Step Three
void distributed_answer_gparent_query(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp, Message_Union_Find msg) {
  diy::Link* l = cp.link();

  std::string child; 
  std::string parent; 
  int gid_child; 

  msg.rec_gparent_query(child, parent, gid_child); 

  // ******* commented 08.29.2020
  // b->set_gid(child, gid_child);

  std::string grandparent = b->parent(parent); 
  int gid_grandparent = b->get_gid(grandparent);

  // Set child's parent to grandparent
  Message_Union_Find send_msg; 

  send_msg.send_gparent(child, grandparent, gid_grandparent); 
  cp.enqueue(l->target(l->find(gid_child)), send_msg); 
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


void total_changes(Block_Union_Find* b, const diy::Master::ProxyWithLink& cp) {
  cp.collectives()->clear();

  cp.all_reduce(b->nchanges, std::plus<int>()); 
  b->nchanges = 0;

  b->nrounds += 1;
  
  #if TRACK_PEAK_MEMORY
    b->update_peak_memory(); 
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
      blocks[0]->time_ss << gids[0]; 

      MPI_Barrier(world); 
      blocks[0]->time_start = MPI_Wtime();
    #endif
  #endif
  

  // Only supports exchange
  bool all_done = false;

  master.foreach(&query_gid); // will also update peak memory
  master.exchange();

  while(!all_done) {
    #if OUTPUT_TIME_EACH_ROUND
      blocks[0]->time_ss << "," << MPI_Wtime() - blocks[0]->time_start; 
    #endif

    exchange_process(master); 

    #if OUTPUT_TIME_EACH_ROUND
      blocks[0]->time_ss << "," << MPI_Wtime() - blocks[0]->time_start; 
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
      if(!filename_time_uf_w.empty()) {

        MPI_Status status;
        MPI_File fh;

        blocks[0]->time_ss<<std::endl;
        const std::string buf = blocks[0]->time_ss.str();

        MPI_File_open(world, filename_time_uf_w.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        
        MPI_File_write_ordered(fh, buf.c_str(), buf.length(), MPI_CHAR, &status);

        MPI_File_close(&fh);
      }
    #endif
  #endif
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
