#ifndef _FTK_DISTRIBUTED_UNION_FIND_EXCHANGE_H
#define _FTK_DISTRIBUTED_UNION_FIND_EXCHANGE_H

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <iostream>

#include <ftk/external/diy/mpi.hpp>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/assigner.hpp>
#include <ftk/external/diy/serialization.hpp>

// Union-find algorithm with distributed-memory parallelism
  // All elements need to be added by invoking the add function at the initialization stage. 

// Reference 
  // Distribute Union-Find in Paper: "Evaluation of connected-component labeling algorithms for distributed-memory systems"

// Add the sparse representation by using Hash map/table 

#define IEXCHANGE 0
#define ISDEBUG   0

namespace ftk {
template <class IdType=std::string>
struct distributed_union_find_exchange
{
  distributed_union_find_exchange() {
    
  }

  // Initialization
  // Add and initialize elements
  void add(IdType i) {
    eles.insert(i); 
    id2parent[i] = i;
  }

  // Operations

  void set_parent(IdType i, IdType par) {
    if(!has(i)) {
      std::cout<< "No element parent. " <<std::endl;
      exit(0); 
    }

    id2parent[i] = par; 
  }


  // Queries

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
  }

  IdType parent(IdType i) {
    if(!has(i)) {
      std::cout<< "No element parent. " <<std::endl;
      exit(0); 
    }

    return id2parent[i]; 
  }

  bool is_root(IdType i) {
    if(!has(i)) {
      std::cout<< "No element. " <<std::endl;
      exit(0); 
    }

    return i == id2parent[i]; 
  }

public:
  std::set<IdType> eles; 

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id2parent;
};

}

// ==========================================================

struct Block : public ftk::distributed_union_find_exchange<std::string> {
  Block(): nchanges(0), distributed_union_find_exchange() { 
    
  }

  // add an element
  void add(std::string ele) {
    distributed_union_find_exchange::add(ele); 
    this->ele2gid[ele] = -1;
  }

  void add_related_element(std::string ele, std::string related_ele) {
      this->related_elements[ele].push_back(related_ele); 
      this->ele2gid[related_ele] = -1; 
  }

public: 
  // map element id to ids of its related elements
  // Can be optimized by ordered the related elements, put related elements on process first and ordered decreingly by ids
  std::map<std::string, std::vector<std::string>> related_elements; 
  std::map<std::string, int> ele2gid; 

  int nchanges; // # of processed unions per round = valid unions (united unions) + passing unions
};

// ==========================================================

struct Message {

// Send message for query
  void send_gid_query(std::string* ele) {
    tag = "gid_query"; 
    strs.push_back(*ele);  
  }

  void send_gid_response(std::string* ele, int* gid) {
    tag = "gid_response"; 

    strs.push_back(*ele); 
    strs.push_back(std::to_string(*gid)); 
  }


  // grandparent query
  void send_gparent_query(std::string* ele, std::string* parent) {
    tag = "gparent_query";

    strs.push_back(*ele); 
    strs.push_back(*parent); 
  }

  void send_gparent_response(std::string* ele, std::string* gparent, int* gid_gparent) {
    tag = "gparent_response";

    strs.push_back(*ele); 
    strs.push_back(*gparent); 
    strs.push_back(std::to_string(*gid_gparent)); 
  }

  void send_union(std::string* ele, std::string* parent, std::string* related_ele, int* gid_related_ele) {
    tag = "union"; 

    strs.push_back(*ele); 
    strs.push_back(*parent); 
    strs.push_back(*related_ele); 
    strs.push_back(std::to_string(*gid_related_ele)); 
  }

  void send_endpoint(std::string* ele, std::string* parent, std::string* related_ele, int* gid_parent) {
    tag = "endpoint"; 

    strs.push_back(*ele); 
    strs.push_back(*parent); 
    strs.push_back(*related_ele); 
    strs.push_back(std::to_string(*gid_parent)); 
  }

  void send_ele_parent_pairs(std::vector<std::pair<std::string, std::string>>& pairs) {
    tag = "ele_parent_pairs"; 

    for(auto& pair : pairs) {
      strs.push_back(pair.first); 
      strs.push_back(pair.second); 
    }
  }


// Receive message

  void rec_gid_query(std::string** ele) {
    *ele = &strs[0]; 
  }

  void rec_gid_response(std::string** ele, std::string** gid_str) {
    *ele = &strs[0]; 
    *gid_str = &strs[1]; 
  }

  void rec_gparent_query(std::string** ele, std::string** parent) {
    *ele = &strs[0]; 
    *parent = &strs[1]; 
  }

  void rec_gparent_response(std::string** ele, std::string** gparent, std::string** gparent_gid) {
    *ele = &strs[0]; 
    *gparent = &strs[1];
    *gparent_gid = &strs[2]; 
  }

  void rec_union(std::string** ele, std::string** parent, std::string** related_ele, std::string** rgid) {
    *ele = &strs[0]; 
    *parent = &strs[1];
    *related_ele = &strs[2]; 
    *rgid = &strs[3]; 
  }

  void rec_endpoint(std::string** ele, std::string** parent, std::string** related_ele, std::string** parent_gid) {
    *ele = &strs[0]; 
    *parent = &strs[1];
    *related_ele = &strs[2]; 
    *parent_gid = &strs[3]; 
  }


  void receive_ele_parent_pairs(std::vector<std::pair<std::string, std::string>>& pairs) {
    for(int i = 0; i < strs.size(); i += 2) {
      pairs.push_back(std::make_pair(strs[i], strs[i+1]));
    }
  }


public:
  // 0 - string
  // 1 - vectors of string
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
    std::string ele = it->first; 

    if(b->has(ele)) {
      it->second = gid; 
    } else {
      for (int i = 0; i < l->size(); ++i) {
        Message msg = Message(); 
        msg.send_gid_query(&ele); 

        cp.enqueue(l->target(i), msg);
      }
    }
  }

  // std::cout<<std::endl; 
}

void answer_gid(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();
  // diy::Master* master = cp.master(); 

  // std::cout<<"Answer: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 

  while(!cp.empty_incoming_queues()) {
    // Save unions from other blocks
    std::vector<int> in; // gids of incoming neighbors in the link
    cp.incoming(in);

    // for all neighbor blocks
    // dequeue data received from this neighbor block in the last exchange
    for (unsigned i = 0; i < in.size(); ++i) {
      int from_gid = in[i]; 

      if(cp.incoming(from_gid)) {

        Message msg; 
        cp.dequeue(from_gid, msg);

        std::string* ele_ptr; 
        msg.rec_gid_query(&ele_ptr); 

        if(b->has(*ele_ptr)) {
          Message send_msg; 
          send_msg.send_gid_response(ele_ptr, &gid); 

          // std::cout<<ele<<"-";
          // std::cout<<gid;
          // std::cout<<std::endl; 

          cp.enqueue(l->target(l->find(from_gid)), send_msg);
          // cp.enqueue(l->target(l->find(in[i])), 1);
        }
      }
    }
  }

  // std::cout<<std::endl; 
}

void save_gid(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 

  // std::cout<<"Save: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 

  while(!cp.empty_incoming_queues()) {
    // Save unions from other blocks
    std::vector<int> in; // gids of incoming neighbors in the link
    cp.incoming(in);

    // for all neighbor blocks
    // dequeue data received from this neighbor block in the last exchange
    for (unsigned i = 0; i < in.size(); ++i) {
      int from_gid = in[i]; 

      if(cp.incoming(from_gid)) {
        Message msg; 
        cp.dequeue(from_gid, msg);

        std::string* ele_ptr; 
        std::string* gid_str; 
        msg.rec_gid_response(&ele_ptr, &gid_str); 

        int gid = std::stoi(*gid_str); 

        // std::cout<<_pair.first<<" - "<<_pair.second<<std::endl; 

        b->ele2gid[*ele_ptr] = gid; 

        // std::cout<<"i-"<<i<<'-'<<in[i]<<std::endl;
        // int t; 
        // cp.dequeue(in[i], t); 
        // std::cout<<t<<std::endl; 
      }
    }
  }

  // std::cout<<std::endl; 
}


void unite_once(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 

  for(auto it = b->eles.begin(); it != b->eles.end(); ++it) {
    std::string ele = *it; 

    // if(ele == "7") {
    //   std::cout<<ele<<": "<<b->related_elements[ele].size()<<std::endl; 
    //   for(std::vector<std::string>::iterator it_vec = b->related_elements[ele].begin(); it_vec != b->related_elements[ele].end(); ++it_vec) {
    //     std::string related_ele = *it_vec; 

    //     std::cout<<ele<<"~~~~"<<related_ele<<std::endl; 
    //   }
    // }

    if(b->is_root(ele)) {
      for(std::vector<std::string>::iterator it_vec = b->related_elements[ele].begin(); it_vec != b->related_elements[ele].end(); ++it_vec) {
        std::string related_ele = *it_vec; 

        int rgid = b->ele2gid[related_ele]; 
        // std::cout<<gid<<" "<<rgid<<std::endl; 
        if(rgid == -1) {
          continue; 
        }

        // Unite with an element with larger id or smaller id is better? 
          // Here is a larger id

        // if(related_ele > ele) {
        // if(std::stoi(related_ele) > std::stoi(ele)) {

        if(rgid < gid || (rgid == gid && related_ele < ele)) {
        // if(rgid > gid || (rgid == gid && std::stoi(related_ele) > std::stoi(ele))) {
          b->set_parent(ele, related_ele); 
          b->nchanges += 1;
          if(ISDEBUG) {
            std::cout<<ele<<" -1> "<<related_ele<<std::endl; 
          }

          b->related_elements[ele].erase(it_vec); 

          break ; 
        }
      }
    }
  }
}


// Distributed path compression - Step One
void distributed_query_gparent(Block* b, const diy::Master::ProxyWithLink& cp) {
  // Across processors

  // Why querying grandparent instead of root? 
    // We can only guarantee the processor of its parent stores its grandparent. Getting root may require additonal communications and lengthen the communication chain. 

  // If the parent and the element are not in the same processor, the element will keep asking its grandparent. 
    // If its parent is already the root, there are no answers are given back; this block will be cooled down. 
    // But a dead lock case: 
      // An extreme case is that if processor A asks the grandparent to processor B, and B also asks A,  
      // then A will ask again, and as well as B. It becomes an endless loop. 

      // Solutions: 
      // Solution 1: (Default)
        // When uniting element to an element with larger ID, consider the rank of processor as well. 
        // Elements of a processor can only be united with processors with higher ranks. 
      // Solution 2: 
        // Assign a rank to each elemet, elements in low rank processor have low rank; elements in high rank processor have high rank. 
        // Ensure no loops
      


  diy::Link* l = cp.link();
  for(auto it = b->eles.begin(); it != b->eles.end(); ++it) {
    std::string ele = *it; 
    std::string parent = b->parent(ele); 

    if(!b->is_root(ele)) {
      if(!b->has(parent)) {
        int gid = b->ele2gid[parent]; 
        // if(gid == -1) {
        //   continue ;
        // }

        Message send_msg; 
        send_msg.send_gparent_query(&ele, &parent); 

        cp.enqueue(l->target(l->find(gid)), send_msg);

        // std::cout<<"Query Parent: "<<ele<<" - "<<parent<<" - "<<gid<<std::endl; 

        // if(ele == "11" || ele == "7") {
        //   std::cout<<ele<<" ~ "<<parent<<std::endl; 
        // }

      }
    }
  }
}

// Distributed path compression - Step Two
void distributed_answer_gparent(Block* b, const diy::Master::ProxyWithLink& cp, Message msg, int from_gid) {
  diy::Link* l = cp.link();

  std::string* ele_ptr; 
  std::string* parent_ptr; 
  
  msg.rec_gparent_query(&ele_ptr, &parent_ptr); 

  if(b->has(*parent_ptr) && !b->is_root(*parent_ptr)) {
    std::string grandparent = b->parent(*parent_ptr); 

    Message send_msg; 
    send_msg.send_gparent_response(ele_ptr, &grandparent, &b->ele2gid[grandparent]); 

    // std::cout<<*ele_ptr<<" - "<<*parent_ptr<<" - "<<grandparent<<" - "<<b->ele2gid[grandparent]<<std::endl; 

    cp.enqueue(l->target(l->find(from_gid)), send_msg); 


    // if(*ele_ptr == "4") {
    //   std::cout<<*ele_ptr<<" ~ "<<*parent_ptr<<" ~ "<<grandparent<<std::endl; 
    // }
  }
}

// Distributed path compression - Step Three
void distributed_save_gparent(Block* b, const diy::Master::ProxyWithLink& cp, Message msg) {
  std::string* ele_ptr; 
  std::string* grandpar_ptr; 
  std::string* grandpar_gid_str; 

  msg.rec_gparent_response(&ele_ptr, &grandpar_ptr, &grandpar_gid_str); 

  int gparent_gid = std::stoi(*grandpar_gid_str); 
  
  // std::cout<<*ele_ptr << " - " << b->parent(*ele_ptr) <<" - "<< *grandpar_ptr<<" - "<<*grandpar_gid_ptr<<std::endl; 

  // if(*ele_ptr == "4") {
  //   std::cout<<*ele_ptr<<" ~ "<<*grandpar_ptr<<std::endl; 
  // }


  b->set_parent(*ele_ptr, *grandpar_ptr); 
  b->nchanges += 1;

  if(ISDEBUG) {
    std::cout<<*ele_ptr<<" -2> "<<*grandpar_ptr<<std::endl; 
  }

  b->ele2gid[*grandpar_ptr] = gparent_gid; 
}


void distributed_pass_unions(Block* b, const diy::Master::ProxyWithLink& cp) {
  diy::Link* l = cp.link();

  for(auto ite_ele = b->eles.begin(); ite_ele != b->eles.end(); ++ite_ele) {
    std::string ele = *ite_ele;

    if(!b->is_root(ele)) {
      auto src = &b->related_elements[ele]; 

      if(src->size() > 0) {
        std::string par = b->parent(ele); 

        if(!b->has(par)) {
          // Communicate with other processes
          int gid = b->ele2gid[par]; 

          for(auto ite_related_ele = src->begin(); ite_related_ele != src->end(); ++ite_related_ele) {
            Message send_msg; 
            std::string related_ele = *ite_related_ele; 
            send_msg.send_union(&ele, &par, &related_ele, &b->ele2gid[related_ele]); 

            cp.enqueue(l->target(l->find(gid)), send_msg); 
          }

          src->clear(); 
        }
      }
    }
  }
}


void local_compress_path(Block* b, const diy::Master::ProxyWithLink& cp) {
  // std::cout<<"Compress Path: "<<std::endl; 
  // std::cout<<"Block ID: "<<cp.gid()<<std::endl; 

  // Local computation
  // Locally compress path
  for(auto it = b->eles.begin(); it != b->eles.end(); ++it) {
    std::string ele = *it; 
    std::string parent = b->parent(ele); 

    if(!b->is_root(ele)) {
      if(b->has(parent)) {
        // if(ele == "4") {
        //   std::cout<<ele<<" ~ "<<parent<<" ~ "<<b->has(parent)<<std::endl; 
        // }

        if(!b->is_root(parent)) {
          std::string grandparent = b->parent(parent); 

          b->set_parent(ele, grandparent); 
          b->nchanges += 1;

          if(ISDEBUG) {
            std::cout<<ele<<" -2> "<<grandparent<<std::endl; 
          }
        }
      }
    }
  }
}


// Tell related elements that the endpoints have been changed
void local_update_endpoint(Block* b, const diy::Master::ProxyWithLink& cp, std::string* ele_ptr, std::string* par_ptr, int* pgid_ptr, std::string* related_ele_ptr) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  for(auto ite = b->related_elements[*related_ele_ptr].begin(); ite != b->related_elements[*related_ele_ptr].end(); ++ite) {

    if(*ite == *ele_ptr) {
      *ite = *par_ptr; 
      b->nchanges += 1;
      
      return ;
    }
  }

  // If we cannot find the union of the related element, there are two possibilities
    // (1) The union has been processed. Then, we don't need to do anything. 
    // (2) The union has been passed to the parent of this related element, also we have two situations: 
      // (1) The union has been successfully passed to its parent, we can notify its parent that this endpoint has changed. 
      // (2) If the union is also a message and has not been received by its parent, we store the edge. 

  // if(*ele_ptr == "9") {
  //   std::cout<<*related_ele_ptr<<" - "<<b->has(*related_ele_ptr) <<" - "<< b->parent(*related_ele_ptr) <<std::endl; 
  // }

  if(!b->is_root(*related_ele_ptr)) {
    std::string parent_related_ele = b->parent(*related_ele_ptr); // the parent of the related element
    int r_p_gid = b->ele2gid[parent_related_ele]; 

    if(r_p_gid == gid) {
      local_update_endpoint(b, cp, ele_ptr, par_ptr, pgid_ptr, &parent_related_ele); 
    } else {
      if(r_p_gid != -1) {
        Message send_msg; 
        send_msg.send_endpoint(ele_ptr, par_ptr, &parent_related_ele, pgid_ptr); 

        // std::cout<<"!!!"<<*ele_ptr<<" - "<<*par_ptr<<" - "<<parent_related_ele<<" - "<<*pgid_ptr<<std::endl; 

        cp.enqueue(l->target(l->find(r_p_gid)), send_msg); 
        // This message should get completed, so we set a change to avoid the algorithm ends. 
        b->nchanges += 1;
      }
    }
  } else { // To remove possibe side-effects, we store the edge. 
    b->related_elements[*related_ele_ptr].push_back(*par_ptr); 
  }
}

void local_pass_unions(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  // Local computation
  // Pass unions of elements in this block to their parents, save these unions
  for(auto ite_ele = b->eles.begin(); ite_ele != b->eles.end(); ++ite_ele) {
    std::string ele = *ite_ele; 

    if(!b->is_root(ele)) {
      auto src = &b->related_elements[ele]; 
      if(src->size() > 0) {
        std::string par = b->parent(ele); 

        if(b->has(par)) {
          // Update directly, if the realted element is in the block

          int p_gid = gid; 

          auto dest = &b->related_elements[par]; 
          dest->insert(
            dest->end(),
            src->begin(),
            src->end()
          );
          b->nchanges += src->size();

          // tell related elements, the end point has changed to its parent
            // Since has changed locally, the message will one send once. 

          for(auto ite_related_ele = src->begin(); ite_related_ele != src->end(); ++ite_related_ele) {
            std::string related_ele = *ite_related_ele; 
            int r_gid = b->ele2gid[related_ele]; 

            if(r_gid == gid) {
              local_update_endpoint(b, cp, &ele, &par, &p_gid, &related_ele); 
            } else {
              if(r_gid != -1) {
                Message send_msg; 
                send_msg.send_endpoint(&ele, &par, &related_ele, &p_gid); 

                cp.enqueue(l->target(l->find(r_gid)), send_msg); 
              }
            }
          }

          src->clear(); 
        }
      }
    }
  }
}

// update unions of related elements
  // Tell related elements that the endpoints have been changed
void distributed_save_union(Block* b, const diy::Master::ProxyWithLink& cp, Message msg) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  // std::cout<<"Save Unions: "<<std::endl; 
  // std::cout<<"Block ID: "<<cp.gid()<<std::endl; 

  std::string* ele_ptr; 
  std::string* par_ptr; 
  std::string* related_ele_ptr; 
  std::string* rgid_str; 
  
  msg.rec_union(&ele_ptr, &par_ptr, &related_ele_ptr, &rgid_str); 
  int rgid = std::stoi(*rgid_str); 

  // std::cout<<*ele_ptr<<" - "<<*par_ptr << " - " << *related_ele_ptr <<" - "<< *gid_ptr<<std::endl; 

  b->related_elements[*par_ptr].push_back(*related_ele_ptr); 
  b->ele2gid[*related_ele_ptr] = rgid; 
  b->nchanges += 1;

  // tell related elements, the end point has changed to its parent

  if(rgid == gid) {
    // Tell related elements that the endpoints have been changed
    local_update_endpoint(b, cp, ele_ptr, par_ptr, &gid, related_ele_ptr); 
  } else {
    if(rgid != -1) {
      Message send_msg; 
      send_msg.send_endpoint(ele_ptr, par_ptr, related_ele_ptr, &gid); 

      cp.enqueue(l->target(l->find(rgid)), send_msg);   
    }
  }

  // std::cout<<std::endl; 
}

void distributed_update_endpoint(Block* b, const diy::Master::ProxyWithLink& cp, Message msg) {
  std::string* ele_ptr; 
  std::string* par_ptr; 
  std::string* related_ele_ptr; 
  std::string* pgid_str; 

  msg.rec_endpoint(&ele_ptr, &par_ptr, &related_ele_ptr, &pgid_str); 
  int pgid = std::stoi(*pgid_str); 

    // std::cout<<*ele_ptr<<" - "<<*par_ptr << " - " << *related_ele_ptr<<std::endl; 

  // Tell related elements that the endpoints have been changed
  b->ele2gid[*par_ptr] = pgid; 
  local_update_endpoint(b, cp, ele_ptr, par_ptr, &pgid, related_ele_ptr); 
}

void distributed_computation(Block* b, const diy::Master::ProxyWithLink& cp) {
  // This message is only needed to send once, and get the feedback
  if(cp.empty_outgoing_queues()) { // reduce sending multiple message, which cannot avoid all
    distributed_query_gparent(b, cp); // Start of distributed pass compression
  }

  // distributed_query_gparent(b, cp);
  distributed_pass_unions(b, cp); 
}

void local_computation(Block* b, const diy::Master::ProxyWithLink& cp) {
  unite_once(b, cp); 
  local_compress_path(b, cp); 
  local_pass_unions(b, cp); 
}

void received_msg(Block* b, const diy::Master::ProxyWithLink& cp, Message msg, int from_gid) {
  // std::cout<<msg.tag<<std::endl; 

  // if(msg.tag == "gid_query") {
  //   answer_gid(b, cp, msg, from_gid); 
  // } else if (msg.tag == "gid_response") {
  //   save_gid(b, cp, msg); 
  // } 

  if(msg.tag == "gparent_query") {
    distributed_answer_gparent(b, cp, msg, from_gid); 
  } else if(msg.tag == "gparent_response") {
    distributed_save_gparent(b, cp, msg); 
  } else if(msg.tag == "union") {
    distributed_save_union(b, cp, msg); 
  } else if(msg.tag == "endpoint") {
    distributed_update_endpoint(b, cp, msg); 
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
  local_computation(b, cp); 
  receive_msg(b, cp); 

  distributed_computation(b, cp); 

  if(ISDEBUG) {
    std::cout<<"================================="<<std::endl; 
  }

  return true; 
}

void iexchange_process(diy::Master& master) {
  master.iexchange(&union_find_iexchange); 
}

void union_find_exchange(Block* b, const diy::Master::ProxyWithLink& cp) {
  local_compress_path(b, cp); 
  local_pass_unions(b, cp); 
  receive_msg(b, cp); 

  distributed_pass_unions(b, cp); 

  total_changes(b, cp); 
}


// Method 0
// void exchange_process(diy::Master& master) {
//   master.foreach(&unite_once);

//   master.foreach(&local_compress_path);

//   master.foreach(&distributed_query_gparent);
  
//   master.exchange();                 
//   // master.foreach(&distributed_answer_gparent_exchange);
//   // std::cout<<"!!!!!distributed_answer_gparent!!!!!"<<std::endl; 
//   master.foreach(&receive_msg); // distributed_answer_gparent
  
//   master.exchange();
//   // master.foreach(&distributed_save_gparent_exchange);
//   // std::cout<<"!!!!!distributed_save_gparent!!!!!"<<std::endl; 
//   master.foreach(&receive_msg); // distributed_save_gparent
  
//   master.foreach(&distributed_pass_unions);
  
//   master.exchange();
//   // master.foreach(&distributed_save_union_exchange);
//   // std::cout<<"!!!!!distributed_save_union!!!!!"<<std::endl; 
//   master.foreach(&receive_msg); // distributed_save_union

//   master.foreach(&local_pass_unions);

//   master.exchange();
//   // master.foreach(&distributed_update_endpoint_exchange);
//   // std::cout<<"!!!!!distributed_update_endpoint!!!!!"<<std::endl; 
//   master.foreach(&receive_msg); // distributed_update_endpoint
//   // master.exchange();
//   // std::cout<<"!!!!!!!!!!"<<std::endl; 
//   // master.foreach(&receive_msg);
//   // master.exchange();
//   // std::cout<<"!!!!!!!!!!"<<std::endl; 
//   // master.foreach(&receive_msg);
//   // master.exchange();

//   master.foreach(&total_changes); // can use a reduce all to check whether every thing is done. 
  
//   master.exchange();
//   // It is possible that some endpoints are still tranferring
//     // While two consecutive exchange() without receiving will clear the buffer
//     // Which means, in-between each pair of exchange should contain at least one receive
//     // Hence, an additional receive is needed
//   master.foreach(&receive_msg);

//   std::cout<<"================================="<<std::endl; 
// }

// Method 1
// gparent query and gparent answer should be synchronized. 
  // otherwise it is possible that one query message is transferring, but the program ends. 
  // Other parts are not required to synchronize. 
void exchange_process(diy::Master& master) {
  master.foreach(&unite_once);
  master.foreach(&local_compress_path);
  master.foreach(&distributed_query_gparent);
  
  master.exchange();                 
  master.foreach(&receive_msg); // distributed_answer_gparent + additional responses

  master.exchange();
  master.foreach(&receive_msg); // distributed_save_gparent + additional responses

  // master.foreach(&union_find_exchange);
  // master.exchange();

  master.foreach(&local_compress_path);
  master.foreach(&local_pass_unions);
  master.foreach(&distributed_pass_unions);
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
  
  master.foreach(&query_gid);
  // master.execute(); 
  master.exchange();
  master.foreach(&answer_gid);  
  master.exchange();
  master.foreach(&save_gid);

  for(int i = 0; i < gids.size(); ++i) {
    Block* b = static_cast<Block*> (master.get(i));
    for(auto& ele : b->ele2gid) {
      if(ele.second < 0) {
        std::cout<<world.rank()<<": gid is negative: "<<ele.first<<" "<<ele.second<<" - Since some elements are missing / not added into blocks. The invalid elements also need to be added. "<<std::endl; 
      }
    }
  }
  
  if(IEXCHANGE) { // for iexchange
    bool all_done = false;
    while(!all_done) {
      iexchange_process(master);  

      // For i-exchange, it is possible that, finally, some elements do not point to their roots
        // Hence, we need to repeat until all elements point to their roots. 

      master.foreach(&total_changes);
      master.exchange();

      int total_changes = master.proxy(master.loaded_block()).read<int>();
      all_done = total_changes == 0;
    }
  } else { // for exchange
    bool all_done = false;
    while(!all_done) {
      exchange_process(master); 

      // int total_changes;
      // for (int i = 0; i < master.size(); i++)
      //     total_changes = master.proxy(i).get<size_t>();

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

// Gather all element-parent information to the root block
bool gather_2_root(Block* b, const diy::Master::ProxyWithLink& cp) {
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

          if(msg.tag != "ele_parent_pairs") {
            std::cout<<"Wrong! Tag is not correct: "<<msg.tag<<std::endl; 
          }

          std::vector<std::pair<std::string, std::string>> pairs; 
          msg.receive_ele_parent_pairs(pairs); 

          // std::cout<<"Received Pairs: "<<pairs.size()<<std::endl; 

          for(auto& pair : pairs) {
            b->add(pair.first); 
            b->set_parent(pair.first, pair.second); 
          }
        }
      }
    }
  } else {
    std::vector<std::pair<std::string, std::string>> local_pairs; 

    for(auto& ele : b->eles) {
      std::string parent = b->parent(ele); 

      local_pairs.push_back(std::make_pair(ele, parent)); 
    }

    Message send_msg; 
    send_msg.send_ele_parent_pairs(local_pairs); 

    // std::cout<<"Sent Pairs: "<<local_pairs.size()<<std::endl; 

    cp.enqueue(l->target(l->find(0)), send_msg); 
  }

  return true;
}

// Get sets of elements
void get_sets(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {
  master.iexchange(&gather_2_root); 

  std::vector<int> gids;                     // global ids of local blocks
  assigner.local_gids(world.rank(), gids);

  if(world.rank() == 0) {
    Block* b = static_cast<Block*> (master.get(0)); 

    std::map<std::string, std::set<std::string>> root2set; 

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

#endif
