// Serilization
// https://github.com/diatomic/diy/blob/master/examples/simple/block.h


// Probibly, can be optimized by "An O(logn) parallel connectivity algorithm"

#define IEXCHANGE 1


#include <vector>
#include <iostream>

#include <ftk/external/diy/mpi.hpp>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/assigner.hpp>
#include <ftk/external/diy/serialization.hpp>

#include <ftk/basic/distributed_union_find.hh>

typedef std::pair<std::string, int> ele_gid;  // stores element and its global block id
typedef std::map<std::string, std::vector<std::string>> r_ele_map; 

typedef std::map<std::string, int> ele2gid_map; 

struct Block : public ftk::distributed_union_find<std::string> {
  Block(): nchanges(0), distributed_union_find() { 
    
  }

public: 
  int nchanges; // # of processed unions per round = valid unions (united unions) + passing unions

  // map element id to ids of its related elements
  // Can be optimized by ordered the related elements, put related elements on process first and ordered decreingly by ids
  r_ele_map related_elements; 
  ele2gid_map ele2gid; 
};

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


public:
  // 0 - string
  // 1 - vectors of string
  std::string tag; 

  std::vector<std::string> strs; 

}; 

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

// void*   create_block()                      { return new Block; }
// void    destroy_block(void* b)              { delete static_cast<Block*>(b); }
// void    save_block(const void* b,
//                    diy::BinaryBuffer& bb)   { diy::save(bb, *static_cast<const Block*>(b)); }
// void    load_block(void* b,
//                    diy::BinaryBuffer& bb)   { diy::load(bb, *static_cast<Block*>(b)); }


// std::string getID(int gid, int j) {
//   return std::to_string(gid) + '_' + std::to_string(j); 
// }


void query_gid(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  // std::cout<<"Query: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 
  
  // Can be optimized more, link by a queue. 
  for(ele2gid_map::iterator it = b->ele2gid.begin(); it != b->ele2gid.end(); ++it) {
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
        // if(rgid == -1) {
        //   continue; 
        // }

        // Unite with an element with larger id or smaller id is better? 
          // Here is a larger id

        // if(related_ele > ele) {
        // if(std::stoi(related_ele) > std::stoi(ele)) {

        if(rgid > gid || (rgid == gid && related_ele > ele)) {
        // if(rgid > gid || (rgid == gid && std::stoi(related_ele) > std::stoi(ele))) {
          b->set_parent(ele, related_ele); 
          b->nchanges += 1;
          std::cout<<ele<<" -1> "<<related_ele<<std::endl; 

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
  std::cout<<*ele_ptr<<" -2> "<<*grandpar_ptr<<std::endl; 

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
        if(ele == "4") {
          std::cout<<ele<<" ~ "<<parent<<" ~ "<<b->has(parent)<<std::endl; 
        }

        if(!b->is_root(parent)) {
          std::string grandparent = b->parent(parent); 

          b->set_parent(ele, grandparent); 
          b->nchanges += 1;
          std::cout<<ele<<" -2> "<<grandparent<<std::endl; 
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
    std::string parent_related_ele = b->parent(*related_ele_ptr); 
    int r_p_gid = b->ele2gid[parent_related_ele]; 

    if(r_p_gid == gid) {
      local_update_endpoint(b, cp, ele_ptr, par_ptr, pgid_ptr, &parent_related_ele); 
    } else {
      Message send_msg; 
      send_msg.send_endpoint(ele_ptr, par_ptr, &parent_related_ele, pgid_ptr); 

      // std::cout<<"!!!"<<*ele_ptr<<" - "<<*par_ptr<<" - "<<parent_related_ele<<" - "<<*pgid_ptr<<std::endl; 

      cp.enqueue(l->target(l->find(r_p_gid)), send_msg); 
      // This message should get completed, so we set a change to avoid the algorithm ends. 
      b->nchanges += 1;
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
              Message send_msg; 
              send_msg.send_endpoint(&ele, &par, &related_ele, &p_gid); 

              cp.enqueue(l->target(l->find(r_gid)), send_msg); 
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
    Message send_msg; 
    send_msg.send_endpoint(ele_ptr, par_ptr, related_ele_ptr, &gid); 

    cp.enqueue(l->target(l->find(rgid)), send_msg); 
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

  std::cout<<"================================="<<std::endl; 

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

  std::cout<<"================================="<<std::endl; 
}



int main(int argc, char* argv[]) {
  diy::mpi::environment     env(argc, argv);
  diy::mpi::communicator    world;
  
  int                       nthreads = 1; 
  
  int                       nblocks = 4*world.size();
  int                       memblocks = -1; //2; // number of blocks to store in memory


  diy::FileStorage          storage("./DIY.XXXXXX");

  // diy::Master               master(world,
  //                                  nthreads,
  //                                  memblocks,
  //                                  &create_block,
  //                                  &destroy_block,
  //                                  &storage,
  //                                  &save_block,
  //                                  &load_block);

  diy::Master               master(world, nthreads);

  srand(time(NULL));

  diy::ContiguousAssigner   assigner(world.size(), nblocks);
  // diy::RoundRobinAssigner     assigner(world.size(), nblocks);

  std::vector<int> gids;                     // global ids of local blocks
  assigner.local_gids(world.rank(), gids);   // get the gids of local blocks for a given process rank 

  // for the local blocks in this processor
  for (unsigned i = 0; i < gids.size(); ++i) {
    int gid = gids[i];

    diy::Link*    link = new diy::Link; 
    for(unsigned j = 0; j < gids.size(); ++j) {
      int neighbor_gid = gids[j]; 
      if(gid != neighbor_gid) {
        diy::BlockID  neighbor;
        neighbor.gid  = neighbor_gid;  // gid of the neighbor block
        neighbor.proc = assigner.rank(neighbor.gid);  // process of the neighbor block
        link->add_neighbor(neighbor); 
      }
    }

    Block* b = new Block;                // create a new block
    // init the block, related to the gid

    // for (unsigned j = 0; j < 3; ++j) {
    //   std::string ele = getID(gid, j); 
    //   b->add(ele); 
    //   b->ele2gid[ele] = gid; 

    //   int ngid_0 = (gid - 1 + nblocks) % nblocks; 
    //   std::string related_ele_0 = getID(ngid_0, j); 
    //   b->related_elements[ele].push_back(related_ele_0); 
    //   if(!b->has(related_ele_0)) {
    //     b->ele2gid[related_ele_0] = -1; 
    //   }

    //   int ngid_1 = (gid + 1) % nblocks; 
    //   std::string related_ele_1 = getID(ngid_1, j); 
    //   b->related_elements[ele].push_back(related_ele_1); 
    //   if(!b->has(related_ele_1)) {
    //     b->ele2gid[related_ele_1] = -1; 
    //   }
    // }

    if(gid == 0) {
      std::string eles[3] = {"0", "2", "3"}; 
      for(int i = 0; i < 3; ++i) {
        b->add(eles[i]); 
        b->ele2gid[eles[i]] = gid;   
      }
      b->related_elements["0"].push_back("1"); 
      b->related_elements["0"].push_back("7"); 
      b->related_elements["2"].push_back("6"); 
      b->related_elements["3"].push_back("1");
      b->related_elements["3"].push_back("4"); 

      b->ele2gid["1"] = -1; 
      b->ele2gid["7"] = -1; 
      b->ele2gid["6"] = -1; 
      b->ele2gid["4"] = -1; 
    }

    if(gid == 1) {
      std::string eles[3] = {"1", "4", "5"}; 
      for(int i = 0; i < 3; ++i) {
        b->add(eles[i]); 
        b->ele2gid[eles[i]] = gid;   
      }
      b->related_elements["1"].push_back("0"); 
      b->related_elements["1"].push_back("3"); 
      b->related_elements["4"].push_back("3"); 
      b->related_elements["4"].push_back("8"); 
      b->related_elements["5"].push_back("9"); 
      b->related_elements["5"].push_back("11"); 

      b->ele2gid["0"] = -1; 
      b->ele2gid["3"] = -1; 
      b->ele2gid["8"] = -1; 
      b->ele2gid["9"] = -1; 
      b->ele2gid["11"] = -1; 
    }

    if(gid == 2) {
      std::string eles[3] = {"6", "7", "10"}; 
      for(int i = 0; i < 3; ++i) {
        b->add(eles[i]); 
        b->ele2gid[eles[i]] = gid;   
      }
      b->related_elements["6"].push_back("2"); 
      b->related_elements["7"].push_back("0"); 
      b->related_elements["10"].push_back("11"); 

      b->ele2gid["2"] = -1; 
      b->ele2gid["0"] = -1; 
      b->ele2gid["11"] = -1; 
    }

    if(gid == 3) {
      std::string eles[3] = {"8", "9", "11"}; 
      for(int i = 0; i < 3; ++i) {
        b->add(eles[i]); 
        b->ele2gid[eles[i]] = gid;   
      }
      b->related_elements["8"].push_back("4"); 
      b->related_elements["9"].push_back("5"); 
      b->related_elements["11"].push_back("5"); 
      b->related_elements["11"].push_back("10");

      b->ele2gid["4"] = -1; 
      b->ele2gid["5"] = -1; 
      b->ele2gid["10"] = -1;  
    }

    // Deadlock case

    // if(gid == 2) {
    //   std::string eles[3] = {"6", "7", "11"}; 
    //   for(int i = 0; i < 3; ++i) {
    //     b->add(eles[i]); 
    //     b->ele2gid[eles[i]] = gid;   
    //   }
    //   b->related_elements["6"].push_back("2"); 
    //   b->related_elements["7"].push_back("0"); 
    //   b->related_elements["11"].push_back("5"); 
    //   b->related_elements["11"].push_back("10");
      
    //   b->ele2gid["2"] = -1; 
    //   b->ele2gid["0"] = -1; 
    //   b->ele2gid["5"] = -1; 
    //   b->ele2gid["10"] = -1; 
    // }

    // if(gid == 3) {
    //   std::string eles[3] = {"8", "9", "10"}; 
    //   for(int i = 0; i < 3; ++i) {
    //     b->add(eles[i]); 
    //     b->ele2gid[eles[i]] = gid;   
    //   }
    //   b->related_elements["8"].push_back("4"); 
    //   b->related_elements["9"].push_back("5"); 
    //   b->related_elements["10"].push_back("11"); 

    //   b->ele2gid["4"] = -1; 
    //   b->ele2gid["5"] = -1; 
    //   b->ele2gid["11"] = -1;  
    // }

    master.add(gid, b, link); 
  }

  master.foreach(&query_gid);
  // master.execute(); 
  master.exchange();
  master.foreach(&answer_gid);  
  master.exchange();
  master.foreach(&save_gid);

  if(IEXCHANGE) { // for iexchange
    
    iexchange_process(master); 
  } else { // for exchange
    bool all_done = false;
    while(!all_done) {
      exchange_process(master); 

      // int total_changes;
      // for (int i = 0; i < master.size(); i++)
      //     total_changes = master.proxy(i).get<size_t>();

      int total_changes = master.proxy(master.loaded_block()).read<int>();
      // int total_changes = master.proxy(master.loaded_block()).get<int>();

      std::cout<<total_changes<<std::endl; 

      all_done = total_changes == 0;
    }
  }

  // master.iexchange(&union_find_iexchange, 16, 1000);
}

