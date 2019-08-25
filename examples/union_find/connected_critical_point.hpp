#include <ftk/basic/distributed_union_find.hh>

struct intersection_t {
  template <class Archive> void serialize(Archive & ar) {
    ar(eid, x[0], x[1], x[2], val);
  }

  float&  operator[](unsigned i)                          { return corner[i]; }
  float   operator[](unsigned i) const                    { return corner[i]; }

  float x[3]; // the spacetime coordinates of the trajectory
  float corner[3]; // the spacetime coordinates of the left corner of the element

  float val; // scalar value at the intersection
  
  std::string eid; // element id
  std::set<std::string> related_elements; 
};

namespace diy
{
  template<>
  struct Serialization<intersection_t>
  {
      static void save(BinaryBuffer& bb, const intersection_t& msg)
      {
          diy::save(bb, msg.x);
          diy::save(bb, msg.corner);
          diy::save(bb, msg.val);
          diy::save(bb, msg.eid);
          diy::save(bb, msg.related_elements);
      }

      static void load(BinaryBuffer& bb, intersection_t& msg)
      {
          diy::load(bb, msg.x);
          diy::load(bb, msg.corner);
          diy::load(bb, msg.val);
          diy::load(bb, msg.eid);
          diy::load(bb, msg.related_elements);
      }
  };
}

// ==========================================================

// DIY Block for distributed critical point tracking
struct Block_Critical_Point : public Block_Union_Find {
  Block_Critical_Point() : intersections(), Block_Union_Find() {

  }


  void get_sets(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results);

public:
  std::map<std::string, intersection_t> intersections; // Points on trajectories of critical points
  std::vector<intersection_t> points; 

};

// ==========================================================

struct Message_Critical_Point : public Message_Union_Find {
  Message_Critical_Point() {

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

  void receive_intersection(std::pair<std::string, intersection_t>& pair) {
    pair.first = strs[0]; 

    intersection_t& I = pair.second; 
    I.x[0] = std::stof(strs[1]); 
    I.x[1] = std::stof(strs[2]); 
    I.x[2] = std::stof(strs[3]); 
    I.val = std::stof(strs[4]); 
    I.eid = strs[5];
  }

};

namespace diy
{
  template<>
  struct Serialization<Message_Critical_Point>
  {
      static void save(BinaryBuffer& bb, const Message_Critical_Point& msg)
      {
          diy::save(bb, msg.tag);
          // diy::save(bb, msg.str);
          diy::save(bb, msg.strs);
      }

      static void load(BinaryBuffer& bb, Message_Critical_Point& msg)
      {
          diy::load(bb, msg.tag);
          // diy::load(bb, msg.str);
          diy::load(bb, msg.strs);
      }
  };
}

// ==========================================================


// Generate sets of elements

// Method 1:
  // Send to p0
void send_2_p0(Block_Critical_Point* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  if(gid != 0) {
    // std::vector<std::pair<std::string, std::string>> local_pairs; 

    auto& target = l->target(l->find(0));
    for(auto& ele : b->eles) {
      std::string parent = b->parent(ele); 

      // local_pairs.push_back(std::make_pair(ele, parent)); 

      std::pair<std::string, std::string> local_pair(ele, parent); 

      Message_Critical_Point send_msg; 
      send_msg.send_ele_parent_pair(local_pair); 

      cp.enqueue(target, send_msg); 
    }

    for(auto& pair : b->intersections) {
      Message_Critical_Point send_msg; 
      send_msg.send_intersection(pair); 

      cp.enqueue(target, send_msg); 
    }
  }
}


// Gather all element-parent information to the first block
bool gather_2_p0(Block_Critical_Point* b, const diy::Master::ProxyWithLink& cp) {
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
          Message_Critical_Point msg; 
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
void get_sets_on_p0(Block_Critical_Point* b, diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {
  master.foreach(&send_2_p0); 
  master.iexchange(&gather_2_p0); 

  if(world.rank() == 0) {

    std::map<std::string, std::set<std::string>> root2set; 

    // #ifdef FTK_HAVE_MPI
    //   std::cout<<"# of elements on proc. " << world.rank() <<" : "<< b->eles.size()<<std::endl; 
    // #endif

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
void send_2_roots(Block_Critical_Point* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  for(auto& ele : b->eles) {
    std::string root = b->parent(ele); 
    if(b->has(root)) continue ; 

    int gid_root = b->get_gid(root); 
    auto& target = l->target(l->find(gid_root)); 

    Message_Critical_Point send_msg_intersection; 
    send_msg_intersection.send_intersection(*(b->intersections.find(ele))); 

    cp.enqueue(target, send_msg_intersection); 
  }
}

bool gather_on_roots(Block_Critical_Point* b, const diy::Master::ProxyWithLink& cp) {

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
        Message_Critical_Point msg; 
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
void get_sets_on_roots(Block_Critical_Point* b, diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {
  master.foreach(&send_2_roots); 
  master.iexchange(&gather_on_roots); 

  // std::vector<int> gids;                     // global ids of local blocks
  // assigner.local_gids(world.rank(), gids);

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
  // Step one: send to redistributed processes
  // Step two: gather on redistributed processes

void send_2_redistributed_processes(Block_Critical_Point* b, const diy::Master::ProxyWithLink& cp) {
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

      Message_Critical_Point send_msg; 
      send_msg.send_ele_parent_pair(local_pair); 

      cp.enqueue(target, send_msg); 

      // ==============

      Message_Critical_Point send_msg_intersection; 
      send_msg_intersection.send_intersection(*(b->intersections.find(ele))); 

      cp.enqueue(target, send_msg_intersection); 

      b->erase_element(ele);  
    }    
  }
}

bool gather_on_redistributed_processes(Block_Critical_Point* b, const diy::Master::ProxyWithLink& cp) {

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
        Message_Critical_Point msg; 
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
void get_sets_redistributed(Block_Critical_Point* b, diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {

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
inline void Block_Critical_Point::get_sets(diy::mpi::communicator& world, diy::Master& master, diy::ContiguousAssigner& assigner, std::vector<std::set<std::string>>& results) {
  // get_sets_on_p0(this, world, master, assigner, results); 
  // get_sets_on_roots(this, world, master, assigner, results); 
  get_sets_redistributed(this, world, master, assigner, results); 
}