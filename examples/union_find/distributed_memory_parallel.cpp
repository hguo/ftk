#include <vector>
#include <iostream>

#include <ftk/external/diy/mpi.hpp>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/assigner.hpp>
#include <ftk/external/diy/serialization.hpp>

#include <ftk/basic/sparse_union_find.hh>

typedef std::pair<std::string, int> ele_gid;  // stores element and its global block id
typedef std::map<std::string, std::vector<std::string>> r_ele_map; 

typedef std::map<std::string, int> ele2gid_map; 

struct Block : public ftk::sparse_union_find<std::string> {
  // map element id to ids of its related elements
  r_ele_map related_elements; 
  ele2gid_map ele2gid; 

  Block(): sparse_union_find() { 
    
  }
};


// void*   create_block()                      { return new Block; }
// void    destroy_block(void* b)              { delete static_cast<Block*>(b); }
// void    save_block(const void* b,
//                    diy::BinaryBuffer& bb)   { diy::save(bb, *static_cast<const Block*>(b)); }
// void    load_block(void* b,
//                    diy::BinaryBuffer& bb)   { diy::load(bb, *static_cast<Block*>(b)); }


std::string getID(int gid, int j) {
  return std::to_string(gid) + '_' + std::to_string(j); 
}


// Serilization
// https://github.com/diatomic/diy/blob/master/examples/simple/block.h


void query_gid(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  // std::cout<<"Query: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 
  
  for(ele2gid_map::iterator it = b->ele2gid.begin(); it != b->ele2gid.end(); ++it) {
    std::string ele = it->first; 

    if(b->has(ele)) {
      it->second = gid; 
    } else {
      // std::cout<<ele<<std::endl; 

      // Seems can be optimized by Reduce or Gather
      for (int i = 0; i < l->size(); ++i) {
        cp.enqueue(l->target(i), ele);
      }
    }
  }

  std::cout<<std::endl; 
}

void answer_gid(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();
  // diy::Master* master = cp.master(); 

  // std::cout<<"Answer: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 

  while(!cp.empty_incoming_queues()) {
    std::vector<int> in; // gids of incoming neighbors in the link
    cp.incoming(in);

    // for all neighbor blocks
    // dequeue data received from this neighbor block in the last exchange
    for (unsigned i = 0; i < in.size(); ++i) {
      if(cp.incoming(in[i])) {
        std::string ele; 
        cp.dequeue(in[i], ele);

        if(b->has(ele)) {
          ele_gid _pair(ele, gid); 

          // std::cout<<ele<<"-";
          // std::cout<<gid;
          // std::cout<<std::endl; 

          cp.enqueue(l->target(l->find(in[i])), _pair);
          // cp.enqueue(l->target(l->find(in[i])), 1);
        }
      }
    }
  }

  std::cout<<std::endl; 
}

void save_gid(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 

  // std::cout<<"Save: "<<std::endl; 
  // std::cout<<"Block ID: "<<gid<<std::endl; 

  while(!cp.empty_incoming_queues()) {
    std::vector<int> in; // gids of incoming neighbors in the link
    cp.incoming(in);

    // for all neighbor blocks
    // dequeue data received from this neighbor block in the last exchange
    for (unsigned i = 0; i < in.size(); ++i) {
      if(cp.incoming(in[i])) {
        ele_gid _pair; 
        cp.dequeue(in[i], _pair);

        // std::cout<<_pair.first<<" - "<<_pair.second<<std::endl; 

        b->ele2gid[_pair.first] = _pair.second; 

        // std::cout<<"i-"<<i<<'-'<<in[i]<<std::endl;
        // int t; 
        // cp.dequeue(in[i], t); 
        // std::cout<<t<<std::endl; 
      }
    }
  }

  std::cout<<std::endl; 
}

void unite_once(Block* b, const diy::Master::ProxyWithLink& cp) {
  for(r_ele_map::iterator it = b->related_elements.begin(); it != b->related_elements.end(); ++it) {
    std::string ele = it->first; 
    if(b->is_root(ele)) {
      for(std::vector<std::string>::iterator it_vec = b->related_elements[ele].begin(); it_vec != b->related_elements[ele].end(); ++it_vec) {
        std::string related_ele = *it_vec; 
        if(related_ele > ele) {
          b->related_elements[ele].erase(it_vec); 

          b->set_parent(ele, related_ele); 

          break ; 
        }
      }
    }
  }
}

void query_grandparent(Block* b, const diy::Master::ProxyWithLink& cp) {
  diy::Link* l = cp.link();

  for(std::vector<std::string>::iterator it = b->eles.begin(); it != b->eles.end(); ++it) {
    std::string ele = *it; 
    if(!b->is_root(ele)) {
      std::string parent = b->parent(ele); 

      if(b->has(parent)) {
        std::string grandparent = b->parent(parent); 
        b->set_parent(ele, grandparent); 
      } else {
        int gid = b->ele2gid[parent]; 
        cp.enqueue(l->target(l->find(gid)), ele);
      }
    }
  }
}

void answer_grandparent(Block* b, const diy::Master::ProxyWithLink& cp) {
  diy::Link* l = cp.link();
  std::vector<int> in; // gids of incoming neighbors in the link
  cp.incoming(in);

  // for all neighbor blocks
  // dequeue data received from this neighbor block in the last exchange
  for (unsigned i = 0; i < in.size(); ++i) {
    if(cp.incoming(in[i])) {
      std::string ele; 
      cp.dequeue(in[i], ele);

      if(b->has(ele) && !b->is_root(ele)) {
        std::pair<std::string, std::string> _pair(ele, b->parent(ele)); 
        cp.enqueue(l->target(l->find(in[i])), _pair);
      }
    }
  }
}

void compress_path(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  std::vector<int> in; // gids of incoming neighbors in the link
  cp.incoming(in);

  // for all neighbor blocks
  // dequeue data received from this neighbor block in the last exchange

  std::cout<<"Compress Path: "<<std::endl; 
  std::cout<<"Block ID: "<<gid<<std::endl; 

  for (unsigned i = 0; i < in.size(); ++i) {
    if(cp.incoming(in[i])) {
      std::pair<std::string, std::string> _pair; 
      cp.dequeue(in[i], _pair); 

      std::cout<<_pair.first<<" - "<<" - "<<b->parent(_pair.first)<<_pair.second<<std::endl; 

      b->set_parent(_pair.first, _pair.second); 
    }
  }
}


// void pass_unions(Block* b, const diy::Master::ProxyWithLink& cp) {
//   auto src = &b->related_elements[ele]; 
//         if(b->has(related_ele)) {
//             // Update directly, if the realted element is in the block

//             auto dest = &b->related_elements[related_ele]; 
//             dest.insert(
//               dest.end(),
//               std::make_move_iterator(src.begin()),
//               std::make_move_iterator(src.end())
//             );
//         } else {
//           // Otherwise, communicate with other processes
//           //...
//           ...;
//         }


//   if(!b->pending_unions.empty()) {
//     std::pair<std::string, std::string> _union = b->pending_unions.pop_back(); 

//     // if two elements are in the same block
//     if(!b->find(_union.first).empty() && !b->find(_union.second).empty()) {
//       b->unite(_union.first, _union.second); 
//     } else {
//       std::string root_first = b->find(_union.first); 


//       if(root_first < )
//     }
//   }

//   // cp.collectives()->clear();
//   // cp.all_reduce(done, std::logical_and<int>());
// }


// void get_unions(Block* b, const diy::Master::ProxyWithLink& cp) {

// }

// void is_done(Block* b, const diy::Master::ProxyWithLink& cp) {
//   // cp.all_reduce...
// }

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
    for (unsigned j = 0; j < 3; ++j) {
      std::string ele = getID(gid, j); 
      b->add(ele); 
      b->ele2gid[ele] = gid; 

      int ngid_0 = (gid - 1 + nblocks) % nblocks; 
      std::string related_ele_0 = getID(ngid_0, j); 
      b->related_elements[ele].push_back(related_ele_0); 
      if(!b->has(related_ele_0)) {
        b->ele2gid[related_ele_0] = -1; 
      }

      int ngid_1 = (gid + 1) % nblocks; 
      std::string related_ele_1 = getID(ngid_1, j); 
      b->related_elements[ele].push_back(related_ele_1); 
      if(!b->has(related_ele_1)) {
        b->ele2gid[related_ele_1] = -1; 
      }
    }

    master.add(gid, b, link); 
  }

  master.foreach(&query_gid);
  master.exchange();
  master.foreach(&answer_gid);  
  master.exchange();
  master.foreach(&save_gid);




  // bool all_done = false;
  // while (!all_done) {
  //   // compute, exchange, compute
    
    master.foreach(&unite_once);
    master.foreach(&query_grandparent);
    master.exchange();                 
    master.foreach(&answer_grandparent);
    master.exchange();
    master.foreach(&compress_path);
    // master.foreach(&pass_unions);
    // master.exchange();
    // master.foreach(&get_unions);

    // master.foreach(&is_done); // can use a reduce all to check whether every thing is done. 
    // master.exchange();

  //   all_done = master.proxy(master.loaded_block()).read<bool>();
  // }

  // if (world.rank() == 0)
  //   std::cout << "Total iterations: " << master.block<Block>(master.loaded_block())->count << std::endl;
}

