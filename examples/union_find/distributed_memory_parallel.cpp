// Serilization
// https://github.com/diatomic/diy/blob/master/examples/simple/block.h


// Probibly, can be optimized by "An O(logn) parallel connectivity algorithm"


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
  Block(): nchanges(0), sparse_union_find() { 
    
  }

public: 
  // map element id to ids of its related elements
  int nchanges;
  r_ele_map related_elements; 
  ele2gid_map ele2gid; 
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

  // std::cout<<std::endl; 
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

  // std::cout<<std::endl; 
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

  // std::cout<<std::endl; 
}

void unite_once(Block* b, const diy::Master::ProxyWithLink& cp) {
  for(auto it = b->eles.begin(); it != b->eles.end(); ++it) {
    std::string ele = *it; 
    if(b->is_root(ele)) {
      for(std::vector<std::string>::iterator it_vec = b->related_elements[ele].begin(); it_vec != b->related_elements[ele].end(); ++it_vec) {
        std::string related_ele = *it_vec; 
        // if(related_ele > ele) {
        std::cout<<related_ele<<" "<<ele<<std::endl; 
        std::cout<<related_ele.empty()<<std::endl; 

        if(std::stoi(related_ele) > std::stoi(ele)) {
          b->related_elements[ele].erase(it_vec); 

          if(related_ele.empty()) {
            std::cout<<"111111111111"<<std::endl; 
          }

          b->set_parent(ele, related_ele); 
          b->nchanges += 1;
          std::cout<<ele<<" -> "<<related_ele<<std::endl; 

          break ; 
        }
      }
    }
  }
}

// Also locally update parent to grandparent by path compression
void query_grandparent(Block* b, const diy::Master::ProxyWithLink& cp) {
  diy::Link* l = cp.link();

  // std::cout<<"Query Grandparent: "<<std::endl; 
  // std::cout<<"Block ID: "<<cp.gid()<<std::endl; 

  for(std::vector<std::string>::iterator it = b->eles.begin(); it != b->eles.end(); ++it) {
    std::string ele = *it; 
    if(!b->is_root(ele)) {
      std::string parent = b->parent(ele); 

      if(parent.empty()) {
        std::cout<<"Error!"<<std::endl; 
      }

      if(b->has(parent)) {
        std::string grandparent = b->parent(parent); 

        if(grandparent.empty()) {
          std::cout<<"22222222222222 "<<parent<<std::endl; 
        }

        b->set_parent(ele, grandparent); 
        b->nchanges += 1;
        std::cout<<ele<<" -> "<<grandparent<<std::endl; 

      } else {
        int gid = b->ele2gid[parent]; 
        std::pair<std::string, std::string> send_pair(ele, parent); 
        cp.enqueue(l->target(l->find(gid)), send_pair);

        // std::cout<<ele<<" - "<<parent<<" - "<<gid<<std::endl; 
      }
    }
  }

  // std::cout<<std::endl; 
}

void answer_grandparent(Block* b, const diy::Master::ProxyWithLink& cp) {
  diy::Link* l = cp.link();
  // std::cout<<"Answer Grandparent: "<<std::endl; 
  // std::cout<<"Block ID: "<<cp.gid()<<std::endl; 

  while(!cp.empty_incoming_queues()) {
    std::vector<int> in; // gids of incoming neighbors in the link
    cp.incoming(in);

    // for all neighbor blocks
    // dequeue data received from this neighbor block in the last exchange
    for (unsigned i = 0; i < in.size(); ++i) {
      if(cp.incoming(in[i])) {
        std::pair<std::string, std::string> rec_pair; 
        cp.dequeue(in[i], rec_pair);
        std::string* ele_ptr = &rec_pair.first; 
        std::string* parent_ptr = &rec_pair.second; 

        if(b->has(*parent_ptr) && !b->is_root(*parent_ptr)) {
          std::string grandparent = b->parent(*parent_ptr); 

          if(grandparent.empty()) {
            std::cout<<"Error"<<std::endl; 
          }

          std::tuple<std::string, std::string, int> send_tuple(*ele_ptr, grandparent, b->ele2gid[grandparent]); 
          cp.enqueue(l->target(l->find(in[i])), send_tuple);

          // std::cout<<*ele_ptr<<" - "<<*parent_ptr<<" - "<<grandparent<<" - "<<b->ele2gid[grandparent]<<std::endl; 
        }
      }
    }
  }

  // std::cout<<std::endl; 
}

void compress_path(Block* b, const diy::Master::ProxyWithLink& cp) {
  // std::cout<<"Compress Path: "<<std::endl; 
  // std::cout<<"Block ID: "<<cp.gid()<<std::endl; 

  while(!cp.empty_incoming_queues()) {
    std::vector<int> in; // gids of incoming neighbors in the link
    cp.incoming(in);

    // for all neighbor blocks
    // dequeue data received from this neighbor block in the last exchange

    for (unsigned i = 0; i < in.size(); ++i) {
      if(cp.incoming(in[i])) {
        std::tuple<std::string, std::string, int> rec_tuple; 
        cp.dequeue(in[i], rec_tuple); 
        std::string* ele_ptr = &std::get<0>(rec_tuple); 
        std::string* grandpar_ptr = &std::get<1>(rec_tuple); 
        int* grandpar_gid_ptr = &std::get<2>(rec_tuple); 

        // std::cout<<*ele_ptr << " - " << b->parent(*ele_ptr) <<" - "<< *grandpar_ptr<<" - "<<*grandpar_gid_ptr<<std::endl; 

        if(grandpar_ptr->empty()) {
          std::cout<<"3333333333"<<std::endl; 
        }

        b->set_parent(*ele_ptr, *grandpar_ptr); 
        b->nchanges += 1;
        std::cout<<*ele_ptr<<" -> "<<*grandpar_ptr<<std::endl; 

        b->ele2gid[*grandpar_ptr] = *grandpar_gid_ptr; 
      }
    }
  }

  // std::cout<<std::endl; 
}


void pass_unions(Block* b, const diy::Master::ProxyWithLink& cp) {
  diy::Link* l = cp.link();

  for(auto ite_ele = b->eles.begin(); ite_ele != b->eles.end(); ++ite_ele) {
    std::string ele = *ite_ele; 
    auto src = &b->related_elements[ele]; 

    // Can be optimized by only its parent is a root, then do unions' passing

    if(src->size() > 0) {
      std::string par = b->parent(ele); 
      if(par.empty()) {
        std::cout<<"Error! "<<std::endl; 
      }


      if(b->has(par)) {
          // Update directly, if the realted element is in the block

          auto dest = &b->related_elements[par]; 
          dest->insert(
            dest->end(),
            std::make_move_iterator(src->begin()),
            std::make_move_iterator(src->end())
          );
      } else {
        // Otherwise, communicate with other processes
        int gid = b->ele2gid[par]; 

        for(auto ite_related_ele = src->begin(); ite_related_ele != src->end(); ++ite_related_ele) {
          std::tuple<std::string, std::string, std::string, int> send_tuple(ele, par, *ite_related_ele, b->ele2gid[*ite_related_ele]); 

          cp.enqueue(l->target(l->find(gid)), send_tuple); 
        }
      }
    }
  }
}

void update_unions2block(Block* b, std::string* ele_ptr, std::string* par_ptr, std::string* related_ele_ptr) {
  for(auto ite = b->related_elements[*related_ele_ptr].begin(); ite != b->related_elements[*related_ele_ptr].end(); ++ite) {

    if(*ite == *ele_ptr) {

      if(par_ptr->empty()) {
        std::cout<<"Error! "<<std::endl; 
      }

      *ite = *par_ptr; 

      break ; 
    }
  }
}

void save_unions(Block* b, const diy::Master::ProxyWithLink& cp) {
  int gid = cp.gid(); 
  diy::Link* l = cp.link();

  std::cout<<"Save Unions: "<<std::endl; 
  std::cout<<"Block ID: "<<cp.gid()<<std::endl; 

  while(!cp.empty_incoming_queues()) {
    std::vector<int> in; // gids of incoming neighbors in the link
    cp.incoming(in);

    // for all neighbor blocks
    // dequeue data received from this neighbor block in the last exchange

    for (unsigned i = 0; i < in.size(); ++i) {
      if(cp.incoming(in[i])) {
        std::tuple<std::string, std::string, std::string, int> rec_tuple; 
        cp.dequeue(in[i], rec_tuple); 

        std::string* ele_ptr = &std::get<0>(rec_tuple); 
        std::string* par_ptr = &std::get<1>(rec_tuple); 
        std::string* related_ele_ptr = &std::get<2>(rec_tuple); 
        int* gid_ptr = &std::get<3>(rec_tuple); 

        std::cout<<*ele_ptr<<" - "<<*par_ptr << " - " << *related_ele_ptr <<" - "<< *gid_ptr<<std::endl; 

        if(related_ele_ptr->empty()) {
          std::cout<<"Error! "<<std::endl; 
        }

        b->related_elements[*par_ptr].push_back(*related_ele_ptr); 
        b->ele2gid[*related_ele_ptr] = *gid_ptr; 

        // tell related elements, the end point has changed to its parent
        std::tuple<std::string, std::string, std::string> send_tuple(*ele_ptr, *par_ptr, *related_ele_ptr); 

        if(*gid_ptr == gid) {
          update_unions2block(b, ele_ptr, par_ptr, related_ele_ptr); 
        } else {
          cp.enqueue(l->target(l->find(*gid_ptr)), send_tuple); 
        }
      }
    }
  }

  std::cout<<std::endl; 
}

// update unions of related elements
void update_unions(Block* b, const diy::Master::ProxyWithLink& cp) {
  std::cout<<"Update Unions: "<<std::endl; 
  std::cout<<"Block ID: "<<cp.gid()<<std::endl; 

  while(!cp.empty_incoming_queues()) {
    std::vector<int> in; // gids of incoming neighbors in the link
    cp.incoming(in);

    // for all neighbor blocks
    // dequeue data received from this neighbor block in the last exchange

    for (unsigned i = 0; i < in.size(); ++i) {
      if(cp.incoming(in[i])) {
        std::tuple<std::string, std::string, std::string> rec_tuple; 
        cp.dequeue(in[i], rec_tuple); 

        std::string* ele_ptr = &std::get<0>(rec_tuple); 
        std::string* par_ptr = &std::get<1>(rec_tuple); 
        std::string* related_ele_ptr = &std::get<2>(rec_tuple); 

        // std::cout<<*ele_ptr<<" - "<<*par_ptr << " - " << *related_ele_ptr<<std::endl; 

        update_unions2block(b, ele_ptr, par_ptr, related_ele_ptr); 
      }
    }
  }

  std::cout<<std::endl; 
}

void total_changes(Block* b, const diy::Master::ProxyWithLink& cp) {
  // for(auto ite = b->related_elements.begin(); ite != b->related_elements.end(); ++ite) {
  //   std::string ele = ite->first; 
  //   for(auto ite_related_ele = ite->second.begin(); ite_related_ele != ite->second.end(); ++ite_related_ele) {
  //     if(*ite_related_ele == ele) {
  //       ite->second.erase(ite_related_ele); 
  //     }
  //   }
  // }

  cp.collectives()->clear();
  cp.all_reduce(b->nchanges, std::plus<int>()); 
  b->nchanges = 0;
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

    master.add(gid, b, link); 
  }

  master.foreach(&query_gid);
  master.exchange();
  master.foreach(&answer_gid);  
  master.exchange();
  master.foreach(&save_gid);

  bool all_done = false;
  while (!all_done) {
    master.foreach(&unite_once);
    master.foreach(&query_grandparent);
    master.exchange();                 
    master.foreach(&answer_grandparent);
    master.exchange();
    master.foreach(&compress_path);
    master.foreach(&pass_unions);
    master.exchange();
    master.foreach(&save_unions);
    master.exchange();
    master.foreach(&update_unions);

    master.foreach(&total_changes); // can use a reduce all to check whether every thing is done. 
    master.exchange();

    int total_changes = master.proxy(master.loaded_block()).read<int>();
    all_done = total_changes == 0;

    std::cout<<total_changes<<"==========================="<<std::endl; 
  }

  // if (world.rank() == 0)
  //   std::cout << "Total iterations: " << master.block<Block>(master.loaded_block())->count << std::endl;
}

