// Serilization
// https://github.com/diatomic/diy/blob/master/examples/simple/block.h


// Probibly, can be optimized by "An O(logn) parallel connectivity algorithm"

#include <vector>
#include <iostream>

#include <ftk/external/diy/mpi.hpp>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/assigner.hpp>
// #include <ftk/external/diy/serialization.hpp>

#include <ftk/basic/distributed_union_find.hh>

// typedef std::pair<std::string, int> ele_gid;  // stores element and its global block id
// typedef std::map<std::string, std::vector<std::string>> r_ele_map; 
// typedef std::map<std::string, int> ele2gid_map; 

// void*   create_block()                      { return new Block; }
// void    destroy_block(void* b)              { delete static_cast<Block*>(b); }
// void    save_block(const void* b,
//                    diy::BinaryBuffer& bb)   { diy::save(bb, *static_cast<const Block*>(b)); }
// void    load_block(void* b,
//                    diy::BinaryBuffer& bb)   { diy::load(bb, *static_cast<Block*>(b)); }


// std::string getID(int gid, int j) {
//   return std::to_string(gid) + '_' + std::to_string(j); 
// }

void get_synthetic_data(std::vector<Block*>& blocks) {
  int ngrid = 4; 

  // for the local blocks in this processor
  for (unsigned i = 0; i < ngrid; ++i) {

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

    if(i == 0) {
      std::string eles[3] = {"0", "2", "3"}; 
      for(int j = 0; j < 3; ++j) {
        b->add(eles[j]); 
        // b->ele2gid[eles[j]] = i;   
      }

      b->add_related_element("0", "1"); 
      b->add_related_element("0", "7");
      b->add_related_element("2", "6"); 
      b->add_related_element("3", "1"); 
      b->add_related_element("3", "4"); 
    }

    if(i == 1) {
      std::string eles[3] = {"1", "4", "5"}; 
      for(int j = 0; j < 3; ++j) {
        b->add(eles[j]); 
        // b->ele2gid[eles[j]] = i;   
      }

      b->add_related_element("1", "0"); 
      b->add_related_element("1", "3");
      b->add_related_element("4", "3"); 
      b->add_related_element("4", "8"); 
      b->add_related_element("5", "9"); 
      b->add_related_element("5", "11"); 
    }

    if(i == 2) {
      std::string eles[3] = {"6", "7", "10"}; 
      for(int j = 0; j < 3; ++j) {
        b->add(eles[j]); 
        // b->ele2gid[eles[j]] = i;   
      }

      b->add_related_element("6", "2"); 
      b->add_related_element("7", "0");
      b->add_related_element("10", "11"); 
    }

    if(i == 3) {
      std::string eles[3] = {"8", "9", "11"}; 
      for(int j = 0; j < 3; ++j) {
        b->add(eles[j]); 
        // b->ele2gid[eles[j]] = i;   
      }

      b->add_related_element("8", "4"); 
      b->add_related_element("9", "5");
      b->add_related_element("11", "5"); 
      b->add_related_element("11", "10"); 
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

    blocks.push_back(b); 
  }
}

int main(int argc, char* argv[]) {

  int                       nthreads = 1; 
  int                       memblocks = -1; //2; // number of blocks to store in memory

  diy::mpi::environment     env(argc, argv);
  // diy::mpi::environment     env(NULL, NULL);
  diy::mpi::communicator    world;


  std::vector<Block*> blocks;
  get_synthetic_data(blocks); 

  int                       nblocks = blocks.size();
  // int                       nblocks = 4*world.size();

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

  // get_connected_components
  std::vector<Block*> local_blocks; 
  std::vector<int> gids;                     // global ids of local blocks
  assigner.local_gids(world.rank(), gids);   // get the gids of local blocks for a given process rank 
  for (unsigned i = 0; i < gids.size(); ++i) {
    int gid = gids[i];
    local_blocks.push_back(blocks[gid]); 
  }

  exec_distributed_union_find(world, master, assigner, local_blocks); 

  std::vector<std::set<std::string>> ele_sets;
  get_sets(world, master, assigner, ele_sets); 
}

