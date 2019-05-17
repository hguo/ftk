#include <vector>
#include <iostream>
#include <thread>

#include <ftk/external/diy/mpi.hpp>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy/assigner.hpp>
#include <ftk/external/diy/serialization.hpp>

#include <ftk/basic/sparse_union_find.hh>


// CPP program to demonstrate multithreading 
// using three different callables. 

// https://www.geeksforgeeks.org/multithreading-in-cpp/

using namespace std; 
  
// A dummy function 
void foo(int Z) 
{ 
    cout<<"1: "; 
    for (int i = 0; i < Z; i++) { 
        cout << "Thread using function"
               " pointer as callable\n"; 
    } 
} 
  
// A callable object 
class thread_obj { 
public: 
    void operator()(int x) 
    { 
        cout<<"2: "; 
        for (int i = 0; i < x; i++) 
            cout << "Thread using function"
                  " object as  callable\n"; 
    } 
}; 
  
int main() 
{ 
    cout << "Threads 1 and 2 and 3 "
         "operating independently" << endl; 
  
    // This thread is launched by using  
    // function pointer as callable 
    thread th1(foo, 3); 
  
    // This thread is launched by using 
    // function object as callable 
    thread th2(thread_obj(), 3); 
  
    // Define a Lambda Expression 
    auto f = [](int x) { 
        cout<<"3: "; 
        for (int i = 0; i < x; i++) 
            cout << "Thread using lambda"
             " expression as callable\n"; 
    }; 
  
    // This thread is launched by using  
    // lamda expression as callable 
    thread th3(f, 3); 
  
    // Wait for the threads to finish 
    // Wait for thread t1 to finish 
    th1.join(); 
  
    // Wait for thread t2 to finish 
    th2.join(); 
  
    // Wait for thread t3 to finish 
    th3.join(); 
  
    return 0; 
} 



// struct Block : public ftk::sparse_quick_union<std::string> {
//   // map element id to ids of its related elements
//   std::map<std::string, std::vector<std::string>> related_elements; 
//   int* t; 
//   int gid; 

//   Block(): sparse_quick_union() { 
    
//   }
// };

// std::string getID(int gid, int j) {
//   return std::to_string(gid) + '_' + std::to_string(j); 
// }


// void test(Block* b, const diy::Master::ProxyWithLink& cp) {
//   if(b->gid == 0) {
//     b->t[1] += 1;
//   }

//   std::cout << "test" << std::endl; 
//   std::cout << b->gid << std::endl; 
//   std::cout << b->t[1] << std::endl;  
// }


// int main(int argc, char* argv[]) {
//   diy::mpi::environment     env(argc, argv);
//   diy::mpi::communicator    world;
  
//   int                       nblocks = 2;
//   int                       nthreads = 1;


//   diy::FileStorage          storage("./DIY.XXXXXX");

//   diy::Master               master(world, nthreads);

//   srand(time(NULL));

//   diy::ContiguousAssigner        assigner(world.size(), nblocks);
//   // diy::RoundRobinAssigner     assigner(world.size(), nblocks);
//   // diy::DynamicAssigner        assigner(world.size(), nblocks);

//   std::vector<int> gids;                     // global ids of local blocks
//   assigner.local_gids(world.rank(), gids);   // get the gids of local blocks for a given process rank 

//   int t[2] = {1, 2}; 

//   // for the local blocks in this processor
//   for (unsigned i = 0; i < gids.size(); ++i) {
//     int gid = gids[i];

//     diy::Link*    link = new diy::Link; 

//     Block* b = new Block;                // create a new block
//     b->t = t; 
//     b->gid = gid; 
    
//     for (unsigned j = 0; j < 3; ++j) {
//       std::string ele = std::to_string(j); 
//       b->add(ele); 
//     }

//     master.add(gid, b, link); 
//   }

//   bool is_done = false;
//   // while (!is_done) {
//   master.foreach(&test);          // callback function executed on each local block
//   // master.exchange(); 
//   // master.foreach(&test);

//     // is_done = master.block(master.loaded_block())->;
//   // }

//   // if (world.rank() == 0)
//   //   std::cout << "Total iterations: " << master.block<Block>(master.loaded_block())->count << std::endl;
// }

