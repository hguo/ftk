#include <vector>
#include <iostream>
#include <thread>

#include <ftk/basic/shared_union_find.hh>


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
