#ifndef _FTK_INTERVAL_HH
#define _FTK_INTERVAL_HH

namespace ftk {

template <typename T>
struct interval {
  T lb = std::numeric_limits<T>::infinity(), 
    ub = std::numeric_limits<T>::infinity();
  
};

}

#endif
