#ifndef _FTK_STORAGE
#define _FTK_STORAGE

#include <iostream>

class ftkStorage {
public: 
  virtual void open(void*) = 0;
  virtual void open(const std::string&) = 0;
  virtual void close() = 0;

  virtual void put(const std::string& key, const std::string& val) = 0;
  virtual std::string get(const std::string& key) = 0;
};

#endif
