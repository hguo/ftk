#ifndef _FTK_STORAGE
#define _FTK_STORAGE

#include <iostream>
#include "ftk/external/json.hh"

namespace ftk {

class storage {
public: 
  virtual ~storage() {}

  virtual bool open(void*) {return false;}
  virtual bool open(const std::string&) = 0;
  virtual void close() = 0;

  virtual void put(const std::string& key, const std::string& val) = 0;
  void put_json(const std::string& key, const nlohmann::json& j) {put(key, j.dump());}
  template <typename T> void put_obj(const std::string& key, const T& val) {
    nlohmann::json j;
    nlohmann::adl_serializer<T>::to_json(j, val);
    put(key, j.dump());
  } 

  virtual std::string get(const std::string& key) = 0;
};

}

#endif
