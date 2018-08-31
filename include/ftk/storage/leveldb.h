#ifndef _FTK_LEVELDB_STORAGE
#define _FTK_LEVELDB_STORAGE

#include "ftk/storage/base.h"
#include <leveldb/db.h>

namespace ftk {

class storage_leveldb : public storage {
public: 
  bool open(void *p) {
    _db = static_cast<leveldb::DB*>(p);
    _external_db = true;
    return true;
  }

  bool open(const std::string& dbname) {
    leveldb::Options options;
    options.create_if_missing = true;
    leveldb::Status status = leveldb::DB::Open(options, dbname.c_str(), &_db);
    return status.ok();
    // assert(status.ok());
  }

  void close() {
    if (_external_db == false)
      delete _db;
  }

  void put(const std::string& key, const std::string& val) {
    _db->Put(leveldb::WriteOptions(), key, val);
  }

  std::string get(const std::string& key) {
    std::string val;
    _db->Get(leveldb::ReadOptions(), key, &val);
    return val;
  }

private:
  leveldb::DB *_db;
  bool _external_db = false;
};

}

#endif
