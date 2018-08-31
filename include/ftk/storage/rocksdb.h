#ifndef _FTK_ROCKSDB_STORAGE
#define _FTK_ROCKSDB_STORAGE

#include "ftk/storage/base.h"
#include <rocksdb/db.h>

namespace ftk {

class storage_rocksdb : public storage {
public: 
  bool open(void *p) {
    _db = static_cast<rocksdb::DB*>(p);
    _external_db = true;
    return true;
  }

  bool open(const std::string& dbname) {
    rocksdb::Options options;
    options.create_if_missing = true;
    rocksdb::Status status = rocksdb::DB::Open(options, dbname.c_str(), &_db);
    // assert(status.ok());
    return status.ok();
  }

  void close() {
    if (_external_db == false)
      delete _db;
  }

  void put(const std::string& key, const std::string& val) {
    _db->Put(rocksdb::WriteOptions(), key, val);
  }

  std::string get(const std::string& key) {
    std::string val;
    _db->Get(rocksdb::ReadOptions(), key, &val);
    return val;
  }

private:
  rocksdb::DB *_db;
  bool _external_db = false;
};

}

#endif
