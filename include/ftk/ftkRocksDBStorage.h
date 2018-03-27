#ifndef _FTK_ROCKSDB_STORAGE
#define _FTK_ROCKSDB_STORAGE

#include <rocksdb/db.h>

class ftkRocksDBStorage : public ftkStorage {
public: 
  void open(void *p) {
    _db = static_cast<rocksdb::DB*>(p);
    _external_db = true;
  }

  void open(const std::string& dbname) {
    rocksdb::Options options;
    options.create_if_missing = true;
    options.compression = rocksdb::kBZip2Compression;
    rocksdb::Status status = rocksdb::DB::Open(options, dbname.c_str(), &db);
    assert(status.ok());
  }

  void close() {
    if (_external_db == false)
      delete _db;
  }

  void put(const std::string& key, const std::string& val) {
    _db->Put(rocksdb::WriteOptions(), key, val);
  }

  std::string get(const std::string& key) {
    _db->Get(rocksdb::ReadOptions(), key, val);
  }

private:
  rocksdb::DB *_db;
  bool _external_db = false;
};

#endif
