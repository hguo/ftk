#ifndef _FTK_STORAGE
#define _FTK_STORAGE

#include <iostream>

#if FTK_USE_ROCKSDB
#include <rocksdb/db.h>
#endif

class ftkStorage {
public: 
  virtual void open(void*) = 0;
  virtual void open(const std::string&) = 0;
  virtual void close() = 0;

  virtual void put(const std::string& key, const std::string& val) = 0;
  virtual std::string get(const std::string& key) = 0;
};

class ftkStorageLevelDB : public ftkStorage {
  // TODO
};

class ftkStorageRocksDB : public ftkStorage {
public: 
  void open(void *p) {
#if FTK_USE_ROCKSDB
    _db = static_cast<rocksdb::DB*>(p);
    _external_db = true;
#else
    assert(false);
#endif
  }

  void open(const std::string& dbname) {
#if FTK_USE_ROCKSDB
    rocksdb::Options options;
    options.create_if_missing = true;
    options.compression = rocksdb::kBZip2Compression;
    rocksdb::Status status = rocksdb::DB::Open(options, dbname.c_str(), &db);
    assert(status.ok());
#else
    assert(false);
#endif
  }

  void close() {
#if FTK_USE_ROCKSDB
    if (_external_db == false)
      delete _db;
#endif
  }

  void put(const std::string& key, const std::string& val) {
#if FTK_USE_ROCKSDB
    _db->Put(rocksdb::WriteOptions(), key, val);
#else 
    assert(false);
#endif
  }

  std::string get(const std::string& key) {
#if FTK_USE_ROCKSDB
    // _db->Get(rocksdb::ReadOptions(), key, val);
#else
    assert(false);
    return "";
#endif
  }

private:
#if FTK_USE_ROCKSDB
  rocksdb::DB *_db;
  bool _external_db = false;
#endif
};

#endif
