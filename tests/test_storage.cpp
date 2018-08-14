#if FTK_USE_LEVELDB
#include "ftk/storage/leveldbStorage.h"
#elif FTK_USE_ROCKSDB
#include "ftk/storage/rocksdbStorage.h"
#else
#include "ftk/storage/dirStorage.h"
#endif

int main(int argc, char **argv) 
{
#if FTK_USE_LEVELDB
  ftk::Storage *store = new ftk::LevelDBStorage;
#elif FTK_USE_ROCKSDB
  ftk::Storage *store = new ftk::RocksDBStorage;
#else
  ftk::Storage *store = new ftk::DirStorage;
#endif
  store->open("/tmp/ftkdb");
  store->put("key", "val");
  delete store;

  return 0;
}
