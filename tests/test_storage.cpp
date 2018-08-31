#if FTK_USE_LEVELDB
#include "ftk/storage/leveldb.h"
#elif FTK_USE_ROCKSDB
#include "ftk/storage/rocksdb.h"
#else
#include "ftk/storage/native.h"
#endif

int main(int argc, char **argv) 
{
#if FTK_USE_LEVELDB
  std::shared_ptr<ftk::storage> s(new ftk::storage_leveldb);
#elif FTK_USE_ROCKSDB
  std::shared_ptr<ftk::storage> s(new ftk::storage_rocksdb);
#else
  std::shared_ptr<ftk::storage> s(new ftk::storage_native);
#endif

  s->open("/tmp/ftkdb");
  s->put("key", "val");

  return 0;
}
