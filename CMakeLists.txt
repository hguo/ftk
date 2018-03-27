cmake_minimum_required (VERSION 3.1.3)
project (ftk)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)

option (FTK_USE_ROCKSDB "Build with RocksDB" OFF)
  
if (FTK_USE_ROCKSDB)
  find_package (RocksDB REQUIRED)
  include_directories (${RocksDB_INCLUDE_DIR})
endif ()

# add_subdirectory (src)
add_subdirectory (tests)
