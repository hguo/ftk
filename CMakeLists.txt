cmake_minimum_required (VERSION 3.1.3)
project (ftk)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory (tests)
add_subdirectory (examples)
