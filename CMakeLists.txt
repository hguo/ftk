cmake_minimum_required (VERSION 3.1.3)
project (ftk)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/hypermesh/include)

option (FTK_BUILD_TESTS "Build the test programs." OFF)
option (FTK_BUILD_EXAMPLES "Build the example programs." OFF)

# If FTK_BUILD_TESTS is on, enable testing.
if (FTK_BUILD_TESTS)
  enable_testing ()
  add_subdirectory (tests)
endif ()

# If FTK_BUILD_EXAMPLES is on, add the examples subdirectory to the build.
if (FTK_BUILD_EXAMPLES)
  add_subdirectory (examples)
endif ()
