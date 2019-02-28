cmake_minimum_required (VERSION 3.1.3)
project (ftk)

set (FTK_MAJOR_VERSION 0)
set (FTK_MINOR_VERSION 0)
set (FTK_PATCH_VERSION 0)
set (FTK_VERSION ${FTK_MAJOR_VERSION}.${FTK_MINOR_VERSION}.${FTK_PATCH_VERSION})

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (Git QUIET)
if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  message (STATUS "Git submodule update")
  execute_process (COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE GIT_SUBMODULE_RESULT)
  if (NOT GIT_SUBMODULE_RESULT EQUAL "0")
    message (FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMODULE_RESULT}, please checkout submodules")
  endif()
endif()

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/hypermesh/include)

option (FTK_BUILD_TESTS "Build the test programs." OFF)
option (FTK_BUILD_EXAMPLES "Build the example programs." OFF)

install (DIRECTORY 
  ${CMAKE_CURRENT_SOURCE_DIR}/include/ftk 
  DESTINATION include)

install (DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/hypermesh/include/hypermesh
  DESTINATION include)

set (FTK_INCLUDE_DIR 
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/hypermesh/include)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FTKConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/FTKConfig.cmake)

include (CMakePackageConfigHelpers)
set (FTK_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
set (FTK_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FTKConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/install/FTKConfig.cmake
  PATH_VARS FTK_INCLUDE_DIR FTK_ROOT_DIR
  INSTALL_DESTINATION lib/cmake
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/install/FTKConfig.cmake
  DESTINATION lib/cmake)

if (FTK_BUILD_TESTS OR FTK_BUILD_EXAMPLES)
  find_package (VTK QUIET)
  if (VTK_FOUND)
    include (${VTK_USE_FILE})
  endif ()
  find_package (NetCDF QUIET)
endif ()

# If FTK_BUILD_TESTS is on, enable testing.
if (FTK_BUILD_TESTS)
  find_package (GTest REQUIRED)
  enable_testing ()
  add_subdirectory (tests)
endif ()

# If FTK_BUILD_EXAMPLES is on, add the examples subdirectory to the build.
if (FTK_BUILD_EXAMPLES)
  add_subdirectory (examples)
endif ()

