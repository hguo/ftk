cmake_minimum_required (VERSION 3.1.3)
project (ftk)

set (FTK_MAJOR_VERSION 0)
set (FTK_MINOR_VERSION 0)
set (FTK_PATCH_VERSION 0)
set (FTK_VERSION ${FTK_MAJOR_VERSION}.${FTK_MINOR_VERSION}.${FTK_PATCH_VERSION})

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)

# include(CMakeDependentOption)
# cmake_dependent_option (FTK_BUILD_THREADS "Build FTK with threading" ON "NOT WIN32" OFF)
# option (FTK_BUILD_MPI "Build FTK with MPI" OFF)
option (FTK_BUILD_TESTS "Build the test programs." OFF)
option (FTK_BUILD_EXAMPLES "Build the example programs." OFF)

install (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ftk DESTINATION include)
install (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/hypermesh DESTINATION include)

set (FTK_INCLUDE_DIR 
  ${CMAKE_CURRENT_SOURCE_DIR}/include)
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

# Threads
if                          (NOT FTK_BUILD_THREADS)
    add_definitions         (-DDIY_NO_THREADS)
else                        (NOT FTK_BUILD_THREADS)
    find_package            (Threads)
    set                     (libraries ${libraries}     ${CMAKE_THREAD_LIBS_INIT})
endif                       (NOT FTK_BUILD_THREADS)

# MPI
if                          (NOT FTK_BUILD_MPI)
    add_definitions         (-DDIY_NO_MPI)
else                        (NOT FTK_BUILD_MPI)
    find_package            (MPI REQUIRED)
    set                     (mpi_libraries              ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
    set                     (libraries ${libraries}     ${mpi_libraries})
endif                       (NOT FTK_BUILD_MPI)

if                          (FTK_BUILD_MPI)
    include_directories     (SYSTEM ${MPI_INCLUDE_PATH})
endif                       (FTK_BUILD_MPI)

if (FTK_BUILD_TESTS OR FTK_BUILD_EXAMPLES)
  find_package (Threads REQUIRED)
  set (CMAKE_THREAD_PREFER_PTHREAD)

  find_package (MPI QUIET)
  if (MPI_FOUND)
  else ()
    add_definitions (-DDIY_NO_MPI)
  endif ()

  find_package (VTK QUIET)
  if (VTK_FOUND)
    include (${VTK_USE_FILE})
  endif ()

  find_package (OpenGL QUIET)
  find_package (Qt5Widgets QUIET)
  find_package (Qt5OpenGL QUIET)

  find_package (NetCDF QUIET)
  if (NETCDF_FOUND)
    include_directories (${NETCDF_INCLUDE_DIR})
  endif ()
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
