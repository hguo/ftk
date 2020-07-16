cmake_minimum_required (VERSION 3.10)
project (ftk)

set (FTK_MAJOR_VERSION 0)
set (FTK_MINOR_VERSION 0)
set (FTK_PATCH_VERSION 1)
set (FTK_VERSION ${FTK_MAJOR_VERSION}.${FTK_MINOR_VERSION}.${FTK_PATCH_VERSION})

include(GoogleTest)
include(GNUInstallDirs)
# include (ExternalProject)

if (FTK_BUILD_PARAVIEW)
  find_package (ParaView REQUIRED)
endif ()

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# dependencies
function (ftk_option name description default)
  set (FTK_USE_${name} ${default} CACHE STRING "${description}")
  set_property (CACHE FTK_USE_${name} PROPERTY STRINGS "TRUE;AUTO;FALSE")
endfunction ()

ftk_option (VTK "Use VTK" AUTO)
ftk_option (Kokkos "Use Kokkos" FALSE) # experimental
ftk_option (Qt5 "Use Qt5" AUTO)
ftk_option (MPI "Use MPI" AUTO)
ftk_option (NETCDF "Use NetCDF" AUTO)
ftk_option (HDF5 "Use HDF5" AUTO)
ftk_option (PNETCDF "Use parallel-netcdf" AUTO)
ftk_option (MPSolve "Use MPSolve" AUTO)
ftk_option (GSL "Use GSL (GNU Scientific Library)" AUTO)
ftk_option (TBB "Use TBB (Intel Thread Building Blocks)" FALSE) # this feature is experimental
ftk_option (PNG "Use PNG" AUTO)
ftk_option (LevelDB "Use LevelDB" AUTO)
ftk_option (RocksDB "Use RocksDB" AUTO)
ftk_option (Boost "Use Boost" AUTO)
ftk_option (OpenMP "Use OpenMP" AUTO)
ftk_option (CUDA "Use CUDA" FALSE) # experimental

find_package (Threads REQUIRED)

if (FTK_USE_VTK STREQUAL AUTO)
  find_package (VTK QUIET)
elseif (FTK_USE_VTK)
  find_package (VTK REQUIRED)
endif ()
if (VTK_FOUND)
  set (FTK_HAVE_VTK TRUE)
  if (VTK_MAJOR_VERSION STRLESS 9)
    include (${VTK_USE_FILE})
  endif ()
endif ()

if (FTK_USE_Kokkos STREQUAL AUTO)
  find_package (Kokkos QUIET)
elseif (FTK_USE_Kokkos)
  find_package (Kokkos REQUIRED)
endif ()
if (Kokkos_FOUND)
  set (FTK_HAVE_KOKKOS TRUE)
  include_directories (${Kokkos_INCLUDE_DIRS})
endif ()

if (FTK_USE_Qt5 STREQUAL AUTO)
  find_package (OpenGL QUIET)
  find_package (Qt5Widgets QUIET)
  find_package (Qt5OpenGL QUIET)
elseif (FTK_USE_Qt5)
  find_package (OpenGL REQUIRED)
  find_package (Qt5Widgets REQUIRED)
  find_package (Qt5OpenGL REQUIRED)
endif ()
if (Qt5OpenGL_FOUND)
  set (FTK_HAVE_QT5 TRUE)
  set (FTK_HAVE_QT TRUE)
  add_definitions("-Wno-deprecated-declarations")
endif ()

if (FTK_USE_MPI STREQUAL AUTO)
  find_package (MPI QUIET)
elseif (FTK_USE_MPI)
  find_package (MPI REQUIRED)
endif ()
if (MPI_FOUND)
  set (FTK_HAVE_MPI TRUE)
  set (CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  include_directories (${MPI_C_INCLUDE_DIRS})
else ()
  add_definitions (-DDIY_NO_MPI)
endif ()

if (FTK_USE_HDF5 STREQUAL AUTO)
  find_package (HDF5 QUIET)
elseif (FTK_USE_HDF5)
  find_package (HDF5 REQUIRED)
endif()
if (HDF5_FOUND)
  set (FTK_HAVE_HDF5 TRUE)
  include_directories (${HDF5_INCLUDE_DIRS})
endif ()

if (FTK_USE_NETCDF STREQUAL AUTO)
  find_package (NetCDF QUIET)
elseif (FTK_USE_NETCDF)
  find_package (NetCDF REQUIRED)
endif ()
if (NETCDF_FOUND)
  set (FTK_HAVE_NETCDF TRUE)
  include_directories (${NETCDF_INCLUDE_DIR})
endif ()

if (FTK_USE_MPSolve STREQUAL AUTO)
  find_package (MPSolve QUIET)
  find_package (GMP QUIET)
elseif (FTK_USE_MPSolve)
  find_package (MPSolve REQUIRED)
  find_package (GMP REQUIRED)
endif ()
if (MPSolve_FOUND AND GMP_FOUND)
  set (FTK_HAVE_MPSOLVE TRUE)
  add_definitions ("-Wno-register")
  include_directories (${MPSolve_INCLUDE_DIR} ${GMP_INCLUDE_DIRS})
endif ()

if (FTK_USE_GSL STREQUAL AUTO)
  find_package (GSL QUIET)
elseif (FTK_USE_GSL)
  find_package (GSL REQUIRED)
endif ()
if (GSL_FOUND)
  set (FTK_HAVE_GSL TRUE)
  include_directories (${GSL_INCLUDE_DIRS})
endif ()

if (FTK_USE_TBB STREQUAL AUTO)
  find_package (TBB QUIET)
elseif (FTK_USE_TBB)
  find_package (TBB REQUIRED)
endif()
if (TBB_FOUND)
  set (FTK_HAVE_TBB TRUE)
  include_directories (${TBB_INCLUDE_DIRS})
endif ()

if (FTK_USE_PNETCDF STREQUAL AUTO)
  find_package (PNetCDF QUIET)
elseif (FTK_USE_PNETCDF)
  find_package (PNetCDF REQUIRED)
endif()
if (PNETCDF_FOUND)
  set (FTK_HAVE_PNETCDF TRUE)
  include_directories (${PNETCDF_INCLUDE_DIR})
endif ()

if (FTK_USE_PNG STREQUAL AUTO)
  find_package (PNG QUIET)
elseif (FTK_USE_PNG)
  find_package (PNG REQUIRED)
endif ()
if (PNG_FOUND)
  set (FTK_HAVE_PNG TRUE)
  include_directories (${PNG_INCLUDE_DIRS})
endif ()

if (FTK_USE_RocksDB STREQUAL AUTO)
  find_package (RocksDB QUIET)
elseif (FTK_USE_RocksDB)
  find_package (RocksDB REQUIRED)
endif ()
if (RocksDB_FOUND)
  set (FTK_HAVE_ROCKSDB TRUE)
  include_directories (${RocksDB_INCLUDE_DIR})
endif ()

if (FTK_USE_OpenMP STREQUAL AUTO)
  find_package (OpenMP QUIET)
elseif (FTK_USE_OpenMP)
  find_package (OpenMP REQUIRED)
endif ()
if (OpenMP_CXX_FOUND)
  set (FTK_HAVE_OPENMP TRUE)
endif ()

if (FTK_USE_CUDA STREQUAL AUTO)
  find_package (CUDA QUIET)
elseif (FTK_USE_CUDA)
  find_package (CUDA REQUIRED)
endif ()
if (CUDA_FOUND)
  set (FTK_HAVE_CUDA TRUE)
  set (CUDA_PROPAGATE_HOST_FLAGS OFF)
  set (CUDA_NVCC_FLAGS "--expt-relaxed-constexpr")
endif ()

set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)
if (FTK_USE_Boost STREQUAL AUTO)
  find_package (Boost QUIET COMPONENTS system)
elseif (FTK_USE_Boost)
  find_package (Boost REQUIRED COMPONENTS system)
endif ()
if (Boost_SYSTEM_FOUND)
  set (FTK_HAVE_BOOST TRUE)
  # add_definitions (-DBOOST_NO_AUTO_PTR) # c++17 deprecated auto_ptr
  include_directories (${Boost_INCLUDE_DIRS})
endif ()

### 

set (FTK_INCLUDE_DIR 
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_BINARY_DIR}/include
)
include_directories (${FTK_INCLUDE_DIR})

option (FTK_BUILD_TESTS "Build tests" OFF)
option (FTK_BUILD_EXAMPLES "Build examples" ON)
option (FTK_BUILD_EXPERIMENTS "Build experimental code" OFF)
option (FTK_BUILD_PARAVIEW "Build ParaView plugins" OFF)
option (FTK_BUILD_PYFTK "Build Python bindings" OFF)
set (FTK_FP_PRECISION "32768" CACHE STRING "Fixed point precision")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/ftk/ftk_config.hh.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/ftk/ftk_config.hh)

install (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ftk DESTINATION include)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/include/ftk/ftk_config.hh DESTINATION include)

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

# If FTK_BUILD_TESTS is on, enable testing.
if (FTK_BUILD_TESTS)
  # find_package (GTest REQUIRED)
  enable_testing ()
  add_subdirectory (tests)
endif ()

# If FTK_BUILD_EXAMPLES is on, add the examples subdirectory to the build.
if (FTK_BUILD_EXAMPLES)
  add_subdirectory (examples)
endif ()

if (FTK_BUILD_EXPERIMENTS)
  add_subdirectory (experiments)
endif ()

if (FTK_BUILD_PARAVIEW)
  add_subdirectory (paraview)
endif ()

if (FTK_BUILD_PYFTK)
  add_subdirectory (python)
endif ()

add_subdirectory (src)
add_subdirectory (standalone)
