file (GLOB_RECURSE ftk_sources *.cpp)

add_library (ftk SHARED ${ftk_sources})

target_link_libraries (ftk ${CMAKE_THREAD_LIBS_INIT})

if (FTK_HAVE_BOOST)
  target_link_libraries (ftk ${Boost_LIBRARIES})
endif ()

if (FTK_HAVE_VTK)
  target_link_libraries (ftk ${VTK_LIBRARIES})
endif ()

if (FTK_HAVE_KOKKOS)
  target_link_libraries (ftk ${Kokkos_LIBRARIES})
endif ()

if (FTK_HAVE_NETCDF)
  target_link_libraries (ftk ${NETCDF_LIBRARY})
endif ()

if (FTK_HAVE_HDF5)
  target_link_libraries (ftk ${HDF5_LIBRARIES})
endif ()

if (FTK_HAVE_MPSOLVE)
  target_link_libraries (ftk ${MPSolve_LIBRARY} ${GMP_LIBRARIES})
endif ()

if (FTK_HAVE_GSL)
  target_link_libraries (ftk ${GSL_LIBRARIES})
endif ()

if (FTK_HAVE_TBB)
  target_link_libraries (ftk ${TBB_LIBRARIES})
endif ()

if (FTK_HAVE_MPI)
  target_link_libraries (ftk ${MPI_C_LIBRARIES})
endif ()

if (FTK_HAVE_PNG)
  target_link_libraries (ftk ${PNG_LIBRARIES})
endif ()

if (FTK_HAVE_ROCKSDB)
  target_link_libraries (ftk ${RocksDB_LIBRARY})
endif ()

if (FTK_HAVE_OPENMP)
  target_link_libraries (ftk OpenMP::OpenMP_CXX)
endif ()
