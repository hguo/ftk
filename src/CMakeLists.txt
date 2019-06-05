file (GLOB_RECURSE ftk_sources *.cpp)

add_library (ftk SHARED ${ftk_sources})

target_link_libraries (ftk ${CMAKE_THREAD_LIBS_INIT})

if (FTK_HAVE_VTK)
  target_link_libraries (ftk ${VTK_LIBRARIES})
endif ()

if (FTK_HAVE_NETCDF)
  target_link_libraries (ftk ${NETCDF_LIBRARIES})
endif ()

if (FTK_HAVE_MPSOLVE)
  target_link_libraries (ftk ${MPSolve_LIBRARY} ${GMP_LIBRARIES})
endif ()
