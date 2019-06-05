add_executable (ex_critical_point_extraction_3d ex_critical_point_extraction_3d.cpp)
target_link_libraries (ex_critical_point_extraction_3d ${CMAKE_THREAD_LIBS_INIT})

if (NETCDF_FOUND)
  target_compile_definitions (ex_critical_point_extraction_3d PUBLIC HAVE_NETCDF)
  target_link_libraries (ex_critical_point_extraction_3d ${NETCDF_LIBRARIES})
endif ()
