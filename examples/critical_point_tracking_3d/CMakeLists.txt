add_executable (ex_critical_point_tracking_3d ex_critical_point_tracking_3d.cpp)
target_link_libraries (ex_critical_point_tracking_3d ${CMAKE_THREAD_LIBS_INIT})

if (VTK_FOUND)
  target_compile_definitions (ex_critical_point_tracking_3d PUBLIC HAVE_VTK)
  target_link_libraries (ex_critical_point_tracking_3d ${VTK_LIBRARIES})
endif ()

if (NETCDF_FOUND)
  target_compile_definitions (ex_critical_point_tracking_3d PUBLIC HAVE_NETCDF)
  target_link_libraries (ex_critical_point_tracking_3d ${NETCDF_LIBRARIES})
endif ()
