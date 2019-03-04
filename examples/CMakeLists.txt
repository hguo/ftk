add_subdirectory (tracking_graph)
add_subdirectory (critical_point_tracking_2d)

if (NETCDF_FOUND)
  add_subdirectory (critical_point_tracking_3d)
endif ()
