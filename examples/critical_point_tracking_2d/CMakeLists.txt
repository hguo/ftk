add_executable (ex_critical_point_tracking_2d ex_critical_point_tracking_2d.cpp)

if (VTK_FOUND)
  target_compile_definitions (ex_critical_point_tracking_2d PUBLIC HAVE_VTK)
  target_link_libraries (ex_critical_point_tracking_2d ${VTK_LIBRARIES})
endif ()
