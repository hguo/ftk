add_executable (example_critical_point_tracking_2d main.cpp)

if (VTK_FOUND)
  message ("VTK Found.")
  target_compile_definitions (example_critical_point_tracking_2d PUBLIC HAVE_VTK)
  target_link_libraries (example_critical_point_tracking_2d ${VTK_LIBRARIES})
endif ()
