add_executable (track_critical_points track_critical_points.cpp)
target_link_libraries (track_critical_points ftk)

add_executable (track_levelsets track_levelsets.cpp)
target_link_libraries (track_levelsets ftk)

if (FTK_BUILD_PYFTK)
  add_subdirectory (pybind11)
  pybind11_add_module (pyftk pyftk.cpp)
endif ()
