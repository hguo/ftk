if (FTK_BUILD_PYFTK)
  include_directories (${CMAKE_CURRENT_SOURCE_DIR}/python/pybind11/include
    ${PYTHON_INCLUDE_DIRS})
endif ()

add_executable (track_critical_points track_critical_points.cpp)
target_link_libraries (track_critical_points PRIVATE ftk)

add_executable (track_levelsets track_levelsets.cpp)
target_link_libraries (track_levelsets PRIVATE ftk)

add_executable (track_contours track_contours.cpp)
target_link_libraries (track_contours PRIVATE ftk)

add_executable (track_tdgl_vortices track_tdgl_vortices.cpp)
target_link_libraries (track_tdgl_vortices PRIVATE ftk)

install (TARGETS track_critical_points track_levelsets track_contours track_tdgl_vortices
  EXPORT FTKTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
