include_directories (../common/trackball)

if (FTK_HAVE_QT5)
  set (sources main.cpp widget.cpp)
else ()
  set (sources main.cpp)
endif ()

add_executable (track_critical_points_xgc ${sources})

if (FTK_HAVE_QT5)
  qt5_use_modules (track_critical_points_xgc Widgets OpenGL)
  target_link_libraries (track_critical_points_xgc ftk trackball ${GLEW_SHARED_LIBRARY})
else ()
  target_link_libraries (track_critical_points_xgc ftk)
endif ()
