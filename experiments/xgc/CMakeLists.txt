include_directories (../common/trackball)

if (FTK_HAVE_QT5)
  set (sources main.cpp widget.cpp)
else ()
  set (sources main.cpp)
endif ()

add_executable (critical_point_tracking_xgc ${sources})

if (FTK_HAVE_QT5)
  qt5_use_modules (critical_point_tracking_xgc Widgets OpenGL)
  target_link_libraries (critical_point_tracking_xgc ftk trackball ${GLEW_SHARED_LIBRARY})
else ()
  target_link_libraries (critical_point_tracking_xgc ftk)
endif ()
