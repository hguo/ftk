set (sources main.cpp)

if (FTK_HAVE_QT5)
  set (CMAKE_AUTOMOC ON)
  list (APPEND sources ../common/trackball.cpp widget.cpp)
endif ()

add_executable (ex_critical_point_tracking_2d ${sources})
target_link_libraries (ex_critical_point_tracking_2d ftk)

if (FTK_HAVE_QT5)
  qt5_use_modules (ex_critical_point_tracking_2d Widgets OpenGL)
  target_link_libraries (ex_critical_point_tracking_2d ${OPENGL_LIBRARIES})
endif ()
