set (sources main.cpp)

if (OpenGL_FOUND AND Qt5OpenGL_FOUND)
  set (CMAKE_AUTOMOC ON)
  list (APPEND sources ../common/trackball.cpp widget.cpp)
endif ()

add_executable (ex_critical_point_tracking_2d ${sources})
target_link_libraries (ex_critical_point_tracking_2d ${CMAKE_THREAD_LIBS_INIT})

if (VTK_FOUND)
  target_compile_definitions (ex_critical_point_tracking_2d PUBLIC HAVE_VTK)
  target_link_libraries (ex_critical_point_tracking_2d ${VTK_LIBRARIES})
endif ()
  
if (OpenGL_FOUND AND Qt5OpenGL_FOUND)
  target_compile_definitions (ex_critical_point_tracking_2d PUBLIC HAVE_QT)
  qt5_use_modules (ex_critical_point_tracking_2d Widgets OpenGL)
  target_link_libraries (ex_critical_point_tracking_2d ${OPENGL_LIBRARIES})
endif ()
