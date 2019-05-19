if (OpenGL_FOUND AND Qt5OpenGL_FOUND AND NETCDF_FOUND)
  set (CMAKE_AUTOMOC ON)
  
  set (sources main.cpp widget.cpp ../common/trackball.cpp)
  add_executable (bout_vis ${sources})
  qt5_use_modules (bout_vis Widgets OpenGL)
  
  target_compile_definitions (bout_vis PUBLIC HAVE_QT HAVE_NETCDF)
  target_link_libraries (bout_vis ${CMAKE_THREAD_LIBS_INIT} ${OPENGL_LIBRARIES} ${NETCDF_LIBRARIES})
endif ()
