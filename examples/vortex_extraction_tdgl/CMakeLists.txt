set (sources 
  main.cpp
  io/BDATReader.cpp
  io/GLGPU_IO_Helper.cpp
)

if (FTK_HAVE_QT5)
  set (CMAKE_AUTOMOC ON)
  list (APPEND sources ../common/trackball.cpp widget.cpp)
endif ()

add_executable (ex_vortex_extraction_tdgl ${sources})
target_link_libraries (ex_vortex_extraction_tdgl ftk)

if (FTK_HAVE_QT5)
  qt5_use_modules (ex_vortex_extraction_tdgl Widgets OpenGL)
  target_link_libraries (ex_vortex_extraction_tdgl ${OPENGL_LIBRARIES})
endif ()
