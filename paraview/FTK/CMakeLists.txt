paraview_add_plugin(FTK
  VERSION FTK_VERSION
  MODULES FTKFilters
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filters/vtk.module")

if (FTK_BUILD_TESTS)
  add_subdirectory (tests)
endif ()
