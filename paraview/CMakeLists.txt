set (_paraview_plugin_default_FTK ON)
set (BUILD_SHARED_LIBS ON)
# option (PARAVIEW_PLUGIN_ENABLE_FTK ON)

paraview_plugin_scan(
  PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/FTK/paraview.plugin"
  PROVIDES_PLUGINS plugins
  REQUIRES_MODULES modules)

foreach (module IN LISTS REQUIRES_MODULES)
  if (NOT TARGET "${module}")
    message (FATAL_ERROR
      "The required ${module} module is not available in the "
      "ParaView build.")
  endif ()
endforeach()

paraview_plugin_build(
  RUNTIME_DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY_SUBDIRECTORY "${PARAVIEW_PLUGIN_SUBDIR}"
  TARGET paraview_plugins
  PLUGINS ${plugins})
