# Find FTK, the feature tracking kit: https://github.com/hguo/ftk

find_path(FTK_INCLUDE_DIR NAMES ftk/tracking_graph.hh)

if (FTK_INCLUDE_DIR)
  set(FTK_FOUND TRUE)
endif ()

if (FTK_FOUND)
  if(NOT FTK_FIND_QUIETLY)
    message(STATUS "Found FTK: ${FTK_LIBRARY}")
  endif(NOT FTK_FIND_QUIETLY)
else(FTK_FOUND)
  if(FTK_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find FTK library.  Please check out https://github.com/hguo/ftk")
  endif(FTK_FIND_REQUIRED)
endif(FTK_FOUND)
