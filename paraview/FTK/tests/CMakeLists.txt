if (PARAVIEW_USE_PYTHON)
  file (GENERATE
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/test_ftkCriticalPointTracker2D.py
    INPUT  ${CMAKE_CURRENT_SOURCE_DIR}/test_ftkCriticalPointTracker2D.py)
  set (_vtk_build_TEST_FILE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  paraview_add_test_python(NO_RT DIRECT_DATA test_ftkCriticalPointTracker2D.py)
endif ()
