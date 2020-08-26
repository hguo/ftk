unset (EXECUTABLE_OUTPUT_PATH)

if (FTK_BUILD_PYFTK)
  include_directories (${CMAKE_CURRENT_SOURCE_DIR}/python/pybind11/include
    ${PYTHON_INCLUDE_DIRS})
endif ()

add_executable (test_quadratic_interpolation test_quadratic_interpolation.cpp)
target_link_libraries (test_quadratic_interpolation ftk)
catch_discover_tests (test_quadratic_interpolation)

add_executable (test_regular_simplex_mesh test_regular_simplex_mesh.cpp)
target_link_libraries (test_regular_simplex_mesh ftk)
# catch_discover_tests (test_regular_simplex_mesh)

add_executable (test_array_stream test_array_stream.cpp)
target_link_libraries (test_array_stream ftk)

add_executable (test_union_find test_union_find.cpp)
target_link_libraries (test_union_find ftk)
catch_discover_tests (test_union_find)

add_executable (test_parallel_vectors test_parallel_vectors.cpp)
target_link_libraries (test_parallel_vectors ftk)
catch_discover_tests (test_parallel_vectors)

add_executable (test_polynomial test_polynomial.cpp)
target_link_libraries (test_polynomial ftk)
catch_discover_tests (test_polynomial)

add_executable (test_hoshen_kopelman test_hoshen_kopelman.cpp)
target_link_libraries (test_hoshen_kopelman ftk)
catch_discover_tests (test_hoshen_kopelman)

add_executable (test_conv test_conv.cpp)
target_link_libraries (test_conv ftk)
catch_discover_tests (test_conv)

add_executable (test_inequality_solvers test_inequality_solvers.cpp)
target_link_libraries (test_inequality_solvers ftk)
catch_discover_tests (test_inequality_solvers)

add_executable (test_matrix test_matrix.cpp)
target_link_libraries (test_matrix ftk)
catch_discover_tests (test_matrix)

add_executable (test_inverse_interpolation test_inverse_interpolation.cpp)
target_link_libraries (test_inverse_interpolation ftk)
catch_discover_tests (test_inverse_interpolation)

add_executable (test_io test_io)
target_link_libraries (test_io ftk)
catch_discover_tests (test_io)

add_executable (test_critical_point_tracking_woven test_critical_point_tracking_woven)
target_link_libraries (test_critical_point_tracking_woven ftk)

add_executable (test_critical_point_tracking_moving_extremum_2d test_critical_point_tracking_moving_extremum_2d)
target_link_libraries (test_critical_point_tracking_moving_extremum_2d ftk)

add_executable (test_critical_point_tracking_moving_extremum_3d test_critical_point_tracking_moving_extremum_3d)
target_link_libraries (test_critical_point_tracking_moving_extremum_3d ftk)

add_executable(test_critical_point_tracking_double_gyre_unstructured test_critical_point_tracking_double_gyre_unstructured.cpp)
target_link_libraries (test_critical_point_tracking_double_gyre_unstructured ftk)

add_executable(test_adios2 test_adios2.cpp)
target_link_libraries (test_adios2 ftk)

# cli test
if (FTK_BUILD_EXECUTABLES)
  if (FTK_HAVE_MPI)
    add_test (NAME test_critical_point_tracking_cli_woven
      COMMAND ${MPIEXEC} -np 1 $<TARGET_FILE:track_critical_points> --synthetic woven --output woven.txt)
  else ()
    add_test (NAME test_critical_point_tracking_cli_woven
      COMMAND $<TARGET_FILE:track_critical_points> --synthetic woven --output woven.txt)
  endif ()
endif ()

# mpi test
if (FTK_HAVE_MPI)
  set (maxp 4)
  foreach (p RANGE 1 ${maxp})
    add_test (NAME test_critical_point_tracking_woven_p${p} 
      COMMAND ${MPIEXEC} -np ${p} $<TARGET_FILE:test_critical_point_tracking_woven>)
  endforeach (p)
  foreach (p RANGE 1 ${maxp})
    add_test (NAME test_critical_point_tracking_moving_extremum_2d_p${p} 
      COMMAND ${MPIEXEC} -np ${p} $<TARGET_FILE:test_critical_point_tracking_moving_extremum_2d>)
  endforeach (p)
  #foreach (p RANGE 1 ${maxp})
  #  add_test (NAME test_critical_point_tracking_moving_extremum_3d_p${p} 
  #    COMMAND ${MPIEXEC} -np ${p} $<TARGET_FILE:test_critical_point_tracking_moving_extremum_3d>)
  #endforeach (p)
else () # serial only
  catch_discover_tests (test_critical_point_tracking_woven)
  catch_discover_tests (test_critical_point_tracking_moving_extremum_2d)
  catch_discover_tests (test_critical_point_tracking_moving_extremum_3d)
endif ()

# pyftk
if (FTK_BUILD_PYFTK)
  add_test (NAME test_critical_point_tracking_py
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_critical_point_tracking.py)

  set_tests_properties (test_critical_point_tracking_py 
    PROPERTIES 
    ENVIRONMENT "PYTHONPATH=${LIBRARY_OUTPUT_PATH}")
endif ()
