unset (EXECUTABLE_OUTPUT_PATH)

if (NOT "${FTK_XGC_TEST_DATA_PATH}" STREQUAL "") 
  set (FTK_TEST_XGC TRUE)
endif ()

if (FTK_HAVE_CUDA AND DETECT_GPU_COUNT)
  set (FTK_TEST_CUDA TRUE)
endif ()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/main.hh.in
  ${CMAKE_CURRENT_BINARY_DIR}/main.hh)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

if (FTK_BUILD_PYFTK)
  include_directories (${CMAKE_CURRENT_SOURCE_DIR}/python/pybind11/include
    ${PYTHON_INCLUDE_DIRS})
endif ()

add_custom_target(copy_test_data ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/data/* ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable (test_quadratic_interpolation test_quadratic_interpolation.cpp)
target_link_libraries (test_quadratic_interpolation libftk)
catch_discover_tests (test_quadratic_interpolation)

add_executable (test_regular_simplex_mesh test_regular_simplex_mesh.cpp)
target_link_libraries (test_regular_simplex_mesh libftk)
# catch_discover_tests (test_regular_simplex_mesh)

add_executable (test_array_group_stream test_array_group_stream.cpp)
target_link_libraries (test_array_group_stream libftk)

add_executable (test_array_stream test_array_stream.cpp)
target_link_libraries (test_array_stream libftk)

add_executable (test_union_find test_union_find.cpp)
target_link_libraries (test_union_find libftk)
catch_discover_tests (test_union_find)

add_executable (test_parallel_vectors test_parallel_vectors.cpp)
target_link_libraries (test_parallel_vectors libftk)
catch_discover_tests (test_parallel_vectors)

add_executable (test_polynomial test_polynomial.cpp)
target_link_libraries (test_polynomial libftk)
catch_discover_tests (test_polynomial)

add_executable (test_hoshen_kopelman test_hoshen_kopelman.cpp)
target_link_libraries (test_hoshen_kopelman libftk)
catch_discover_tests (test_hoshen_kopelman)

add_executable (test_conv test_conv.cpp)
target_link_libraries (test_conv libftk)
catch_discover_tests (test_conv)

add_executable (test_inequality_solvers test_inequality_solvers.cpp)
target_link_libraries (test_inequality_solvers libftk)
catch_discover_tests (test_inequality_solvers)

add_executable (test_matrix test_matrix.cpp)
target_link_libraries (test_matrix libftk)
catch_discover_tests (test_matrix)

add_executable (test_inverse_interpolation test_inverse_interpolation.cpp)
target_link_libraries (test_inverse_interpolation libftk)
catch_discover_tests (test_inverse_interpolation)

add_executable (test_io test_io)
target_link_libraries (test_io libftk)
catch_discover_tests (test_io)

add_executable (test_mesh test_mesh.cpp)
target_link_libraries (test_mesh libftk)
catch_discover_tests (test_mesh)

if (FTK_TEST_XGC)
  message("Adding XGC-specific tests, data_path = ${FTK_XGC_TEST_DATA_PATH}")

  # generate synthetic data
  add_custom_target (synthesize_xgc_data ALL
    # BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/xgc.synthetic.0000.h5
    COMMAND ${CMAKE_COMMAND} -E env FTK_XGC_TEST_DATA_PATH=${FTK_XGC_TEST_DATA_PATH}
    # ${PYTHON_EXECUTABLE}
    # ${Python_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/synthesize_xgc_data.py)

  add_executable (test_xgc_mesh test_xgc_mesh.cpp)
  target_link_libraries (test_xgc_mesh libftk)
  catch_discover_tests (test_xgc_mesh)
  # PROPERTIES ENVIRONMENT FTK_XGC_TEST_DATA_PATH=${FTK_XGC_TEST_DATA_PATH})
  
  add_executable (test_xgc_filament_tracking test_xgc_filament_tracking.cpp)
  target_link_libraries (test_xgc_filament_tracking libftk)
  catch_discover_tests (test_xgc_filament_tracking)
endif ()

#add_executable (test_critical_point_tracking_woven test_critical_point_tracking_woven)
#target_link_libraries (test_critical_point_tracking_woven libftk)

add_executable (test_critical_point_tracking_woven_unstructured test_critical_point_tracking_woven_unstructured.cpp)
target_link_libraries (test_critical_point_tracking_woven_unstructured libftk)

add_executable (test_critical_point_tracking_moving_extremum_2d_unstructured test_critical_point_tracking_moving_extremum_2d_unstructured.cpp)
target_link_libraries (test_critical_point_tracking_moving_extremum_2d_unstructured libftk)

add_executable (test_critical_point_tracking_moving_extremum_3d_unstructured test_critical_point_tracking_moving_extremum_3d_unstructured.cpp)
target_link_libraries (test_critical_point_tracking_moving_extremum_3d_unstructured libftk)

add_executable (test_critical_point_tracking_moving_extremum_2d test_critical_point_tracking_moving_extremum_2d)
target_link_libraries (test_critical_point_tracking_moving_extremum_2d libftk)

add_executable (test_critical_point_tracking_moving_extremum_3d test_critical_point_tracking_moving_extremum_3d)
target_link_libraries (test_critical_point_tracking_moving_extremum_3d libftk)

add_executable(test_adios2 test_adios2.cpp)
target_link_libraries (test_adios2 libftk)

add_executable (test_critical_point_tracking_merger_2d test_critical_point_tracking_merger_2d)
target_link_libraries (test_critical_point_tracking_merger_2d libftk)
catch_discover_tests (test_critical_point_tracking_merger_2d)

add_executable (test_critical_point_tracking_double_gyre test_critical_point_tracking_double_gyre)
target_link_libraries (test_critical_point_tracking_double_gyre libftk)
catch_discover_tests (test_critical_point_tracking_double_gyre)

# in situ adios2 test
if (FTK_HAVE_ADIOS2)
  add_executable (heat2d
    heat2d/HeatTransfer.cpp
    heat2d/IO_adios2.cpp
    heat2d/Settings.cpp
    heat2d/heatSimulation.cpp)
  target_link_libraries (heat2d libftk)
  add_test (NAME test_adios2_heat2d_insitu
    COMMAND ${MPIEXEC} -n 12 $<TARGET_FILE:heat2d> sim.bp 4 3 5 10 200 1 : -n 2 $<TARGET_FILE:ftk> --adios-config heat2d_adios2.xml --adios-name SimulationOutput -f cp --input sim.bp --var T --output-type traced --output heat.txt --stream)
endif ()

function (add_mpi_test target maxp args)
  add_executable (${target} ${target}.cpp)
  target_link_libraries (${target} libftk)
  if (FTK_HAVE_MPI)
    foreach (p RANGE 1 ${maxp})
      add_test (NAME ${target}_p${p}
        COMMAND ${MPIEXEC} -np ${p} $<TARGET_FILE:${target}> ${args})
    endforeach (p)
  else () # serial only
    catch_discover_tests (${target})
  endif ()
endfunction ()

add_mpi_test (test_critical_point_tracking_woven 4 "")
add_mpi_test (test_critical_point_tracking_double_gyre_unstructured 4 "")

# cli test
if (FTK_BUILD_EXECUTABLES)
  if (FTK_HAVE_MPI)
    add_test (NAME test_critical_point_tracking_cli_woven
      COMMAND ${MPIEXEC} -np 1 $<TARGET_FILE:ftk> -f cp --synthetic woven --output woven.txt)
  else ()
    add_test (NAME test_critical_point_tracking_cli_woven
      COMMAND $<TARGET_FILE:ftk> -f cp --synthetic woven --output woven.txt)
  endif ()
endif ()

# mpi test
if (FTK_HAVE_MPI)
  set (maxp 4)
  #foreach (p RANGE 1 ${maxp})
  #  add_test (NAME test_critical_point_tracking_woven_p${p} 
  #    COMMAND ${MPIEXEC} -np ${p} $<TARGET_FILE:test_critical_point_tracking_woven>)
  #endforeach (p)
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

  catch_discover_tests (test_critical_point_tracking_woven_unstructured)
  # catch_discover_tests (test_critical_point_tracking_moving_extremum_2d_unstructured) # unmature for now
endif ()

# pyftk
if (FTK_BUILD_PYFTK)
  add_test (NAME test_critical_point_tracking_py
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_critical_point_tracking.py)
  
  add_test (NAME test_numeric_py
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_numeric.py)

  set_tests_properties (
    test_critical_point_tracking_py 
    test_numeric_py
    PROPERTIES 
    ENVIRONMENT "PYTHONPATH=${LIBRARY_OUTPUT_PATH}")
endif ()

if (NOT FTK_HAVE_MPI AND PARAVIEW_USE_PYTHON AND Python_FOUND) # testing paraview plugins
  file (GENERATE
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/test_paraview_critical_point_tracking.py
    INPUT  ${CMAKE_CURRENT_SOURCE_DIR}/test_paraview_critical_point_tracking.py.in)
  add_test (NAME test_paraview_critical_point_tracker
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/test_paraview_critical_point_tracking.py)
  
  file (GENERATE
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/test_paraview_levelset_tracking.py
    INPUT  ${CMAKE_CURRENT_SOURCE_DIR}/test_paraview_levelset_tracking.py.in)
  add_test (NAME test_paraview_levelset_tracking
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/test_paraview_levelset_tracking.py)
  
  set_tests_properties (
    test_paraview_critical_point_tracker
    test_paraview_levelset_tracking
    PROPERTIES
    ENVIRONMENT "PYTHONPATH=${ParaView_PREFIX_PATH}/${PARAVIEW_PYTHONPATH}")
endif ()
