add_subdirectory (gtest)

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/gtest/include)
set (GTEST_BOTH_LIBRARIES gtest gtest_main)

add_executable (test_matrix test_matrix.cpp)
target_link_libraries (test_matrix ftk ${GTEST_BOTH_LIBRARIES})

add_executable (test_conv test_conv.cpp)
target_link_libraries (test_conv ftk ${GTEST_BOTH_LIBRARIES})

add_executable (test_inequality_solvers test_inequality_solvers.cpp)
target_link_libraries (test_inequality_solvers ftk ${GTEST_BOTH_LIBRARIES})

add_executable (test_polynomial test_polynomial.cpp)
target_link_libraries (test_polynomial ftk ${GTEST_BOTH_LIBRARIES})

add_executable (test_inverse_interpolation test_inverse_interpolation.cpp)
target_link_libraries (test_inverse_interpolation ftk ${GTEST_BOTH_LIBRARIES})

add_executable (test_parallel_vectors test_parallel_vectors.cpp)
target_link_libraries (test_parallel_vectors ftk ${GTEST_BOTH_LIBRARIES})

add_executable (test_quadratic_interpolation test_quadratic_interpolation.cpp)
target_link_libraries (test_quadratic_interpolation ftk ${GTEST_BOTH_LIBRARIES})

add_executable (test_union_find test_union_find.cpp)
target_link_libraries (test_union_find ftk ${GTEST_BOTH_LIBRARIES})

add_executable (test_regular_simplex_mesh test_regular_simplex_mesh.cpp)
target_link_libraries (test_regular_simplex_mesh ftk ${GTEST_BOTH_LIBRARIES})

gtest_discover_tests (test_matrix)
gtest_discover_tests (test_conv)
gtest_discover_tests (test_polynomial)
gtest_discover_tests (test_inequality_solvers)
gtest_discover_tests (test_inverse_interpolation)
gtest_discover_tests (test_parallel_vectors)
gtest_discover_tests (test_quadratic_interpolation)
gtest_discover_tests (test_union_find)
