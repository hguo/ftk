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

add_executable (test_critical_point_tracking test_critical_point_tracking)
target_link_libraries (test_critical_point_tracking ftk)
catch_discover_tests (test_critical_point_tracking)
