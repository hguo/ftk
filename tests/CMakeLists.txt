include_directories (${GTEST_INCLUDE_DIR})

add_executable (test_matrix test_matrix.cpp)
target_link_libraries (test_matrix ${GTEST_BOTH_LIBRARIES})

add_executable (test_polynomial test_polynomial.cpp)
target_link_libraries (test_polynomial ${GTEST_BOTH_LIBRARIES})

add_executable (test_inverse_interpolation test_inverse_interpolation.cpp)
target_link_libraries (test_inverse_interpolation ${GTEST_BOTH_LIBRARIES})

add_executable (test_quadratic_interpolation test_quadratic_interpolation.cpp)
target_link_libraries (test_quadratic_interpolation ${GTEST_BOTH_LIBRARIES})

gtest_discover_tests (test_matrix)
gtest_discover_tests (test_polynomial)
gtest_discover_tests (test_inverse_interpolation)
gtest_discover_tests (test_quadratic_interpolation)
