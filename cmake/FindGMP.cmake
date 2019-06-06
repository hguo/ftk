# Try to find the GNU Multiple Precision Arithmetic Library (GMP)
# See http://gmplib.org/

if (GMP_INCLUDE_DIRS AND GMP_LIBRARIES)
  set(GMP_FIND_QUIETLY TRUE)
endif (GMP_INCLUDE_DIRS AND GMP_LIBRARIES)

find_path(GMP_INCLUDE_DIRS
  NAMES
  gmp.h
  PATHS
  $ENV{GMP_DIR}
  ${INCLUDE_INSTALL_DIR}
)

find_library(GMP_LIBRARIES gmp PATHS $ENV{GMP_DIR} ${LIB_INSTALL_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG
                                  GMP_INCLUDE_DIRS GMP_LIBRARIES)
mark_as_advanced(GMP_INCLUDE_DIRS GMP_LIBRARIES)
