INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(DECAF_CXX_INCLUDE_DIR decaf.hpp HINTS
  ${DECAF_PREFIX}/include/decaf
  /usr/include/decaf
  /usr/local/include/decaf
  /opt/local/include/decaf
  /sw/include/decaf
)
FIND_LIBRARY(DECAF_CXX_DATA_MODEL_LIBRARY NAMES bredala_datamodel HINTS
  ${DECAF_PREFIX}/lib
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  /sw/lib
)

SET(DECAF_CXX_TRANSPORT_LIBRARY "")

FIND_LIBRARY(DECAF_CXX_TRANSPORT_MPI_LIBRARY NAMES bredala_transport_mpi HINTS
  ${DECAF_PREFIX}/lib
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  /sw/lib
)

IF(DECAF_CXX_TRANSPORT_MPI_LIBRARY)
  LIST(APPEND DECAF_CXX_TRANSPORT_LIBRARY ${DECAF_CXX_TRANSPORT_MPI_LIBRARY})
ENDIF(DECAF_CXX_TRANSPORT_MPI_LIBRARY)

FIND_LIBRARY(DECAF_CXX_TRANSPORT_CCI_LIBRARY NAMES bredala_transport_cci HINTS
  ${DECAF_PREFIX}/lib
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  /sw/lib
)

FIND_LIBRARY(DECAF_CXX_RUNTIME_LIBRARY NAMES decaf HINTS
    ${DECAF_PREFIX}/lib
    /usr/lib
    /usr/local/lib
    /opt/local/lib
    /sw/lib
)

IF(DECAF_CXX_TRANSPORT_CCI_LIBRARY)
  LIST(APPEND DECAF_CXX_TRANSPORT_LIBRARY ${DECAF_CXX_TRANSPORT_CCI_LIBRARY})
ENDIF(DECAF_CXX_TRANSPORT_CCI_LIBRARY)



FIND_PATH(DECAF_C_INCLUDE_DIR bredala.h HINTS
  ${DECAF_PREFIX}/include/decaf/C
  /usr/include/decaf/C
  /usr/local/include/decaf/C
  /opt/local/include/decaf/C
  /sw/include/decaf/C
)
FIND_LIBRARY(DECAF_C_DATA_MODEL_LIBRARY NAMES bca HINTS
  ${DECAF_PREFIX}/lib
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  /sw/lib
)

FIND_LIBRARY(DECAF_C_RUNTIME_LIBRARY NAMES dca HINTS
  ${DECAF_PREFIX}/lib
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  /sw/lib
)

FIND_LIBRARY(MANALA_LIBRARY NAMES manala HINTS
    ${DECAF_PREFIX}/lib
    /usr/lib
    /usr/local/lib
    /opt/local/lib
    /sw/lib
)

find_package_handle_standard_args(Decaf DEFAULT_MSG
    DECAF_CXX_INCLUDE_DIR
    DECAF_CXX_DATA_MODEL_LIBRARY
    DECAF_CXX_TRANSPORT_LIBRARY
    DECAF_C_INCLUDE_DIR
    DECAF_C_DATA_MODEL_LIBRARY
    DECAF_C_RUNTIME_LIBRARY
    MANALA_LIBRARY)

if(DECAF_CXX_INCLUDE_DIR)
    SET(DECAF_FOUND 1 CACHE BOOL "Found decaf libraries")
    STRING(REGEX REPLACE "include/decaf" "include" DECAF_CXX_INCLUDE_DIR ${DECAF_CXX_INCLUDE_DIR})
    STRING(REGEX REPLACE "include/decaf/C" "include" DECAF_C_INCLUDE_DIR ${DECAF_C_INCLUDE_DIR})
endif(DECAF_CXX_INCLUDE_DIR)


MARK_AS_ADVANCED(
  DECAF_CXX_INCLUDE_DIR
  DECAF_CXX_DATA_MODEL_LIBRARY
  DECAF_CXX_DATA_TRANSPORT_LIBRARY
  DECAF_CXX_RUNTIME_LIBRARY
  DECAF_C_INCLUDE_DIR 
  DECAF_C_DATA_MODEL_LIBRARY 
  DECAF_C_RUNTIME_LIBRARY
  DECAF_FOUND
)


