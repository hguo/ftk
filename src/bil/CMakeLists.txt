set (BIL_SOURCES 
  bil.c
  bil_block.c
  bil_misc.c
  bil_pio.c
  bil_sched.c
  bil_timing.c
)

include_directories (${CMAKE_SOURCE_DIR}/include/ftk/external/bil)

if (NOT FTK_USE_PNETCDF)
  add_definitions (-DDISABLE_PNETCDF)
endif ()

add_library (bil ${BIL_SOURCES})
target_link_libraries (bil ${MPI_LIBRARIES}) # ${PNETCDF_LIBRARY})

install (TARGETS bil
  EXPORT FTKTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
