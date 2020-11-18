if (FTK_BUILD_PYFTK)
  include_directories (${PYTHON_INCLUDE_DIRS})
endif ()

add_subdirectory (cli)

set (ftk_sources
  numeric/polynomial_solver.cpp
  io/tdgl/BDATReader.cpp
  io/tdgl/glpp/GL_post_process.cpp
  io/tdgl/glpp/paramfile.cpp
  io/tdgl/glpp/fileutils.cpp
  io/tdgl/glpp/stringutils.cpp
  io/tdgl/glpp/datafile.cpp
  io/tdgl/tdgl.cpp
  io/tdgl/GLGPU_IO_Helper.cpp)

if (FTK_HAVE_CUDA)
  set (ftk_cuda_sources
    filters/critical_point_tracer_3d_regular.cu
    filters/tdgl_vortex_tracker_3d_regular.cu
    filters/critical_point_tracer_2d_regular.cu
    filters/contour_tracker_3d_regular.cu)
endif ()

add_library (libftk SHARED ${ftk_sources} ${ftk_cuda_sources})

target_link_libraries (libftk ${CMAKE_THREAD_LIBS_INIT})

if (FTK_HAVE_ADIOS2)
  target_link_libraries (libftk ${ADIOS2_LIBRARIES})
endif ()

if (FTK_HAVE_BOOST)
  target_link_libraries (libftk ${Boost_LIBRARIES})
endif ()

if (FTK_HAVE_CGAL)
  target_link_libraries (libftk CGAL::CGAL)
endif ()

if (FTK_HAVE_GMP)
  target_link_libraries (libftk ${GMP_LIBRARIES})
endif ()

if (FTK_HAVE_VTK)
  target_link_libraries (libftk ${VTK_LIBRARIES})
endif ()

if (FTK_HAVE_KOKKOS)
  target_link_libraries (libftk ${Kokkos_LIBRARIES})
endif ()

if (FTK_HAVE_NETCDF)
  target_link_libraries (libftk ${NETCDF_LIBRARY})
endif ()

if (FTK_HAVE_HDF5)
  target_link_libraries (libftk ${HDF5_LIBRARIES})
endif ()

if (FTK_HAVE_MPSOLVE)
  target_link_libraries (libftk ${MPSolve_LIBRARY})
endif ()

if (FTK_HAVE_GSL)
  target_link_libraries (libftk ${GSL_LIBRARIES})
endif ()

if (FTK_HAVE_TBB)
  target_link_libraries (libftk ${TBB_LIBRARIES})
endif ()

if (FTK_HAVE_MPI)
  target_link_libraries (libftk ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif ()

if (FTK_HAVE_PNG)
  target_link_libraries (libftk ${PNG_LIBRARIES})
endif ()

if (FTK_HAVE_ROCKSDB)
  target_link_libraries (libftk ${RocksDB_LIBRARY})
endif ()

if (FTK_HAVE_OPENMP)
  target_link_libraries (libftk OpenMP::OpenMP_CXX)
endif ()

#if (FTK_HAVE_CUDA)
#  target_link_libraries (libftk ${CUDA_LIBRARIES})
#endif ()

#if (FTK_HAVE_PYBIND11)
#  target_link_libraries (libftk ${PYTHON_LIBRARY})
#endif ()

install (TARGETS libftk
  EXPORT FTKTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
