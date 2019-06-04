file (GLOB_RECURSE ftk_sources *.cpp)

add_library (ftk SHARED ${ftk_sources})
# target_link_libraries (ftk)
