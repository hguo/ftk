add_subdirectory (pybind11)
pybind11_add_module (pyftk pyftk.cpp)

# target_link_libraries (pyftk PRIVATE ftk) # no need to link
