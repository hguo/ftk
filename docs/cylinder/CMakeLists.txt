cmake_minimum_required(VERSION 3.17)
project (ftk-cylinder)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (FTK REQUIRED)

add_executable (ftk-cylinder cylinder.cpp)
target_link_libraries (ftk-cylinder FTK::libftk)
