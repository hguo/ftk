add_executable              (shared_memory_parallel     shared_memory_parallel.cpp)
target_link_libraries       (shared_memory_parallel     ${libraries})

add_executable              (distributed_memory_parallel     distributed_memory_parallel.cpp)
target_link_libraries       (distributed_memory_parallel     ${libraries})