add_executable              (shared_memory_parallel     shared_memory_parallel.cpp)
target_link_libraries       (shared_memory_parallel     ${libraries})

add_executable              (distributed_memory_parallel_exchange     distributed_memory_parallel_exchange.cpp)
target_link_libraries       (distributed_memory_parallel_exchange     ${libraries})

add_executable              (distributed_memory_parallel_iexchange     distributed_memory_parallel_iexchange.cpp)
target_link_libraries       (distributed_memory_parallel_iexchange     ${libraries})