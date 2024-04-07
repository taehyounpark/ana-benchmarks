cmake_minimum_required(VERSION 3.12)

project(ana-benchmarks)

add_subdirectory(queryosity/)
add_subdirectory(HepQuery/)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

function(add_benchmark TASK_NAME)
    add_executable(${TASK_NAME} tasks/${TASK_NAME}/task.cxx)
    target_compile_features(${TASK_NAME} PRIVATE cxx_std_17)
    target_include_directories(${TASK_NAME} PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${ROOT_INCLUDE_DIR}
        ${ROOT_LIBRARY_DIR}
    )
    target_link_libraries(${TASK_NAME} HepQuery ${ROOT_LIBRARIES})
endfunction(add_benchmark)

add_benchmark(task_1)
add_benchmark(task_2)
add_benchmark(task_3)
add_benchmark(task_4)
add_benchmark(task_6)
add_benchmark(task_7)
add_benchmark(task_8)

set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "export PATH=$PATH:" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} "\n")
file(APPEND ${SETUP} "source HepQuery/setup.sh\n")