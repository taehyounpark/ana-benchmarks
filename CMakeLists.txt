cmake_minimum_required(VERSION 3.12)

project(ana-benchmarks)

add_subdirectory(ana/)
add_subdirectory(rootana/)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# set(CMAKE_BUILD_TYPE Debug)

add_executable(task_1 tasks/task_1/task.cxx)
target_compile_features(task_1 PRIVATE cxx_std_17)
target_include_directories(task_1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ANA_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_1 ana rootana ${ROOT_LIBRARIES})

add_executable(task_2 tasks/task_2/task.cxx)
target_compile_features(task_2 PRIVATE cxx_std_17)
target_include_directories(task_2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ANA_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_2 ana rootana ${ROOT_LIBRARIES})

add_executable(task_3 tasks/task_3/task.cxx)
target_compile_features(task_3 PRIVATE cxx_std_17)
target_include_directories(task_3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ANA_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_3 ana rootana ${ROOT_LIBRARIES})

add_executable(task_4 tasks/task_4/task.cxx)
target_compile_features(task_4 PRIVATE cxx_std_17)
target_include_directories(task_4 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ANA_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_4 ana rootana ${ROOT_LIBRARIES})

add_executable(task_5 tasks/task_5/task.cxx)
target_compile_features(task_5 PRIVATE cxx_std_17)
target_include_directories(task_5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ANA_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_5 ana rootana ${ROOT_LIBRARIES})

add_executable(task_6 tasks/task_6/task.cxx)
target_compile_features(task_6 PRIVATE cxx_std_17)
target_include_directories(task_6 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ANA_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_6 ana rootana ${ROOT_LIBRARIES})

add_executable(task_8 tasks/task_8/task.cxx)
target_compile_features(task_8 PRIVATE cxx_std_17)
target_include_directories(task_8 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ANA_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_8 ana rootana ${ROOT_LIBRARIES})

set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "source ana/setup.sh\n")
file(APPEND ${SETUP} "source rootana/setup.sh\n")
file(APPEND ${SETUP} "export PATH=$PATH:"${CMAKE_RUNTIME_OUTPUT_DIRECTORY})