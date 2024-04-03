cmake_minimum_required(VERSION 3.12)

project(ana-benchmarks)

add_subdirectory(queryosity/)
add_subdirectory(HepQuery/)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

# add_executable(task_1 tasks/task_1/task.cxx)
# target_compile_features(task_1 PRIVATE cxx_std_17)
# target_include_directories(task_1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
# target_link_libraries(task_1 HepQuery ${ROOT_LIBRARIES})

# add_executable(task_2 tasks/task_2/task.cxx)
# target_compile_features(task_2 PRIVATE cxx_std_17)
# target_include_directories(task_2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
# target_link_libraries(task_2 HepQuery ${ROOT_LIBRARIES})

# add_executable(task_3 tasks/task_3/task.cxx)
# target_compile_features(task_3 PRIVATE cxx_std_17)
# target_include_directories(task_3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
# target_link_libraries(task_3 HepQuery ${ROOT_LIBRARIES})

# add_executable(task_4 tasks/task_4/task.cxx)
# target_compile_features(task_4 PRIVATE cxx_std_17)
# target_include_directories(task_4 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
# target_link_libraries(task_4 HepQuery ${ROOT_LIBRARIES})

# add_executable(task_5 tasks/task_5/task.cxx)
# target_compile_features(task_5 PRIVATE cxx_std_17)
# target_include_directories(task_5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
# target_link_libraries(task_5 HepQuery ${ROOT_LIBRARIES})

add_executable(task_6 tasks/task_6/task.cxx)
target_compile_features(task_6 PRIVATE cxx_std_17)
target_include_directories(task_6 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_6 HepQuery ${ROOT_LIBRARIES})

add_executable(task_7 tasks/task_7/task.cxx)
target_compile_features(task_7 PRIVATE cxx_std_17)
target_include_directories(task_7 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_7 HepQuery ${ROOT_LIBRARIES})

add_executable(task_8 tasks/task_8/task.cxx)
target_compile_features(task_8 PRIVATE cxx_std_17)
target_include_directories(task_8 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR})
target_link_libraries(task_8 HepQuery ${ROOT_LIBRARIES})

set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "export PATH=$PATH:" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} "\n")
file(APPEND ${SETUP} "source HepQuery/setup.sh\n")