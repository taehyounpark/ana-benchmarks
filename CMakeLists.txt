cmake_minimum_required(VERSION 3.12)

project(ana-benchmarks)

add_subdirectory(ana/)
add_subdirectory(rootana/)

set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "source ana/setup.sh\n")
file(APPEND ${SETUP} "source rootana/setup.sh\n")