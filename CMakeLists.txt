cmake_minimum_required(VERSION 3.14)
project(biogeosim)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(lib)
add_subdirectory(src)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/tests/lib/Catch2/extras)
include(CTest)
include(Catch)
enable_testing()

add_subdirectory(tests)
