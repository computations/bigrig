cmake_minimum_required(VERSION 3.14)
project(bigrig)

add_subdirectory(lib)
add_subdirectory(src)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/tests/lib/Catch2/extras)
include(CTest)
include(Catch)
enable_testing()

add_subdirectory(tests)
