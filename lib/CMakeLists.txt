add_subdirectory(coraxlib EXCLUDE_FROM_ALL)
add_subdirectory(logger EXCLUDE_FROM_ALL)
add_subdirectory(CLI11 EXCLUDE_FROM_ALL)

target_include_directories(biogeosim_obj PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/pcg/include")
target_include_directories(biogeosim_test PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/pcg/include")
target_include_directories(biogeosim PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/pcg/include")
