add_executable(biogeosim_test
  dist.cpp
)

target_link_libraries(biogeosim_test PRIVATE biogeosim_obj Catch2 Catch2WithMain)

set_target_properties(biogeosim_test
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/bin"
)

target_include_directories(biogeosim_test PUBLIC "${CMAKE_SOURCE_DIR}/lib/pcg/include")
