add_executable(bigrig_test
  dist.cpp
  tree.cpp
  split.cpp
  model.cpp
  period.cpp
  node.cpp
)

option(RIGOROUS_TESTS "Peform rigorous tests" OFF)

target_link_libraries(bigrig_test PRIVATE bigrig_obj Catch2 Catch2WithMain)

set_target_properties(bigrig_test
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/bin"
)

target_include_directories(bigrig_test PUBLIC "${CMAKE_SOURCE_DIR}/lib/pcg/include")

target_compile_definitions(bigrig_test PRIVATE
  "D_RIGOROUS=$<BOOL:${RIGOROUS_TESTS}>")

catch_discover_tests(bigrig_test
  EXTRA_ARGS --skip-benchmarks
)
