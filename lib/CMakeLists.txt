include(FetchContent)

FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/nlohmann/json
  GIT_SHALLOW TRUE
  GIT_TAG v3.12.0
  EXCLUDE_FROM_ALL
)

FetchContent_Declare(
  csv  
  GIT_REPOSITORY https://github.com/vincentlaucsb/csv-parser.git
  GIT_SHALLOW TRUE
  GIT_TAG 2.3.0
  EXCLUDE_FROM_ALL
)

FetchContent_MakeAvailable(json csv)

add_subdirectory(coraxlib EXCLUDE_FROM_ALL)
add_subdirectory(logger EXCLUDE_FROM_ALL)
add_subdirectory(CLI11 EXCLUDE_FROM_ALL)
add_subdirectory(yaml-cpp EXCLUDE_FROM_ALL)
