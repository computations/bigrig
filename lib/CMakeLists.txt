include(FetchContent)

FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/nlohmann/json
  GIT_TAG v3.12.0)

FetchContent_GetProperties(json)
if(NOT json_POPULATED)
  FetchContent_MakeAvailable(json)
endif()

add_subdirectory(coraxlib EXCLUDE_FROM_ALL)
add_subdirectory(logger EXCLUDE_FROM_ALL)
add_subdirectory(CLI11 EXCLUDE_FROM_ALL)
add_subdirectory(yaml-cpp EXCLUDE_FROM_ALL)
