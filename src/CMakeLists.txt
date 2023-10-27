add_library(biogeosim_obj OBJECT
    dist.cpp
    util.cpp
    node.cpp
    tree.cpp
)

target_link_libraries(biogeosim_obj PUBLIC corax logger)
target_include_directories(biogeosim_obj PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(biogeosim
    main.cpp
)

target_link_libraries(biogeosim 
  biogeosim_obj 
  logger 
  CLI11 
  yaml-cpp 
  nlohmann_json::nlohmann_json)

set_target_properties(biogeosim
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
)

target_include_directories(biogeosim_obj PUBLIC "${CMAKE_SOURCE_DIR}/lib/pcg/include")
target_include_directories(biogeosim PUBLIC "${CMAKE_SOURCE_DIR}/lib/pcg/include")
