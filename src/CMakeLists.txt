add_executable(biogeosim
    main.cpp
    dist.cpp
    util.cpp
    node.cpp
    tree.cpp
)

target_link_libraries(biogeosim corax logger CLI11)

set_target_properties(biogeosim
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
)
