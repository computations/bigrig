add_executable(biogeosim
    main.cpp
    cli.cpp
)

target_link_libraries(biogeosim corax logger)

set_target_properties(biogeosim
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
)
