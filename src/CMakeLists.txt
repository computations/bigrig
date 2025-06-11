add_library(bigrig_obj OBJECT
    dist.cpp
    model.cpp
    node.cpp
    tree.cpp
    split.cpp
)

add_library(bigrig_interface_obj OBJECT
  clioptions.cpp
  io.cpp
)

find_package(OpenMP)

target_compile_features(bigrig_obj PUBLIC cxx_std_23)
target_compile_features(bigrig_interface_obj PUBLIC cxx_std_23)

target_compile_options(bigrig_obj PUBLIC -fdiagnostics-color=always)
target_compile_options(bigrig_interface_obj PUBLIC -fdiagnostics-color=always)

if(OpenMP_FOUND) 
  target_link_libraries(bigrig_obj PUBLIC OpenMP::OpenMP_CXX)
endif()

option(ENABLE_BMI "Enable Bit Manipulation Instructions (BMI) instructions" OFF)
if(ENABLE_BMI)
  target_compile_options(bigrig_obj PUBLIC -mbmi)
endif()

option(ENABLE_BMI2 "Enable Bit Manipulation Instructions 2 (BMI2) instructions"
  OFF)
if(ENABLE_BMI2)
  target_compile_options(bigrig_obj PUBLIC -mbmi2)
endif()

option(ENABLE_ABM "Enable Advanced Bit Manipulation (ABM) instructions" OFF)
if(ENABLE_ABM)
  target_compile_options(bigrig_obj PUBLIC -mabm)
endif()

option(ENABLE_NATIVE_FLAGS "Enable Bit Manipulation Instructions (BMI) instructions" ON)
if(ENABLE_NATIVE_FLAGS)
  target_compile_options(bigrig_obj PUBLIC -march=native)
endif()

target_link_libraries(bigrig_obj PUBLIC corax logger)
target_include_directories(bigrig_obj PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(bigrig_interface_obj PUBLIC logger CLI11 yaml-cpp corax
   nlohmann_json::nlohmann_json csv)
target_include_directories(bigrig_obj PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_options(bigrig_obj PRIVATE -Wall -Wextra)

target_compile_options(bigrig_obj PUBLIC $<$<CONFIG:Debug>:-fsanitize=undefined
  --coverage>)
target_link_options(bigrig_obj PUBLIC $<$<CONFIG:Debug>:-fsanitize=undefined
  --coverage>)

target_compile_options(bigrig_interface_obj PRIVATE -Wall -Wextra)

add_executable(bigrig
    main.cpp
)

target_compile_options(bigrig PUBLIC -fdiagnostics-color=always)

target_link_libraries(bigrig 
  bigrig_obj 
  bigrig_interface_obj
  logger 
)

set_target_properties(bigrig
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
)

target_compile_features(bigrig PUBLIC cxx_std_23)
target_compile_options(bigrig PRIVATE -Wall -Wextra)

target_include_directories(bigrig_interface_obj PUBLIC "${CMAKE_SOURCE_DIR}/lib/pcg/include")
target_include_directories(bigrig_obj PUBLIC "${CMAKE_SOURCE_DIR}/lib/pcg/include")
target_include_directories(bigrig PUBLIC "${CMAKE_SOURCE_DIR}/lib/pcg/include")
