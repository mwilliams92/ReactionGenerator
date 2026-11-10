cmake_minimum_required(VERSION 3.12)
project(BeamMC VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BEAMMC_ENABLE_ROOT "Build with ROOT histogram/tree output" ON)

add_executable(BeamMC
    BeamMC.cpp
    src/TwoBodyKinematics.cpp
    target_input.dat
    beam_input.txt
    degrader_input.txt
    reaction_input.txt
)

target_include_directories(BeamMC PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(BeamMC PRIVATE -Wall -Wextra -Wpedantic)

if(BEAMMC_ENABLE_ROOT)
    find_package(ROOT REQUIRED COMPONENTS Core Hist Tree)
    target_compile_definitions(BeamMC PRIVATE BEAMMC_WITH_ROOT)
    target_include_directories(BeamMC PRIVATE ${ROOT_INCLUDE_DIRS})
    target_link_libraries(BeamMC PRIVATE ROOT::Core ROOT::Hist ROOT::Tree)
endif()

install(TARGETS BeamMC RUNTIME DESTINATION bin)