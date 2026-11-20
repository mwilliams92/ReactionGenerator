
cmake_minimum_required(VERSION 3.16)
project(BeamMC LANGUAGES CXX)

option(BEAMMC_WITH_ROOT "Enable ROOT histogram and TTree output" OFF)

add_executable(BeamMC
    "${CMAKE_CURRENT_SOURCE_DIR}/BeamMC.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ExcitationPopulation.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/ReactionCrossSection.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/TwoBodyKinematics.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeState.cpp"
)

target_compile_features(BeamMC PRIVATE cxx_std_11)

target_include_directories(BeamMC PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# Optional ROOT support
if(BEAMMC_WITH_ROOT)
    find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree)
    target_compile_definitions(BeamMC PRIVATE BEAMMC_WITH_ROOT)
    target_link_libraries(BeamMC PRIVATE
        ROOT::Core ROOT::RIO ROOT::Hist ROOT::Tree
    )
endif()

# Copy input and SRIM data into the build directory
foreach(data_dir IN ITEMS User_Input SRIM_Files)
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${data_dir}"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
endforeach()

# Create output directory without copying its contents
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Output_Files")