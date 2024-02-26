project(OpenSim_analyzeImplementation)

cmake_minimum_required(VERSION 3.2)

# Settings.
# ---------
set(TARGET analyzeImplementation CACHE TYPE STRING)

# OpenSim uses C++11 language features.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find and hook up to OpenSim.
# ----------------------------
find_package(OpenSim REQUIRED PATHS "${OPENSIM_INSTALL_DIR}")

# Configure this project.
# -----------------------
file(GLOB SOURCE_FILES *.h *.cpp)

add_executable(${TARGET} ${SOURCE_FILES})

target_link_libraries(${TARGET} ${OpenSim_LIBRARIES})

# This block copies the additional files into the running directory
# For example vtp, obj files. Add to the end for more extentions
file(GLOB DATA_FILES *.osim *.sto *.vtp *.obj)
foreach(dataFile ${DATA_FILES})
    add_custom_command(
        TARGET ${TARGET}
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy
        ${dataFile}
        ${CMAKE_BINARY_DIR})
endforeach(dataFile)

# Copy project specific files
# file(COPY osimModel.osim;subject01_run5_cycle1_mocoSolution.sto DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
