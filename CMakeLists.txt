cmake_minimum_required(VERSION 3.10)
project(NUFEB-dev)

set(CMAKE_CXX_STANDARD 11)

# Set source directories
set(SRC_DIR_NUFEB ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SRC_DIR_LAMMPS ${CMAKE_CURRENT_SOURCE_DIR}/lammps_stable_29Oct2020/src)

include_directories(${SRC_DIR_NUFEB})
include_directories(${SRC_DIR_LAMMPS})

# Add source files from the first directory
file(GLOB SOURCES_NUFEB "${SRC_DIR_NUFEB}/*.cpp" "${SRC_DIR_NUFEB}/*.h")

# Add source files from the second directory
file(GLOB SOURCES_LAMMPS "${SRC_DIR_LAMMPS}/*.cpp" "${SRC_DIR_LAMMPS}/*.h")

# Combine source files from both directories
set(SOURCES ${SOURCES_NUFEB} ${SOURCES_LAMMPS})

# Create executable
add_executable(NUFEB-dev ${SOURCES})