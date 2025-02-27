cmake_minimum_required(VERSION 3.10)

project(dt)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find OpenMesh
find_package(OpenMesh REQUIRED)
include_directories(${OPENMESH_INCLUDE_DIRS})

# Add executable
add_executable(dt main.cpp)

# Link libraries
target_link_libraries(dt ${EIGEN3_LIBS} ${OPENMESH_INCLUDE_DIRS} OpenMeshCore OpenMeshTools)

# Set the output directory for the build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
