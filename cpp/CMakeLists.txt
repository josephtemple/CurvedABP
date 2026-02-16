cmake_minimum_required(VERSION 3.20)
project(CurvedABP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

include_directories(${HDF5_INCLUDE_DIRS})

add_executable(core src/core.cpp)

target_link_libraries(core PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
