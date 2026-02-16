cmake_minimum_required(VERSION 3.20)
project(CurvedABP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

add_executable(core 
  src/core.cpp
  src/vec2.cpp
  src/state.cpp
  src/dynamics.cpp
  src/h5sim.cpp
)

target_include_directories(core PRIVATE
    ${HDF5_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(core PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
