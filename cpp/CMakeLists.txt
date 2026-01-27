cmake_minimum_required(VERSION 3.15)
project(abp_core)

include(FetchContent)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.11.1
)

FetchContent_MakeAvailable(pybind11)


pybind11_add_module(abp_core src/core.cpp)
