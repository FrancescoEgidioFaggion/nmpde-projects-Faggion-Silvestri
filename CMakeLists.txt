cmake_minimum_required(VERSION 3.12.0)

project(serial LANGUAGES CXX)

include(../common/cmake-common.cmake)

add_executable(coupled_executable
    src/main.cpp
    src/Stokes.cpp
    src/Darcy.cpp
)

deal_ii_setup_target(coupled_executable)

target_include_directories(coupled_executable PRIVATE ${PROJECT_SOURCE_DIR}/src)
