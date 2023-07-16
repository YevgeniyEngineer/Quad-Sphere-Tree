cmake_minimum_required(VERSION 3.18)

project(spherical_grid)

add_library(${PROJECT_NAME} INTERFACE spherical_grid.hpp)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE ${PROJECT_NAME})