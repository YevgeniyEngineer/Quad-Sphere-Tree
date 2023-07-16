cmake_minimum_required(VERSION 3.18)

project(quad_sphere_tree)

add_library(${PROJECT_NAME} INTERFACE quad_sphere_tree.hpp)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE ${PROJECT_NAME})