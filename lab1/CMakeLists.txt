cmake_minimum_required(VERSION 3.22)
project(lab1)

set(CMAKE_CXX_STANDARD 20)
include_directories(/usr/include/eigen3)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(lab1 main.cpp)
target_link_libraries (lab1 Eigen3::Eigen)
