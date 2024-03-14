cmake_minimum_required(VERSION 3.22)
project(lab4_cpp)

set(CMAKE_CXX_STANDARD 20)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(lab4_cpp main.cpp)
target_link_libraries (lab4_cpp Eigen3::Eigen)
