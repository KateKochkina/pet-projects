cmake_minimum_required(VERSION 3.17)
project(1)

set(CMAKE_CXX_STANDARD 17)

add_executable(1 main.cpp euler.cpp multistep_methods.cpp runge_kutta.cpp utils.cpp)