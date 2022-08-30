cmake_minimum_required(VERSION 3.17)
project(5)

set(CMAKE_CXX_STANDARD 17)

add_executable(5 main.cpp Gauss_Method.cpp SLAE_solve.cpp Integral_equation_solver.cpp Test_functions.cpp)