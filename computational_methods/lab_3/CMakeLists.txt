cmake_minimum_required(VERSION 3.17)
project(3)

set(CMAKE_CXX_STANDARD 17)

add_executable(3 main.cpp constants.h matrix.h norm.h print.h gaussMethod.h qrMethod.h inverseMatrix.h input.h memory.h vector.h eigen.h inter.cpp)