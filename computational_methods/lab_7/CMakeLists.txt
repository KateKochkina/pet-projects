cmake_minimum_required(VERSION 3.17)
project(2)

set(CMAKE_CXX_STANDARD 17)

add_executable(2 main.cpp utils.cpp initial_and_boundary_conditions.cpp function_bound_conditions.cpp methods.cpp task_checking_order.cpp task_monotony.cpp task_test_1_2.cpp task_test_3.cpp task_test_var_17.cpp task_quantity.cpp)