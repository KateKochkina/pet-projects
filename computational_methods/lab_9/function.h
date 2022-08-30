#pragma once

#include "utils.h"

const double PI = 3.14159265358979323846;

double f_zero(double x_1, double x_2);

double f_four(double x_1, double x_2);

double u_bound_const_0(double x);

double u_bound_const_1(double x);

double u_bound_const_2(double x);

double u_bound_const_minus_1(double x);

double u_bound_test_2(double x);

double u_bound_test_3_left(double x);

double u_bound_test_3_right(double x);

// вариант 9
double f_var_9(double x_1, double x_2);

double u_bound_var_9_gamma_1(double x);

double u_bound_var_9_gamma_2(double x);

//double u_bound_var_9_gamma_3(double x);

double u_bound_var_9_gamma_4(double x);

// вариант 10
double f_var_10(double x_1, double x_2);

//double u_bound_var_10_gamma_1(double x);

double u_bound_var_10_gamma_2(double x);

double u_bound_var_10_gamma_3(double x);

double u_bound_var_10_gamma_4(double x);

// вариант 11
double f_var_11(double x_1, double x_2);

double u_bound_var_11_gamma_1(double x);

double u_bound_var_11_gamma_2(double x);

double u_bound_var_11_gamma_3(double x);

double u_bound_var_11_gamma_4(double x);

// вариант 17
double f_var_17(double x_1, double x_2);

double u_bound_var_17_gamma_1(double x);

double u_bound_var_17_gamma_2(double x);

double u_bound_var_17_gamma_3(double x);

double u_bound_var_17_gamma_4(double x);

// проверка порядка
double ff(double x_1, double x_2);
double uu_gamma(double x);
