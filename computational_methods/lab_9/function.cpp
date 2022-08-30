#include "function.h"

double f_zero(double x_1, double x_2) {
    return 0;
}

double f_four(double x_1, double x_2) {
    return -4;
}

double u_bound_const_0(double x) {
    return 0;
}

double u_bound_const_1(double x) {
    return 1;
}

double u_bound_const_2(double x) {
    return 2;
}

double u_bound_const_minus_1(double x) {
    return -1;
}

double u_bound_test_2(double x) {
    return 1 + x;
}

double u_bound_test_3_left(double x) {
    return x * x;
}

double u_bound_test_3_right(double x) {
    return 1 + x * x;
}

// вариант 9
double f_var_9(double x_1, double x_2) {
    return -(x_1 * x_1 + x_2 * x_2) * sin(x_1 * x_2);
}

double u_bound_var_9_gamma_1(double x) {
    return x;
}

double u_bound_var_9_gamma_2(double x) {
    return -x * cos(x);
}

//double u_bound_var_9_gamma_3(double x) {
//    return 0;
//}

double u_bound_var_9_gamma_4(double x) {
    return -sin(x);
}

// вариант 10
double f_var_10(double x_1, double x_2) {
    return -2 * (x_2 - x_1);
}

//double u_bound_var_10_gamma_1(double x) {
//    return 0;
//}

double u_bound_var_10_gamma_2(double x) {
    return x * (x - 1);
}

double u_bound_var_10_gamma_3(double x) {
    return x * x;
}

double u_bound_var_10_gamma_4(double x) {
    return 4 * x - x * x;
}

// вариант 11
double f_var_11(double x_1, double x_2) {
    double denom = (x_2 + 1);
    return -(2 * x_1 + 2) / (denom * denom * denom);
}

double u_bound_var_11_gamma_1(double x) {
    return x + 1;
}

double u_bound_var_11_gamma_2(double x) {
    return (x + 1) / 2;
}

double u_bound_var_11_gamma_3(double x) {
    return -1 / (x + 1);
}

double u_bound_var_11_gamma_4(double x) {
    return 1 / (1 + x);
}

// вариант 17
double f_var_17(double x_1, double x_2) {
    return 2 - 2 * x_1;
}

double u_bound_var_17_gamma_1(double x) {
    return 1 - x;
}

double u_bound_var_17_gamma_2(double x) {
    return 2 - 2 * x;
}

double u_bound_var_17_gamma_3(double x) {
    return 1 + x * x;
}

double u_bound_var_17_gamma_4(double x) {
    return -1 - x * x;
}

// проверка порядка
double ff(double x_1, double x_2) {
    return 100.0 * sin(x_1) * cos(x_2);
}

double uu_gamma(double x) {
    return 0.0;
}
