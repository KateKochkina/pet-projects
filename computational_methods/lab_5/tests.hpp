#pragma once

#include <cmath>

#include "utils.h"

double function_test1(double x) {
    return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

double function_test1_d(double x) {
    return 0.121495 - 1.5119 * x + 5.9535 * std::pow(x, 2) - 9.28 * std::pow(x, 3) + 5. * std::pow(x, 4);
}

double function_test1_dd(double x) {
    return -1.5119 + 11.907 * x - 27.84 * std::pow(x, 2) + 20. * std::pow(x, 3);
}

double function_test2(double x) {
    return std::sqrt(x + 1.) - 1.;
}

double function_test2_d(double x) {
    return 1. / (2. * std::sqrt(x + 1.));
}

double function_test2_dd(double x) {
    return -0.25 * std::pow(x + 1., -1.5);
}

double function_test3(double x) {
    return 35. * std::pow(x, 3) - 67. * std::pow(x, 2) - 3. * x + 3.;
}

double function_test3_d(double x) {
    return 105. * std::pow(x, 2) - 134. * x - 3.;
}

double function_test3_dd(double x) {
    return 210. * x - 134.;
}

std::vector<double> function_test4(const std::vector<double> &x) {
    return { x[0] * x[0] - x[1] * x[1] - 15., x[0] * x[1] + 4. };
}

static double first4(double x, double y) {
    return (4. + x * y) / (x * x + y * y);
}

static double second4(double x, double y) {
    return (-15. + x * x - y * y) / (2. * (x * x + y * y));
}

std::vector<double> function_test4_newton_an(const std::vector<double> &x) {
    auto a = first4(x[0], x[1]);
    auto b = second4(x[0], x[1]);
    return { x[1] * a + x[0] * b,
             x[0] * a - x[1] * b };
}

std::vector<double> function_test5(const std::vector<double> &x) {
    return { x[0] * x[0] + x[1] * x[1] + x[0] + x[1] - 8.,
             x[0] * x[0] + x[1] * x[1] + x[0] * x[1] - 7. };
}

static double denominator(double x, double y) {
    return (x - y) * (-1. + 2. * x + 2. * y);
}

std::vector<double> function_test5_newton_an(const std::vector<double> &x) {
    auto den = denominator(x[0], x[1]);
    return { (-1. + x[0]) * (-7. + x[0] + x[0] * x[0] + 2. * x[1] - x[1] * x[1]) / den,
             (-1. + x[1]) * (7. - 2. * x[0] + x[0] * x[0] - x[1] - x[1] * x[1]) / den };
}

double function_variant20(double x) {
    return sinh((std::sqrt(13.) * x * x * x - 9. * x - 5. - std::sqrt(17.)) / 10.) +
           tan((x * x + x + std::pow(2., 1. / 3.)) / (3. * x - 5.)) + 0.6;
}

std::vector<double> system_variant20(const std::vector<double> &x) {
    return {6. * x[0] * x[0] - 7. * x[1] * x[1] - 65.,
            13. * x[0] * x[1] + 41.};
}

static double denominator20(double x, double y) {
    return 156. * x * x + 182. * y * y;
}

std::vector<double> system_variant20_newton_an(const std::vector<double> &x) {
    auto den = denominator20(x[0], x[1]);
    return { (574. * x[1] - 845. * x[0]) / den + x[0] / 2,
             (492. * x[0] + 845. * x[1]) / den + x[1] / 2 };
}