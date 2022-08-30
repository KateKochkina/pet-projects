#pragma once

#include <cmath>

#include "TVector.h"

const double PI = 3.1415926536;

TVector<double> function_test1(double t, TVector<double> x) {
    return {2 * x[0] + x[1] * x[1] - 1,
            6 * x[0] - x[1] * x[1] + 1};
}

TVector<double> function_test2(double t, TVector<double> x) {
    return {1 - x[0] * x[0] - x[1] * x[1],
            2 * x[0]};
}

TVector<double> function_test3(double t, TVector<double> x) {
    double sigma = 10.0;
    double r = 28.0;
    double b = 8.0 / 3.0;
    return {sigma * (x[1] - x[0]),
            x[0] * (r - x[2]) - x[1],
            x[0] * x[1] - b * x[2]};
}

TVector<double> function_pendulum(double t, TVector<double> x) {
    double k = 20.0;
    double m = 0.3;
    return {x[1], - k / m * x[0]};
}

TVector<double> function_analysis(double t, TVector<double> x) {
    return {x[1], t * cos(t)};
}

TVector<double> function_nonlinear_oscillatory_circuit(double t, TVector<double> x) {
    double C0 = 1e-7;
    double L = 1e7;
    double R = 28.0;
    double V0 = 0.1;
    double f = 0.25;
    double C = C0 * sqrt(1 + x[0] / 0.6);
    return {x[1],
            (V0 * sin(2 * PI * f * t) - x[0] - R * C * x[1]) / (L * C)};
}