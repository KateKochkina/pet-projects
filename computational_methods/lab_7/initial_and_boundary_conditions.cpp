#include "initial_and_boundary_conditions.h"

double u_init_zero(double x, Problem_Data data) {
    return 0.0;
}

double u_init_const(double x, Problem_Data data) {
    double u0 = 0.5;
    return u0;
}

double u_init_test1(double x, Problem_Data data) {
    double u0 = 0.5;
    return u0 + x * (data.L - x);
}

double u_init_my(double x, Problem_Data data) {
    double omega = 5.0;
    return 2 * sin(omega * x);
}

double u_init_mon(double x, Problem_Data data) {
    if (x <= 0.5 && x >= 0) {
        return 0;
    } else {
        return 1;
    }
}

//функции гран. усл.
double u_bound_zero(double t) {
    return 0.0;
}

double u_bound_const(double t) {
    double u0 = 0.5;
    return u0;
}

double u_bound_test3(double t) {
    double sigma = 2.0;
    double kappa0 = 0.5;
    double c = 5.0;
    double u0 = pow(sigma * c * c / kappa0, 1 / sigma);
    return u0 * pow(t, 1 / sigma);
}

double u_bound_my(double t) {
    double omega = 5.0;
    return 2 * exp(-omega * t);
}

double P_1(double t) {
    double t0 = 0.5;
    double Q = 10.0;
    if (t >= 0.0 && t < t0) {
        return Q;
    } else {
        return 0.0;
    }
}

double P_2(double t) {
    double t0 = 0.5;
    double Q = 10.0;
    if (t >= 0.0 && t < t0) {
        return 2 * Q * t;
    } else {
        return 0.0;
    }
}

//коэффициент теплопроводности
double K_x(double x) {
    double k1 = 0.5;
    double k2 = 1.5;
    double x1 = 0.5;
    double x2 = 0.6;
    if (x >= 0.0 && x <= x1) {
        return k1;
    } else if (x > x1&& x < x2) {
        return k1 * (x - x2) / (x1 - x2) + k2 * (x - x1) / (x2 - x1);
    } else {
        return k2;
    }
}

double K_u(double u) {
    double kappa0 = 0.5;
    double sigma = 2.0;
    return kappa0 * pow(u, sigma);
}

double K_u_vars(double u) {
    double alpha = 2.0;
    double beta = 1.5;
    double gamma = 2.0;
    return alpha + beta * pow(u, gamma);
}

double K_my(double x) {
    double omega = 5.0;
    return 1 / omega;
}

double K_mon(double x) {
    double c = 5.0;
    return c;
}
