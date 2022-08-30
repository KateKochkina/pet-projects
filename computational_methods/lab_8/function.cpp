#include "function.h"

double f1(double x) {
    return sin(PI * x);
}

double g(double x) {
    return 0.0;
}

double f2(double x) {
    return x * (1 - x);
}

//17 var
double f17(double x) {
    return (x + 0.5) * (x + 0.5);
}

double g17(double x) {
    return (x + 1.0) * sin(x);
}

double phi17(double t) {
    return 1.0 + 0.4 * t;
}

double psi17(double t) {
    return 2.25;
}

//10 var
double f10(double x) {
    return (x + 1) * sin(PI * x);
}

double g10(double x) {
    return x * (x + 1);
}

double phi10(double t) {
    return 0.0;
}

double psi10(double t) {
    return 0.5 * t;
}

double f10_deriv(double x) {
    return PI * (2 * cos(PI * x) - PI * (x + 1) * sin(PI * x));
}

//11 var
double f11(double x) {
    return (1 - x) * cos(PI * x / 2);
}

double g11(double x) {
    return 2 * x + 1;
}

double phi11(double t) {
    return 2 * t + 1;
}

double psi11(double t) {
    return 0.0;
}

//доп
double f_dop_a(double x) {
    if (x >= -1 && x <= 1)
        return 1.0;
    else
        return 0.0;
}

double f_dop_b(double x) {
    if (x >= -1 && x <= 0)
        return 1 + x;
    else if (x > 0 && x <= 1)
        return 1 - x;
    else
        return 0.0;
}

double g_dop_old(double x) {
    return 0.0;
}

double f_dop(double x) {
    if ((x >= -1.0 / 3) && (x <= 1.0 / 3)) {
        return 1.0;
    } else {
        return 0.0;
    }
}

double g_dop(double x) {
    if ((-1.0 / 3.0 + 0.00000001 <= x) && (x <= 1.0 / 3.0 + 0.00000001)) {
        return 1.0;
    } else {
        return 0.0;
    }
}

double phi_dop(double t) {
    return asin(sin(t));
}

double f_deriv_zero(double x) {
    return 0.0;
}
