#pragma once

#include "utils.h"

const double PI = 3.14159265358979323846;

double f1(double x);
double g(double x);
double f2(double x);

//17 var
double f17(double x);
double g17(double x);
double phi17(double t);
double psi17(double t);

//10 var
double f10(double x);
double g10(double x);
double phi10(double t);
double psi10(double t);

//11 var
double f11(double x);
double g11(double x);
double phi11(double t);
double psi11(double t);

//доп
double f_dop_a(double x);
double f_dop_b(double x);
double g_dop_old(double x);

double f_dop(double x);
double g_dop(double x);
double phi_dop(double t);
double f_deriv_zero(double x);