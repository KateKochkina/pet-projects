#pragma once

#include "utils.h"

double u_init_zero(double x, Problem_Data data);

double u_init_const(double x, Problem_Data data);

double u_init_test1(double x, Problem_Data data);

double u_init_my(double x, Problem_Data data);

double u_init_mon(double x, Problem_Data data);

//функции гран. усл.
double u_bound_zero(double t);

double u_bound_const(double t);

double u_bound_test3(double t);

double u_bound_my(double t);

double P_1(double t);

double P_2(double t);

//коэффициент теплопроводности
double K_x(double x);

double K_u(double u);

double K_u_vars(double u);

double K_my(double x);

double K_mon(double x);
