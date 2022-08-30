#pragma once

#include "function_bound_conditions.h"
#include "initial_and_boundary_conditions.h"

typedef void(*Bound_cond_func) (TVector<double>&, TVector<double>&, TVector<double>&,
    TVector<double>&, TVector<double>&, Problem_Data, Options, double, TVector<double>&);

void function_copy_pasta(const std::string &file_name, Problem_Data data, Options options, TVector<double>& y_prev,
    Bound_cond_func &bound_cond_func_left, Bound_cond_func &bound_cond_func_right);

void linear_heat_eq(const std::string &file_name, Problem_Data data, Options options, TVector<double>& y_prev,
    TVector<double>& y_next, TVector<double>& A, TVector<double>& B,
    TVector<double>& C, TVector<double>& D, TVector<double>& a);

void quasilinear_heat_eq(const std::string &file_name, Problem_Data data, Options options,
    TVector<double>& y_prev, TVector<double>& y_next,
    TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D, TVector<double>& a);

void quasilinear_heat_eq_non_linear(const std::string &file_name, Problem_Data data, Options options,
    TVector<double>& y_prev, TVector<double>& y_next,
    TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D, TVector<double>& a);

void quasilinear_heat_eq_non_linear_iter(const std::string &file_name, Problem_Data data, Options options,
    TVector<double>& y_prev, TVector<double>& y_next,
    TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D, TVector<double>& a, double tol);