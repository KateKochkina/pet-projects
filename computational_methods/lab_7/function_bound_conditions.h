#pragma once

#include "utils.h"

void bound_cond_left_1(TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D,
    TVector<double>& a, Problem_Data data, Options opt, double t_prev, TVector<double>& y_prev);

void bound_cond_left_2(TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D,
    TVector<double>& a, Problem_Data data, Options opt, double t_prev, TVector<double>& y_prev);

void bound_cond_right_1(TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D,
    TVector<double>& a, Problem_Data data, Options opt, double t_prev, TVector<double>& y_prev);

void bound_cond_right_2(TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D,
    TVector<double>& a, Problem_Data data, Options opt, double t_prev, TVector<double>& y_prev);