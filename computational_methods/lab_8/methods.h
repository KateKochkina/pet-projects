#pragma once

#include "function.h"

void scheme_cross(std::string file_name, Problem_data data, Options opt,
    TVector<double> &y_prev, TVector<double> &y_current, TVector<double> &y_next, func DD_init_func = nullptr);

void order_test(std::string file_name, Problem_data data, Options opt, double q,
    TVector<double> &y_prev, TVector<double> &y_current, TVector<double> &y_next);