#pragma once

#include "function.h"

void method_alternating_directions(std::string file_name, Problem_data data, Options opt,
    TVector<TVector<double>> &y_prev, TVector<TVector<double>> &y_midd, TVector<TVector<double>> &y_current,
    TVector<double> &A_1, TVector<double> &B_1, TVector<double> &C_1, TVector<double> &D_1,
    TVector<double> &A_2, TVector<double> &B_2, TVector<double> &C_2, TVector<double> &D_2);