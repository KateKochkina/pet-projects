#pragma once

#include "TVector.h"
#include "utils.h"

TVector<TVector<double>> runge_kutta_4_step_3(TVector<double> func(double, TVector<double>), double tau,
                                              TVector<double> starting_point, Interval time_interval);

void adams_bashforth(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                     TVector<double> starting_point, Interval time_interval);

void predictor_corrector(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                         TVector<double> starting_point, Interval time_interval);