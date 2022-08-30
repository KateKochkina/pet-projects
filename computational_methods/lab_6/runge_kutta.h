#pragma once

#include "TVector.h"
#include "utils.h"

void runge_kutta_2_const_step(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                              TVector<double> starting_point, Interval time_interval);

void runge_kutta_2_variable_step(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                                 TVector<double> starting_point, Interval time_interval, Options options);

void runge_kutta_4_const_step(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                              TVector<double> starting_point, Interval time_interval);

void runge_kutta_4_variable_step(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                                 TVector<double> starting_point, Interval time_interval, Options options);