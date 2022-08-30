#pragma once

#include "TVector.h"
#include "utils.h"
typedef TVector<double> (*MYFUNC)(double, TVector<double>);
void euler_explicit(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                    TVector<double> starting_point, Interval time_interval);

void euler_implicit(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                    TVector<double> starting_point, Interval time_interval);

void symmetric_scheme(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                      TVector<double> starting_point, Interval time_interval);