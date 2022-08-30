#pragma once

#include "TVector.h"

struct Interval {
    Interval(double left_border, double right_border);

    double left_border;
    double right_border;
};

struct Options {
    Options();
    Options(double fac_min, double fac_max, double fac, bool step_recount);

    double fac_min;
    double fac_max;
    double fac;
    bool step_recount;
};

void export_point(const std::string& file_name, const double step, const TVector<double> point);

void export_tau(const std::string& file_name, const double step, double tau, double error);

void clear_file(const std::string& file_name);

TVector<double> function_euler_implict(double t, TVector<double> x, TVector<double> starting_point,
                                       TVector<double> func(double, TVector<double>), double tau);

TVector<double> function_symmetric_scheme(double t, TVector<double> x, TVector<double> starting_point,
                                          TVector<double> func(double, TVector<double>), double tau);

TVector<TVector<double>> get_jacobian(TVector<double> func(double, TVector<double>),
                                      TVector<double> func_meth(double, TVector<double>, TVector<double>, TVector<double> (*)(double, TVector<double>), double),
                                      TVector<double>& x, TVector<double>& starting_point, double t, double tau, double precision = 1e-6);

TVector<TVector<double>> inverse_jacobian(TVector<TVector<double>>& jacobian);

TVector<double> newton_method(TVector<double> func(double, TVector<double>),
                              TVector<double> func_meth(double, TVector<double>, TVector<double>, TVector<double> (*)(double, TVector<double>), double),
                              TVector<double>& starting_point, double t, double tau, double precision = 1e-6);