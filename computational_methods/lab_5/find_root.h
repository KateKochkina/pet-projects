#pragma once

#include "utils.h"

//! ONE-DIMENSIONAL
Meta bisection_method(const std::function<double (double)> &f, const Interval &interval, double precision = 1e-6);

Meta newton_method(const std::function<double (double)> &f, const std::function<double (double)> &f_d,
        double prev, const Interval &interval, double precision = 1e-6);

Meta newton_method_numerical(const std::function<double (double)> &f,
        double prev, const Interval &interval, double precision = 1e-6);

//! requires function value at interval points
Meta newton_method_modification(const std::function<double (double)> &f, const std::function<double (double)> &f_d,
        const std::function<double (double)> &f_dd, double prev, const Interval &interval, double precision = 1e-6);

//! requires function value at interval points
Meta newton_method_modification_numerical(const std::function<double (double)> &f, double prev,
        const Interval &interval, double precision = 1e-6);

//! requires function value at interval points
double get_chord_method_approx(const Interval &interval);

//! TWO-DIMENSIONAL
Meta newton_method(const std::function<std::vector<double> (std::vector<double>)> &inv_jacobian,
        std::vector<double> prev, const std::vector<Interval> &intervals, double precision = 1e-6);

Meta newton_method_numerical(const std::function<std::vector<double>(std::vector<double>)> &f,
        std::vector<double> prev, const std::vector<Interval> &intervals, double precision = 1e-6);

void export_newton(const std::string &file_name, 
        const std::function<std::vector<double>(std::vector<double>)> &f,
        const std::vector<Interval> &intervals);

Meta bisection_method_test(const std::function<double (double)> &f,
                           const Interval &interval, double precision = 1e-6);

Meta newton_method_test(const std::function<double (double)> &f,
                        double prev, const Interval &interval, double precision = 1e-6);