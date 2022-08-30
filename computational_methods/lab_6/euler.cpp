#include <iostream>

#include "euler.h"

void euler_explicit(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                    TVector<double> starting_point, Interval time_interval) {

    TVector<double> desired_point(starting_point.get_size());
    double step = time_interval.left_border;

    std::string filename = file_name + "_euler_exp.dat";
    clear_file(filename);
    export_point(filename, step, starting_point);
    while (step < time_interval.right_border) {
        desired_point = starting_point + tau * func(step, starting_point);
        step += tau;
        export_point(filename, step, desired_point);
        starting_point = desired_point;
    }
}

void euler_implicit(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                    TVector<double> starting_point, Interval time_interval) {

    TVector<double> desired_point(starting_point.get_size());
    double step = time_interval.left_border;

    std::string filename = file_name + "_euler_imp.dat";
    clear_file(filename);
    export_point(filename, step, starting_point);
    while (step < time_interval.right_border) {
        desired_point = newton_method(func, function_euler_implict, starting_point, step, tau);
        step += tau;
        export_point(filename, step, desired_point);
        starting_point = desired_point;
    }
}

void symmetric_scheme(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                      TVector<double> starting_point, Interval time_interval) {

    TVector<double> desired_point(starting_point.get_size());
    double step = time_interval.left_border;

    std::string filename = file_name + "_sym_scheme.dat";
    clear_file(filename);
    export_point(filename, step, starting_point);
    while (step < time_interval.right_border) {
        desired_point = newton_method(func, function_symmetric_scheme, starting_point, step, tau);
        step += tau;
        export_point(filename, step, desired_point);
        starting_point = desired_point;
    }
}