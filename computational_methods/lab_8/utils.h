#pragma once

#include "TVector.h"

typedef double(*func) (double);

struct Problem_data {
    double x_left;
    double x_right;
    double L;
    double T;
    double a;
    func init_func;
    func init_deriv;
    func left_bound_func;
    func right_bound_func;

    Problem_data(double x_left_, double x_right_, double T_, double a_, func init_func_, func init_deriv_,
        func left_bound_func_, func right_bound_func_);
};

struct Options {
    double h;
    double tau;
    size_t size_x;
    size_t size_t;

    Options(double h_, double tau_, Problem_data data);
};

void export_point(const std::string &file_name, const TVector<double> y);

void export_mesh(const std::string &file_name, Options options, double x_left = 0.0);

void clear_file(const std::string &file_name);

double DD(func f, double x, Options opt);