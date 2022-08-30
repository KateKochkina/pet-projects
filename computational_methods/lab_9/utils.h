#pragma once

#include "TVector.h"

typedef double(*func) (double);
typedef double(*func_right_side) (double, double);

struct Problem_data {
    double x_1_left;
    double x_1_right;
    double x_2_left;
    double x_2_right;
    double L_x_1, L_x_2;
    double T;
    func_right_side f;
    func left_bound_func;
    func right_bound_func;
    func down_bound_func;
    func up_bound_func;
    TVector<int> bound;

    Problem_data(double x_1_left_, double x_1_right_, double x_2_left_, double x_2_right_, double T_,
        func_right_side f_, func left_bound_func_, func right_bound_func_, func down_bound_func_, func up_bound_func_, size_t size_bound = 4);
};

struct Options {
    double h_x_1, h_x_2;
    double tau;
    size_t size_x_1;
    size_t size_x_2;
    size_t size_t;

    Options(double h_x_1_, double h_x_2_, double tau_, Problem_data data);
};


void export_point(const std::string &file_name, const TVector<TVector<double>> y);

void export_mesh(const std::string &file_name, Options options);

void clear_file(const std::string &file_name);

double DD_1(TVector<double> y, size_t i, Options opt);

double DD_2(TVector<double> y, size_t j, Options opt);

void right_tridiag_run(TVector<double>& x, TVector<double>& a, TVector<double>& b, TVector<double>& c, TVector<double>& d);