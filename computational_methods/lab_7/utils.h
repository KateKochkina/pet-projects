#pragma once

#include "TVector.h"

enum Bound_Cond {
    BOUND_COND_OF_THE_FIRST_TYPE,
    BOUND_COND_OF_THE_SECOND_TYPE,
};

struct Problem_Data;

typedef double(*func1D) (double);
typedef double(*func_init) (double, Problem_Data);

struct Problem_Data {
    double L; //длина стержня
    double c; //теплоемкость
    double rho; //плотность
    double T; //время

    //тип гран. условий
    Bound_Cond left_bound_cond_type;
    Bound_Cond right_bound_cond_type;

    //конкретные функции гран. условий
    func1D left_bound_cond_func;
    func1D right_bound_cond_func;
    func_init init_cond_func;

    //функция коэффициента теплопроводности
    func1D K;

    Problem_Data(double L_, double c_, double rho_, double T_, Bound_Cond left_bound_cond_type_,
        Bound_Cond right_bound_cond_type_, func1D left_bound_cond_func_,
        func1D right_bound_cond_func_, func_init init_cond_func_, func1D K_);

    Problem_Data(Bound_Cond left_bound_cond_type_,
        Bound_Cond right_bound_cond_type_, func1D left_bound_cond_func_,
        func1D right_bound_cond_func_, func_init init_cond_func_, func1D K_);

};

struct Options {
    double h;
    double tau;
    double sigma;
    size_t size_x;
    size_t size_t;

    Options(double h_, double tau_, double sigma_, Problem_Data data);
};


void export_point(const std::string &file_name, const TVector<double> y);

void export_mesh(const std::string &file_name, Options options);

void clear_file(const std::string &file_name);

void right_tridiag_run(TVector<double>& x, TVector<double>& a, TVector<double>& b, TVector<double>& c, TVector<double>& d);
