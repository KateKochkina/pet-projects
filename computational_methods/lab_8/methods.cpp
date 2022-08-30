#include <iostream>
#include <string>
#include "methods.h"

void scheme_cross(std::string file_name, Problem_data data, Options opt,
                  TVector<double>& y_prev, TVector<double>& y_current,
                  TVector<double>& y_next, func DD_init_func) {
    clear_file(file_name + "_mesh.dat");
    export_mesh(file_name + "_mesh.dat", opt, data.x_left);

    //начальные условия - решение на нулевом временном слое
    y_prev[0] = data.left_bound_func(0.0);
    y_prev[opt.size_x - 1] = data.right_bound_func(0.0);
    for (size_t i = 1; i < opt.size_x - 1; ++i) {
        y_prev[i] = data.init_func(data.x_left + i * opt.h);
    }
    clear_file(file_name + "_sol.dat");
    export_point(file_name + "_sol.dat", y_prev);

    //решение на первом слое
    double coeff = data.a * data.a * opt.tau * opt.tau / 2.0;
    y_current[0] = data.left_bound_func(opt.tau);
    y_current[opt.size_x - 1] = data.right_bound_func(opt.tau);
    if (DD_init_func == nullptr) {
        for (size_t i = 1; i < opt.size_x - 1; ++i) {
            y_current[i] = y_prev[i] + opt.tau * data.init_deriv(data.x_left + i * opt.h) + coeff * DD(data.init_func, data.x_left + i * opt.h, opt);
            //std::cout << std::fixed << std::setprecision(20) << " i = " << data.x_left + i * opt.h << "\nderiv = " << data.init_deriv(data.x_left + i * opt.h) << "\n";
        }
    } else {
        for (size_t i = 1; i < opt.size_x - 1; ++i) {
            y_current[i] = y_prev[i] + opt.tau * data.init_deriv(data.x_left + i * opt.h) + coeff * DD_init_func(data.x_left + i * opt.h);
        }
    }
    export_point(file_name + "_sol.dat", y_current);

    coeff = data.a * data.a * opt.tau * opt.tau / (opt.h * opt.h);
    for (size_t j = 2; j < opt.size_t; ++j) {
        y_next[0] = data.left_bound_func(j * opt.tau);
        y_next[opt.size_x - 1] = data.right_bound_func(j * opt.tau);
        for (size_t n = 1; n < opt.size_x - 1; ++n) {
            y_next[n] = 2.0 * (1.0 - coeff) * y_current[n] + coeff * (y_current[n + 1] + y_current[n - 1]) - y_prev[n];
        }
        export_point(file_name + "_sol.dat", y_next);
        y_prev = y_current;
        y_current = y_next;
    }
}

//проверка порядка
void order_test(std::string file_name, Problem_data data, Options opt, double q,
    TVector<double>& y_prev, TVector<double>& y_current, TVector<double>& y_next) {
    size_t old_size_x = opt.size_x;
    for (size_t i = 0; i < 5; ++i) {
        scheme_cross(file_name + "_t_" + std::to_string(opt.tau) + "_h_" + std::to_string(opt.h), data, opt, y_prev, y_current, y_next);
        opt.tau = opt.tau / (q);
        opt.h = opt.h / (q);
        opt.size_x = (int)trunc(data.L / opt.h) + 1;
        opt.size_t = (int)trunc(data.T / opt.tau) + 1;
        y_prev.resize(opt.size_x);
        y_current.resize(opt.size_x);
        y_next.resize(opt.size_x);
    }
    y_prev.resize(old_size_x);
    y_current.resize(old_size_x);
    y_next.resize(old_size_x);
}
