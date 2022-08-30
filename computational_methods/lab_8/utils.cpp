#include <fstream>

#include "utils.h"

Problem_data::Problem_data(double x_left_, double x_right_, double T_, double a_, func init_func_, func init_deriv_,
    func left_bound_func_, func right_bound_func_) {
    x_left = x_left_;
    x_right = x_right_;
    L = x_right_ - x_left_;
    T = T_;
    a = a_;
    init_func = init_func_;
    init_deriv = init_deriv_;
    left_bound_func = left_bound_func_;
    right_bound_func = right_bound_func_;
}

Options::Options(double h_, double tau_, Problem_data data) {
    h = h_;
    tau = tau_;
    size_x = (int)trunc(data.L / h) + 1;
    size_t = (int)trunc(data.T / tau) + 1;
}

void export_point(const std::string &file_name, const TVector<double> y) {
    std::ofstream file("/Users/artem/Desktop/study/вычи/3/export/" + file_name, std::ios::app);
    for (size_t i = 0; i < y.get_size(); ++i) {
        file << std::fixed << std::setprecision(20) << y[i] << "    ";
    }
    file << std::endl;
    file.close();
    file.clear();
}

void export_mesh(const std::string &file_name, Options options, double x_left) {
    std::ofstream file("/Users/artem/Desktop/study/вычи/3/export/" + file_name, std::ios::app);
    for (size_t i = 0; i < options.size_x; ++i) {
        file << std::fixed << std::setprecision(10) << x_left + options.h * i << "    ";
    }
    file << std::endl;

    for (size_t j = 0; j < options.size_t; ++j) {
        file << std::fixed << std::setprecision(10) << options.tau * j << "    ";
    }
    file.close();
    file.clear();
}

void clear_file(const std::string &file_name) {
    std::ofstream file("/Users/artem/Desktop/study/вычи/3/export/" + file_name, std::ios::trunc);
    file.close();
    file.clear();
}

double DD(func f, double x, Options opt) {
    return (f(x - opt.h) - 2 * f(x) + f(x + opt.h)) / (opt.h * opt.h);
}