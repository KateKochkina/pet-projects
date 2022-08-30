#include <fstream>

#include "utils.h"

Problem_data::Problem_data(double x_1_left_, double x_1_right_, double x_2_left_, double x_2_right_, double T_,
    func_right_side f_, func left_bound_func_, func right_bound_func_, func down_bound_func_, func up_bound_func_, size_t size_bound) {

    x_1_left = x_1_left_;
    x_1_right = x_1_right_;
    x_2_left = x_2_left_;
    x_2_right = x_2_right_;
    L_x_1 = x_1_right - x_1_left;
    L_x_2 = x_2_right - x_2_left;
    T = T_;
    f = f_;
    left_bound_func = left_bound_func_;
    right_bound_func = right_bound_func_;
    down_bound_func = down_bound_func_;
    up_bound_func = up_bound_func_;
    bound.resize(size_bound);
}

Options::Options(double h_x_1_, double h_x_2_, double tau_, Problem_data data) {
    h_x_1 = h_x_1_;
    h_x_2 = h_x_2_;
    tau = tau_;
    size_x_1 = (int)trunc(data.L_x_1 / h_x_1) + 1;
    size_x_2 = (int)trunc(data.L_x_2 / h_x_2) + 1;
    size_t = (int)trunc(data.T / tau) + 1;
}

void export_point(const std::string &file_name, const TVector<TVector<double>> y) {
    std::ofstream file("/Users/artem/Desktop/study/вычи/4/export/" + file_name, std::ios::app);
    for (size_t i = 0; i < y.get_size(); ++i) {
        for (size_t j = 0; j < y[i].get_size(); ++j) {
            file << std::fixed << std::setprecision(20) << y[i][j] << "    ";
        }
        file << std::endl;
    }
    //file << std::endl;
    file.close();
    file.clear();
}

void export_mesh(const std::string &file_name, Options options) {
    std::ofstream file("/Users/artem/Desktop/study/вычи/4/export/" + file_name, std::ios::app);
    for (size_t i = 0; i < options.size_x_1; ++i) {
        file << std::fixed << std::setprecision(10) << options.h_x_1 * i << "    ";
    }
    file << std::endl;

    for (size_t j = 0; j < options.size_x_2; ++j) {
        file << std::fixed << std::setprecision(10) << options.h_x_2 * j << "    ";
    }
    file << std::endl;

    for (size_t j = 0; j < options.size_t; ++j) {
        file << std::fixed << std::setprecision(10) << options.tau * j << "    ";
    }
    file.close();
    file.clear();
}

void clear_file(const std::string &file_name) {
    std::ofstream file("/Users/artem/Desktop/study/вычи/4/export/" + file_name, std::ios::trunc);
    file.close();
    file.clear();
}

double DD_1(TVector<double> y, size_t i, Options opt) {
    return (y[i - 1] - 2 * y[i] + y[i + 1]) / (opt.h_x_1 * opt.h_x_1);
}

double DD_2(TVector<double> y, size_t j, Options opt) {
    return (y[j - 1] - 2 * y[j] + y[j + 1]) / (opt.h_x_2 * opt.h_x_2);
}

//double DD(func f, double x, Options opt) {
//    return (f(x - opt.h_x_2) - 2 * f(x) + f(x + opt.h_x_2)) / (opt.h_x_2 * opt.h_x_2);
//}

void right_tridiag_run(TVector<double>& x, TVector<double>& a, TVector<double>& b, TVector<double>& c, TVector<double>& d) {
    size_t size = x.get_size();
    TVector<double> alpha(size - 1);
    TVector<double> beta(size - 1);

    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];
    for (auto i = 1; i < size - 1; ++i) {
        alpha[i] = c[i] / (b[i] - a[i] * alpha[i - 1]);
        beta[i] = (d[i] + a[i] * beta[i - 1]) / (b[i] - a[i] * alpha[i - 1]);
    }
    //std::cout << alpha[size - 2] << " " << beta[size - 2] << "\n";
    x[size - 1] = (d[size - 1] + a[size - 1] * beta[size - 2]) / (b[size - 1] - a[size - 1] * alpha[size - 2]);
    for (int i = size - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
}