#include <fstream>

#include "utils.h"

Problem_Data::Problem_Data(double L_, double c_, double rho_, double T_, Bound_Cond left_bound_cond_type_,
    Bound_Cond right_bound_cond_type_, func1D left_bound_cond_func_,
    func1D right_bound_cond_func_, func_init init_cond_func_, func1D K_) {

    L = L_;
    c = c_;
    rho = rho_;
    T = T_;
    left_bound_cond_type = left_bound_cond_type_;
    right_bound_cond_type = right_bound_cond_type_;
    left_bound_cond_func = left_bound_cond_func_;
    right_bound_cond_func = right_bound_cond_func_;
    init_cond_func = init_cond_func_;
    K = K_;
}

Problem_Data::Problem_Data(Bound_Cond left_bound_cond_type_,
    Bound_Cond right_bound_cond_type_, func1D left_bound_cond_func_,
    func1D right_bound_cond_func_, func_init init_cond_func_, func1D K_) {

    L = 1.0;
    c = 1.0;
    rho = 0.75;
    T = 1.0;
    left_bound_cond_type = left_bound_cond_type_;
    right_bound_cond_type = right_bound_cond_type_;
    left_bound_cond_func = left_bound_cond_func_;
    right_bound_cond_func = right_bound_cond_func_;
    init_cond_func = init_cond_func_;
    K = K_;
}

Options::Options(double h_, double tau_, double sigma_, Problem_Data data) {
    h = h_;
    tau = tau_;
    sigma = sigma_;
    size_x = (int)trunc(data.L / h) + 1;
    size_t = (int)trunc(data.T / tau) + 1;
}

void export_point(const std::string &file_name, const TVector<double> y) {
    std::ofstream file("/Users/artem/Desktop/study/вычи/2/export/" + file_name, std::ios::app);
    for (size_t i = 0; i < y.get_size(); ++i) {
        file << std::fixed << std::setprecision(20) << y[i] << "    ";
    }
    file << std::endl;
    file.close();
    //file.clear();
}

void export_mesh(const std::string &file_name, Options options) {
    std::ofstream file("/Users/artem/Desktop/study/вычи/2/export/" + file_name, std::ios::app);
    for (size_t i = 0; i < options.size_x; ++i) {
        file << std::fixed << std::setprecision(10) << options.h * i << "    ";
    }
    file << std::endl;

    for (size_t j = 0; j < options.size_t; ++j) {
        file << std::fixed << std::setprecision(10) << options.tau * j << "    ";
    }
    file.close();
    file.clear();
}

void clear_file(const std::string &file_name) {
    std::ofstream file("/Users/artem/Desktop/study/вычи/2/export/" + file_name, std::ios::trunc);
    file.close();
    file.clear();
}


void right_tridiag_run(TVector<double>& x, TVector<double>& a, TVector<double>& b, TVector<double>& c, TVector<double>& d) {
    size_t size = x.get_size();
    TVector<double> alpha(size);
    TVector<double> beta(size);

    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];
    for (auto i = 1; i < size - 1; ++i) {
        alpha[i] = c[i] / (b[i] - a[i] * alpha[i - 1]);
        beta[i] = (d[i] + a[i] * beta[i - 1]) / (b[i] - a[i] * alpha[i - 1]);
    }
    x[size - 1] = (d[size - 1] + a[size - 1] * beta[size - 2]) / (b[size - 1] - a[size - 1] * alpha[size - 2]);
    for (int i = size - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
}