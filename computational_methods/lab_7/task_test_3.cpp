#include "lab_tasks.h"

void task_test_3() {
    std::string file_name = "task_test_3/test_3";
    double L = 10.0;
    double c = 1.0;
    double rho = 1.0;
    double T = 2.0;
    Problem_Data data(L, c, rho, T, BOUND_COND_OF_THE_FIRST_TYPE, BOUND_COND_OF_THE_SECOND_TYPE,
        u_bound_test3, u_bound_zero, u_init_zero, K_u);
    double h = 0.2;
    double tau = 2e-4;
    double sigma = 1.0;
    Options opt(h, tau, sigma, data);

    TVector<double> y_prev(opt.size_x);
    TVector<double> y_next(opt.size_x);

    TVector<double> A(opt.size_x);
    TVector<double> B(opt.size_x);
    TVector<double> C(opt.size_x);
    TVector<double> D(opt.size_x);

    TVector<double> a(opt.size_x);

    quasilinear_heat_eq(file_name, data, opt, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_3 linear complete\n";
    quasilinear_heat_eq_non_linear(file_name, data, opt, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_3 non linear complete\n";
}