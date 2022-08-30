#include "lab_tasks.h"

void task_test_monotony() {
    std::string file_name = "task_test_monotony/test";
    Problem_Data data(1, 1, 1, 1, BOUND_COND_OF_THE_SECOND_TYPE, BOUND_COND_OF_THE_SECOND_TYPE,
                      u_bound_zero, u_bound_zero, u_init_mon, K_mon);
    double h = 0.01;
    double tau = 5e-3;
    double sigma = 1.0;
    Options opt(h, tau, sigma, data);

    TVector<double> y_prev(opt.size_x);
    TVector<double> y_next(opt.size_x);

    TVector<double> A(opt.size_x);
    TVector<double> B(opt.size_x);
    TVector<double> C(opt.size_x);
    TVector<double> D(opt.size_x);

    TVector<double> a(opt.size_x);

    linear_heat_eq(file_name + "_s_1", data, opt, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_monotony sigma 1 complete\n";

    opt.sigma = 0.5;

    linear_heat_eq(file_name + "_s_0.5", data, opt, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_monotony sigma 0.5 complete\n";
}