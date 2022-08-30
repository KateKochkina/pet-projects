#include "lab_tasks.h"

void task_quantity() {
    std::string file_name = "task_quantity/test";
    Problem_Data data(BOUND_COND_OF_THE_SECOND_TYPE, BOUND_COND_OF_THE_FIRST_TYPE, P_1,
        u_bound_const, u_init_const, K_u_vars);

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

    for (size_t i = 0; i < 3; ++i) {
        double tol = 0.01;
        quasilinear_heat_eq_non_linear_iter(file_name, data, opt, y_prev, y_next, A, B, C, D, a, tol / pow(10,i));
        std::cout << "tol = " << tol / pow(10, i) << " complete\n";
    }
}
