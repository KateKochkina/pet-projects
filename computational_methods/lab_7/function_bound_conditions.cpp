#include "function_bound_conditions.h"

void bound_cond_left_1(TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D,
    TVector<double>& a, Problem_Data data, Options opt, double t_prev, TVector<double>& y_prev) {

    A[0] = 0.0;
    B[0] = 1.0;
    C[0] = 0.0;
    D[0] = data.left_bound_cond_func(t_prev + opt.tau);
}

void bound_cond_left_2(TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D,
    TVector<double>& a, Problem_Data data, Options opt, double t_prev, TVector<double>& y_prev) {

    A[0] = 0.0;
    B[0] = data.c * data.rho * opt.h * opt.h + 2 * opt.sigma * opt.tau * a[1];
    C[0] = 2 * opt.sigma * opt.tau * a[1];
    D[0] = 2 * opt.sigma * opt.tau * opt.h * data.left_bound_cond_func(t_prev + opt.tau) +
        2 * (1 - opt.sigma) * opt.tau * (opt.h * data.left_bound_cond_func(t_prev) + a[1] * (y_prev[1] - y_prev[0])) +
        data.c * data.rho * opt.h * opt.h * y_prev[0];
}

void bound_cond_right_1(TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D,
    TVector<double>& a, Problem_Data data, Options opt, double t_prev, TVector<double>& y_prev) {

    A[opt.size_x - 1] = 0.0;
    B[opt.size_x - 1] = 1.0;
    C[opt.size_x - 1] = 0.0;
    D[opt.size_x - 1] = data.right_bound_cond_func(t_prev + opt.tau);
}

void bound_cond_right_2(TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D,
    TVector<double>& a, Problem_Data data, Options opt, double t_prev, TVector<double>& y_prev) {

    A[opt.size_x - 1] = 2 * opt.sigma * opt.tau * a[opt.size_x - 1];
    B[opt.size_x - 1] = data.c * data.rho * opt.h * opt.h + 2 * opt.sigma * opt.tau * a[opt.size_x - 1];
    C[opt.size_x - 1] = 0.0;
    D[opt.size_x - 1] = 2 * opt.sigma * opt.tau * opt.h * data.right_bound_cond_func(t_prev + opt.tau) +
        2 * (1 - opt.sigma) * opt.tau * (opt.h * data.right_bound_cond_func(t_prev) - a[opt.size_x - 1] *
        (y_prev[opt.size_x - 1] - y_prev[opt.size_x - 2])) +
        data.c * data.rho * opt.h * opt.h * y_prev[opt.size_x - 1];
}

