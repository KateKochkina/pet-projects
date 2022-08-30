#include "lab_tasks.h"

void task_test_var_17() {
    //LINEAR HEAT EQUATION
    std::string file_name = "task_test_var_17/test_var_17_Kx";
    Problem_Data data17_x(BOUND_COND_OF_THE_SECOND_TYPE, BOUND_COND_OF_THE_FIRST_TYPE, P_1,
        u_bound_const, u_init_const, K_x);

    double h = 0.01;
    double tau = 5e-3;
    double sigma = 1.0;
    Options opt_linear(h, tau, sigma, data17_x);

    TVector<double> y_prev(opt_linear.size_x);
    TVector<double> y_next(opt_linear.size_x);

    TVector<double> A(opt_linear.size_x);
    TVector<double> B(opt_linear.size_x);
    TVector<double> C(opt_linear.size_x);
    TVector<double> D(opt_linear.size_x);

    TVector<double> a(opt_linear.size_x);

    linear_heat_eq(file_name, data17_x, opt_linear, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_var17_Kx complete\n";

    //QUASILINEAR HEAT EQUATION
    file_name = "task_test_var_17/test_var_17_Ku";
    Problem_Data data17_u(BOUND_COND_OF_THE_SECOND_TYPE, BOUND_COND_OF_THE_FIRST_TYPE, P_1,
        u_bound_const, u_init_const, K_u_vars);

    Options opt_quasiliner(h, tau, sigma, data17_u);

    y_prev.resize(opt_quasiliner.size_x);
    y_next.resize(opt_quasiliner.size_x);

    A.resize(opt_quasiliner.size_x);
    B.resize(opt_quasiliner.size_x);
    C.resize(opt_quasiliner.size_x);
    D.resize(opt_quasiliner.size_x);

    a.resize(opt_quasiliner.size_x);

    quasilinear_heat_eq(file_name, data17_u, opt_quasiliner, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_var17_Ku linear complete\n";
    quasilinear_heat_eq_non_linear(file_name, data17_u, opt_quasiliner, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_var17_Ku non linear complete\n";
}