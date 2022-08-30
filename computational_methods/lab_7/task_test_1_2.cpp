#include "lab_tasks.h"

void task_test_1_2() {
    std::string file_name1 = "task_test_1_2/test_1";
    Problem_Data data(BOUND_COND_OF_THE_FIRST_TYPE, BOUND_COND_OF_THE_FIRST_TYPE, 
        u_bound_const, u_bound_const, u_init_test1, K_x);
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

    linear_heat_eq(file_name1 + "_s_1", data, opt, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_1 sigma 1 complete\n";

    std::string file_name2 = "task_test_1_2/test_2";
    Problem_Data data2(BOUND_COND_OF_THE_SECOND_TYPE, BOUND_COND_OF_THE_SECOND_TYPE,
        u_bound_zero, u_bound_zero, u_init_test1, K_x);
    linear_heat_eq(file_name2 + "_s_1", data2, opt, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_2 sigma 1 complete\n";

    opt.h = 0.01;
    opt.tau = 5e-3;
    opt.sigma = 0.5;

    linear_heat_eq(file_name1 + "_s_0.5", data, opt, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_1 sigma 0.5 complete\n";

    linear_heat_eq(file_name2 + "_s_0.5", data2, opt, y_prev, y_next, A, B, C, D, a);
    std::cout << "test_2 sigma 0.5 complete\n";
}