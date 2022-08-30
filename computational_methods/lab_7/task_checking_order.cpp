#include "lab_tasks.h"

const double PI = 3.14159265358979323846;

void task_checking_order() {
    std::string file_name = "task_checking_order/my_test";
    double omega = 5.0;
    Problem_Data data(PI / (2 * omega), 1.0, 1.0, 1.0, BOUND_COND_OF_THE_FIRST_TYPE, BOUND_COND_OF_THE_FIRST_TYPE,
        u_bound_zero, u_bound_my, u_init_my, K_my);

    double h = 0.01;
    double tau = 0.005;
    double sigma = 1.0;
    Options opt(h, tau, sigma, data);

    TVector<double> y_prev(opt.size_x);
    TVector<double> y_next(opt.size_x);

    TVector<double> A(opt.size_x);
    TVector<double> B(opt.size_x);
    TVector<double> C(opt.size_x);
    TVector<double> D(opt.size_x);

    TVector<double> a(opt.size_x);

    double q = 2;
//    std::cout << "checking_order\n";
//    for (auto i = 0; i < 5; ++i) {
//    	opt.sigma = 1.0;
//        linear_heat_eq(file_name + "_t_" + std::to_string(opt.tau) + "_h_" + std::to_string(opt.h) + "_s_" + std::to_string(opt.sigma),
//            data, opt, y_prev, y_next, A, B, C, D, a);
//        std::cout << "s = 1, h = " + std::to_string(opt.h) + ", t = " + std::to_string(opt.tau) + "\n";
//
//    	opt.sigma = 0.0;
//        linear_heat_eq(file_name + "_t_" + std::to_string(opt.tau) + "_h_" + std::to_string(opt.h) + "_s_" + std::to_string(opt.sigma),
//            data, opt, y_prev, y_next, A, B, C, D, a);
//        std::cout << "s = 0, h = " + std::to_string(opt.h) + ", t = " + std::to_string(opt.tau) + "\n";
//
////    	opt.tau = opt.tau / (q * q);
////    	opt.h = opt.h / q;
//    	opt.tau = opt.tau / 2.0;
//    	opt.h = opt.h / sqrt(2.0);
//    	opt.size_x = (int)trunc(data.L / opt.h) + 1;
//    	opt.size_t = (int)trunc(data.T / opt.tau) + 1;
//    	y_prev.resize(opt.size_x);
//        y_next.resize(opt.size_x);
//        A.resize(opt.size_x);
//        B.resize(opt.size_x);
//        C.resize(opt.size_x);
//        D.resize(opt.size_x);
//        a.resize(opt.size_x);
//    }

    opt.h = 0.01;
    opt.tau = 0.005;
    opt.sigma = 0.5;
    opt.size_x = (int)trunc(data.L / opt.h) + 1;
    opt.size_t = (int)trunc(data.T / opt.tau) + 1;
    y_prev.resize(opt.size_x);
    y_next.resize(opt.size_x);
    A.resize(opt.size_x);
    B.resize(opt.size_x);
    C.resize(opt.size_x);
    D.resize(opt.size_x);
    a.resize(opt.size_x);
    for (auto i = 0; i < 5; ++i) {
        
        linear_heat_eq(file_name + "_t_" + std::to_string(opt.tau) + "_h_" + std::to_string(opt.h) + "_s_" + std::to_string(opt.sigma),
            data, opt, y_prev, y_next, A, B, C, D, a);
        std::cout << "s = 0.5, h = " + std::to_string(opt.h) + ", t = " + std::to_string(opt.tau) + "\n";

//        opt.tau = opt.tau / (q * q);
//        opt.h = opt.h / (q * q);
        opt.tau = opt.tau / 2.0;
        opt.h = opt.h / 2.0;
        opt.size_x = (int)trunc(data.L / opt.h) + 1;
        opt.size_t = (int)trunc(data.T / opt.tau) + 1;
        y_prev.resize(opt.size_x);
        y_next.resize(opt.size_x);
        A.resize(opt.size_x);
        B.resize(opt.size_x);
        C.resize(opt.size_x);
        D.resize(opt.size_x);
        a.resize(opt.size_x);
    }
}