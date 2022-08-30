#include "lab_tasks.h"

void task_test_var_17() {
    double x_left = 0;
    double L = 1.0;
    double T = 1.0;
    double a = 1.0;
    Problem_data data17(x_left, L, T, a, f17, g17, phi17, psi17);
    double h = 0.1;
    double tau = 0.0025;
    Options opt17(h, tau, data17);

    TVector<double> y_prev(opt17.size_x);
    TVector<double> y_current(opt17.size_x);
    TVector<double> y_next(opt17.size_x);

    std::string file_name = "task_test_var_17/test_var_17";
    scheme_cross(file_name, data17, opt17, y_prev, y_current, y_next);

//    double q = 2.0;
//    order_test(file_name, data17, opt17, q, y_prev, y_current, y_next);

    std::cout << "var 17 ok!\n";
}