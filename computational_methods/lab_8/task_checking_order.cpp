#include "lab_tasks.h"

void task_test_checking_order() {
    double x_left = 0;
    double L = 1.0;
    double T = 1.0;
    double a = 1.0;
    Problem_data data1(x_left, L, T, a, f1, g, g, g);
    Problem_data data2(x_left, L, T, a, f2, g, g, g);
    Problem_data data17(x_left, L, T, a, f17, g17, phi17, psi17);

    double h = 0.1;
    double tau = 0.0025;
    Options opt1(h, tau, data1);
    Options opt2(h, tau, data2);
    Options opt17(h, tau, data17);

    TVector<double> y_prev(opt1.size_x);
    TVector<double> y_current(opt1.size_x);
    TVector<double> y_next(opt1.size_x);

    std::string file_name = "task_test_checking_order/test_1";
    double q = 2.0;
    order_test(file_name, data1, opt1, q, y_prev, y_current, y_next);

    //file_name = "task_test_checking_order/test2";
    //order_test(file_name, data2, opt2, q, y_prev, y_current, y_next);

    //file_name = "task_test_checking_order/var17";
    //order_test(file_name, data17, opt17, q, y_prev, y_current, y_next);

    std::cout << "checking order ok!\n";
}