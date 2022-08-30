#include "lab_tasks.h"

void task_test_2() {
    double x_left = 0;
    double L = 1.0;
    double T = 1.0;
    double a = 1.0;
    Problem_data data2(x_left, L, T, a, f2, g, g, g);

    double h = 0.001;
    double tau = 0.001;
    Options opt2(h, tau, data2);

    TVector<double> y_prev(opt2.size_x);
    TVector<double> y_current(opt2.size_x);
    TVector<double> y_next(opt2.size_x);

    std::string file_name = "task_test_2/test_2";
    scheme_cross(file_name, data2, opt2, y_prev, y_current, y_next);

    std::cout << "test 2 ok!\n";
}