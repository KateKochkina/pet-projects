#include "lab_tasks.h"

void task_test_1() {
    double x_left = 0;
    double L = 1.0;
    double T = 1.0;
    double a = 1.0;
    Problem_data data1(x_left, L, T, a, f1, g, g, g);

    double h = 0.001;
    double tau = 0.001;
    Options opt1(h, tau, data1);

    TVector<double> y_prev(opt1.size_x);
    TVector<double> y_current(opt1.size_x);
    TVector<double> y_next(opt1.size_x);

    std::string file_name = "task_test_1/test_1";
    scheme_cross(file_name, data1, opt1, y_prev, y_current, y_next);

    std::cout << "test 1 ok!\n";
}