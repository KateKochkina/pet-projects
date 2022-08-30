#include "lab_tasks.h"

void task_dop1() {
    double x_left = -2.0;
    double x_right = 2.0;
    double T = 1.0;
    double a = 1.0;
    Problem_data data(x_left, x_right, T, a, f_dop, g, g, g);
    double h = 0.1;
    double tau = 0.01;
    Options opt01(h, tau, data);
    h = 0.1;
    tau = 0.05;
    Options opt05(h, tau, data);
    h = 0.1;
    tau = 0.075;
    Options opt075(h, tau, data);
    h = 0.1;
    tau = 0.1;
    Options opt1(h, tau, data);

    TVector<double> y_prev(opt01.size_x);
    TVector<double> y_current(opt01.size_x);
    TVector<double> y_next(opt01.size_x);

    std::string file_name = "dop1/test_dop1_01";
    scheme_cross(file_name, data, opt01, y_prev, y_current, y_next);

    file_name = "dop1/test_dop1_05";
    y_prev.resize(opt05.size_x);
    y_current.resize(opt05.size_x);
    y_next.resize(opt05.size_x);
    scheme_cross(file_name, data, opt05, y_prev, y_current, y_next);

    file_name = "dop1/test_dop1_075";
    y_prev.resize(opt075.size_x);
    y_current.resize(opt075.size_x);
    y_next.resize(opt075.size_x);
    scheme_cross(file_name, data, opt075, y_prev, y_current, y_next);

    file_name = "dop1/test_dop1_1";
    y_prev.resize(opt1.size_x);
    y_current.resize(opt1.size_x);
    y_next.resize(opt1.size_x);
    scheme_cross(file_name, data, opt1, y_prev, y_current, y_next);

    // производная аналитически
    file_name = "dop1/test_dop1_anal_01";
    y_prev.resize(opt01.size_x);
    y_current.resize(opt01.size_x);
    y_next.resize(opt01.size_x);
    scheme_cross(file_name, data, opt01, y_prev, y_current, y_next, f_deriv_zero);

    file_name = "dop1/test_dop1_anal_05";
    y_prev.resize(opt05.size_x);
    y_current.resize(opt05.size_x);
    y_next.resize(opt05.size_x);
    scheme_cross(file_name, data, opt05, y_prev, y_current, y_next, f_deriv_zero);

    file_name = "dop1/test_dop1_anal_075";
    y_prev.resize(opt075.size_x);
    y_current.resize(opt075.size_x);
    y_next.resize(opt075.size_x);
    scheme_cross(file_name, data, opt075, y_prev, y_current, y_next, f_deriv_zero);

    file_name = "dop1/test_dop1_anal_1";
    y_prev.resize(opt1.size_x);
    y_current.resize(opt1.size_x);
    y_next.resize(opt1.size_x);
    scheme_cross(file_name, data, opt1, y_prev, y_current, y_next, f_deriv_zero);
}


void task_dop2() {
    double x_left = -1.0;
    double x_right = 1.0;
    double T = 10.0;
    double a = 1.0;
    Problem_data data(x_left, x_right, T, a, g, g_dop, g, g);
    double h = 0.1;
    double tau = 0.01;
    Options opt01(h, tau, data);
    h = 0.1;
    tau = 0.05;
    Options opt05(h, tau, data);
    h = 0.1;
    tau = 0.075;
    Options opt075(h, tau, data);
    h = 0.01;
    tau = 0.01;
    Options opt1(h, tau, data);

    TVector<double> y_prev(opt01.size_x);
    TVector<double> y_current(opt01.size_x);
    TVector<double> y_next(opt01.size_x);

    std::string file_name = "dop2/test_dop2_01";
    scheme_cross(file_name, data, opt01, y_prev, y_current, y_next);

    file_name = "dop2/test_dop2_05";
    y_prev.resize(opt05.size_x);
    y_current.resize(opt05.size_x);
    y_next.resize(opt05.size_x);
    scheme_cross(file_name, data, opt05, y_prev, y_current, y_next);

    file_name = "dop2/test_dop2_075";
    y_prev.resize(opt075.size_x);
    y_current.resize(opt075.size_x);
    y_next.resize(opt075.size_x);
    scheme_cross(file_name, data, opt075, y_prev, y_current, y_next);

    file_name = "dop2/test_dop2_1";
    y_prev.resize(opt1.size_x);
    y_current.resize(opt1.size_x);
    y_next.resize(opt1.size_x);
    scheme_cross(file_name, data, opt1, y_prev, y_current, y_next);
}


void task_dop3() {
    double x_left = 0.0;
    double x_right = 4.0 * PI;
    double T = 20.0;
    double a = 1.0;
    Problem_data data(x_left, x_right, T, a, g, g, phi_dop, g);
    double h = 0.1;
    double tau = 0.01;
    Options opt01(h, tau, data);
    h = 0.1;
    tau = 0.05;
    Options opt05(h, tau, data);
    h = 0.1;
    tau = 0.075;
    Options opt075(h, tau, data);
    h = 0.01;
    tau = 0.01;
    Options opt1(h, tau, data);

    TVector<double> y_prev(opt01.size_x);
    TVector<double> y_current(opt01.size_x);
    TVector<double> y_next(opt01.size_x);

    std::string file_name = "dop3/test_dop3_01";
    scheme_cross(file_name, data, opt01, y_prev, y_current, y_next);

    file_name = "dop3/test_dop3_05";
    y_prev.resize(opt05.size_x);
    y_current.resize(opt05.size_x);
    y_next.resize(opt05.size_x);
    scheme_cross(file_name, data, opt05, y_prev, y_current, y_next);

    file_name = "dop3/test_dop3_075";
    y_prev.resize(opt075.size_x);
    y_current.resize(opt075.size_x);
    y_next.resize(opt075.size_x);
    scheme_cross(file_name, data, opt075, y_prev, y_current, y_next);

    file_name = "dop3/test_dop3_1";
    y_prev.resize(opt1.size_x);
    y_current.resize(opt1.size_x);
    y_next.resize(opt1.size_x);
    scheme_cross(file_name, data, opt1, y_prev, y_current, y_next);
}
