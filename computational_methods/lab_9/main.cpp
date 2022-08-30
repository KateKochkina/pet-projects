#include <iostream>
#include <string>

#include "methods.h"

int main() {
    double x_left = 0.0;
    double x_right = 1.0;
    double y_left = 0.0;
    double y_right = 1.0;
    double T = 1.0;

    //Problem_data data(x_left, x_right, y_left, y_right, T, f_zero, u_bound_const_0, u_bound_const_0, u_bound_const_0, u_bound_const_0);
    //Problem_data data(x_left, x_right, y_left, y_right, T, f_zero, u_bound_const_minus_1, u_bound_const_1, u_bound_test_2, u_bound_test_2);

    std::string file_name = "test_3/test_3";
    Problem_data data(x_left, x_right, y_left, y_right, T, f_four, u_bound_const_0, u_bound_const_2,
                      u_bound_test_3_left, u_bound_test_3_right);

//    x_right = 0.5;
//    std::string file_name = "var_17/var_17";
//    Problem_data data(x_left, x_right, y_left, y_right, T, f_var_17,
//                      u_bound_var_17_gamma_3, u_bound_var_17_gamma_4,
//                      u_bound_var_17_gamma_1, u_bound_var_17_gamma_2);

//    x_right = 3.14159265359;
//    y_right = 3.14159265359;
//    std::string file_name = "co/checking_order";
//    Problem_data data(x_left, x_right, y_left, y_right, T, ff, uu_gamma, uu_gamma, uu_gamma, uu_gamma);

    /* *
     * g1 = down, g2 = up, g3 = left, g4 = right
     * */
    data.bound[0] = 1;
    data.bound[1] = 1;
    data.bound[2] = 0;
    data.bound[3] = 0;


    for (int n = 25; n < 201; n *= 2) {
        double hx = (x_right - x_left) / n;
        double hy = (y_right - y_left) / n;
        double tau = T / n;

        Options opt(hx, hy, tau, data);

        TVector<TVector<double>> y_prev(opt.size_x_1);
        for (size_t i = 0; i < opt.size_x_1; ++i) {
            y_prev[i].resize(opt.size_x_2);
        }
        TVector<TVector<double>> y_midd(opt.size_x_1);
        for (size_t i = 0; i < opt.size_x_1; ++i) {
            y_midd[i].resize(opt.size_x_2);
        }
        TVector<TVector<double>> y_current(opt.size_x_1);
        for (size_t i = 0; i < opt.size_x_1; ++i) {
            y_current[i].resize(opt.size_x_2);
        }

        TVector<double> A_1(opt.size_x_1);
        TVector<double> B_1(opt.size_x_1);
        TVector<double> C_1(opt.size_x_1);
        TVector<double> D_1(opt.size_x_1);

        TVector<double> A_2(opt.size_x_2);
        TVector<double> B_2(opt.size_x_2);
        TVector<double> C_2(opt.size_x_2);
        TVector<double> D_2(opt.size_x_2);

        method_alternating_directions(file_name, data, opt, y_prev, y_midd, y_current,
                                      A_1, B_1, C_1, D_1, A_2, B_2, C_2, D_2);
    }

    //for (auto i = 0; i < 5; ++i) {
    //    std::cout << "h1 = " << opt.size_x_1 << " h2 = " << opt.size_x_2 << " t = " << opt.size_t << "\n";
    //    method_alternating_directions(file_name + std::to_string(opt.tau), data, opt, y_prev, y_midd, y_current,
    //        A_1, B_1, C_1, D_1, A_2, B_2, C_2, D_2);
    //    std::cout << "cry cry\n";
    //    opt.h_x_1 /= 2;
    //    opt.h_x_2 /= 2;
    //    opt.tau /= 2;
    //    opt.size_x_1 = (int)trunc(data.L_x_1 / opt.h_x_1) + 1;
    //    opt.size_x_2 = (int)trunc(data.L_x_2 / opt.h_x_2) + 1;
    //    opt.size_t = (int)trunc(data.T / opt.tau) + 1;

    //    y_prev.resize(opt.size_x_1);
    //    for (size_t i = 0; i < opt.size_x_1; ++i) {
    //        y_prev[i].resize(opt.size_x_2);
    //    }
    //    y_midd.resize(opt.size_x_1);
    //    for (size_t i = 0; i < opt.size_x_1; ++i) {
    //        y_midd[i].resize(opt.size_x_2);
    //    }
    //    y_current.resize(opt.size_x_1);
    //    for (size_t i = 0; i < opt.size_x_1; ++i) {
    //        y_current[i].resize(opt.size_x_2);
    //    }

    //    A_1.resize(opt.size_x_1);
    //    B_1.resize(opt.size_x_1);
    //    C_1.resize(opt.size_x_1);
    //    D_1.resize(opt.size_x_1);

    //    A_2.resize(opt.size_x_2);
    //    B_2.resize(opt.size_x_2);
    //    C_2.resize(opt.size_x_2);
    //    D_2.resize(opt.size_x_2);
    //}

    //for (auto i = 0; i < 50; ++i) {
    //    method_alternating_directions(file_name + std::to_string(opt.tau), data, opt, y_prev, y_midd, y_current,
    //        A_1, B_1, C_1, D_1, A_2, B_2, C_2, D_2);
    //    opt.tau -= 0.00001;
    //    opt.size_t = (int)trunc(data.T / opt.tau) + 1;
    //}

    std::cout << "All done!\n";
    return 0;
}
