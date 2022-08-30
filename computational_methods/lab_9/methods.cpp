#include <iostream>
#include <string>
#include "methods.h"

double func_anal(double x_1, double x_2) {
    // test 3
    return x_1 * x_1 + x_2 * x_2;
    // проверка порядка
//    return 100.0 * exp(-1.0) * sin(x_1) * cos(x_2);
}

void method_alternating_directions(std::string file_name, Problem_data data, Options opt,
            TVector<TVector<double>> &y_prev, TVector<TVector<double>> &y_midd, TVector<TVector<double>> &y_current,
            TVector<double> &A_1, TVector<double> &B_1, TVector<double> &C_1, TVector<double> &D_1,
            TVector<double> &A_2, TVector<double> &B_2, TVector<double> &C_2, TVector<double> &D_2) {

//    clear_file(file_name + "_mesh.dat");
//    export_mesh(file_name + "_mesh.dat", opt);

    // начальный слой
    for (size_t i = 0; i < opt.size_x_1; ++i) {
        for (size_t j = 0; j < opt.size_x_2; ++j) {
            y_prev[i][j] = 1.0;
                /*(data.x_1_left + i * opt.h_x_1) * (data.x_1_left + i * opt.h_x_1) + (data.x_2_left + j * opt.h_x_2);*/
                /*(data.x_1_left + i * opt.h_x_1 + 1) * (-(data.x_2_left + j * opt.h_x_2) / 2 + 1);*/
                /*cos(PI * (data.x_1_left + i * opt.h_x_1)) * sin(PI * (data.x_2_left + j * opt.h_x_2));*/
                /*100;*/
        }
    }

//    clear_file(file_name + "_sol_0.dat");
//    export_point(file_name + "_sol_0.dat", y_prev);

    TVector<TVector<double>> phi(opt.size_x_1);
    for (size_t i = 0; i < opt.size_x_1; ++i) {
        phi[i].resize(opt.size_x_2);
        for (size_t j = 0; j < opt.size_x_2; ++j) {
            phi[i][j] = data.f(data.x_1_left + i * opt.h_x_1, data.x_2_left + j * opt.h_x_2);
        }
    }

    // заполнение коэф. для 1й прогонки
    A_1[0] = 0.0;
    B_1[0] = -(1 - data.bound[0]) + 2.0 * data.bound[0] * (1 / (opt.h_x_1 * opt.h_x_1) + 1 / opt.tau);
    C_1[0] = 2.0 * data.bound[0] / (opt.h_x_1 * opt.h_x_1);

    for (size_t i = 1; i < opt.size_x_1 - 1; ++i) {
        A_1[i] = 1 / (opt.h_x_1 * opt.h_x_1);
        B_1[i] = 2.0 * (1 / (opt.h_x_1 * opt.h_x_1) + 1 / opt.tau);
        C_1[i] = 1 / (opt.h_x_1 * opt.h_x_1);
    }

    A_1[opt.size_x_1 - 1] = 2.0 * data.bound[1] / (opt.h_x_1 * opt.h_x_1);;
    B_1[opt.size_x_1 - 1] = -(1 - data.bound[1]) + 2.0 * data.bound[1] * (1 / (opt.h_x_1 * opt.h_x_1) + 1 / opt.tau);
    C_1[opt.size_x_1 - 1] = 0.0;

    // заполнение коэф. для 2й прогонки
    A_2[0] = 0.0;
    B_2[0] = -(1 - data.bound[2]) + 2.0 * data.bound[2] * (1 / (opt.h_x_2 * opt.h_x_2) + 1 / opt.tau);
    C_2[0] = 2.0 * data.bound[2] / (opt.h_x_2 * opt.h_x_2);

    for (size_t j = 1; j < opt.size_x_2 - 1; ++j) {
        A_2[j] = 1 / (opt.h_x_2 * opt.h_x_2);
        B_2[j] = 2.0 * (1 / (opt.h_x_2 * opt.h_x_2) + 1 / opt.tau);
        C_2[j] = 1 / (opt.h_x_2 * opt.h_x_2);
    }

    A_2[opt.size_x_2 - 1] = 2.0 * data.bound[3] / (opt.h_x_2 * opt.h_x_2);
    B_2[opt.size_x_2 - 1] = -(1 - data.bound[3]) + 2.0 * data.bound[3] * (1 / (opt.h_x_2 * opt.h_x_2) + 1 / opt.tau);
    C_2[opt.size_x_2 - 1] = 0;

    TVector<double> vec_1(opt.size_x_1); // вектор для прогонки
    TVector<double> vec_2(opt.size_x_2); // вектор для определения D_2
    size_t iter = 0;
    double err;
    double err_mean = 0.0;
    double eps = 0.0001;
    for (size_t k = 1; k < opt.size_t; ++k) {
//    do {
//        iter++;
        // промежуточный слой
        if (data.bound[2] == 0) {
            for (size_t i = 0; i < opt.size_x_1; ++i) {
                y_midd[i][0] = data.down_bound_func(data.x_1_left + i * opt.h_x_1);
            }
        } else {
            D_1[0] = -(1 - data.bound[0]) * data.left_bound_func(data.x_2_left)
                + data.bound[0] * (2 / opt.tau * y_prev[0][0] + /*DD_2(y_prev[0], 1, opt)*/2 * (y_prev[0][1] - y_prev[0][0]) / (opt.h_x_2 * opt.h_x_2) + phi[0][0]
                    + 2 * data.left_bound_func(data.x_2_left) / opt.h_x_1 + 2 * data.down_bound_func(data.x_1_left) / opt.h_x_2);
            for (size_t i = 1; i < opt.size_x_1 - 1; ++i) {
                D_1[i] = 2 / opt.tau * y_prev[i][0] + /*DD_2(y_prev[i], 1, opt)*/ 2 * (y_prev[i][1] - y_prev[i][0]) / (opt.h_x_2 * opt.h_x_2) + phi[i][0]
                    + 2 * data.down_bound_func(data.x_1_left + i * opt.h_x_1) / opt.h_x_2;
            }
            D_1[opt.size_x_1 - 1] = -(1 - data.bound[1]) * data.right_bound_func(data.x_2_left)
                + data.bound[1] * (2 / opt.tau * y_prev[opt.size_x_1 - 1][0] + /*DD_2(y_prev[opt.size_x_1 - 1], 1, opt)*/2 * (y_prev[opt.size_x_1 - 1][1] - y_prev[opt.size_x_1 - 1][0]) / (opt.h_x_2 * opt.h_x_2)
                    + phi[opt.size_x_1 - 1][0] + 2 * data.right_bound_func(data.x_2_left) / opt.h_x_1
                    + 2 * data.down_bound_func(data.x_1_right) / opt.h_x_2);

            right_tridiag_run(vec_1, A_1, B_1, C_1, D_1);
            for (size_t i = 0; i < opt.size_x_1; ++i) {
                y_midd[i][0] = vec_1[i];
            }
        }

        for (size_t j = 1; j < opt.size_x_2 - 1; ++j) {
            D_1[0] = -(1 - data.bound[0]) * data.left_bound_func(data.x_2_left + j * opt.h_x_2)
                + data.bound[0] * (2 / opt.tau * y_prev[0][j] + DD_2(y_prev[0], j, opt) + phi[0][j]
                    + 2 * data.left_bound_func(data.x_2_left + j * opt.h_x_2) / opt.h_x_1);
            for (size_t i = 1; i < opt.size_x_1 -1 ; ++i) {
                D_1[i] = 2 / opt.tau * y_prev[i][j] + DD_2(y_prev[i], j, opt) + phi[i][j];
            }
            D_1[opt.size_x_1 - 1] = -(1 - data.bound[1]) * data.right_bound_func(data.x_2_left + j * opt.h_x_2)
                + data.bound[1] * (2 / opt.tau * y_prev[opt.size_x_1 - 1][j] + DD_2(y_prev[opt.size_x_1 - 1], j, opt) +
                    phi[opt.size_x_1 - 1][j] + 2 * data.right_bound_func(data.x_2_left + j * opt.h_x_2) / opt.h_x_1);

            right_tridiag_run(vec_1, A_1, B_1, C_1, D_1);
            for (size_t i = 0; i < opt.size_x_1; ++i) {
                y_midd[i][j] = vec_1[i];
            }
        }
        //for (size_t j = 0; j < opt.size_x_2; ++j) {
        //    y_midd[0][j] = data.left_bound_func(data.x_2_left + j * opt.h_x_2);
        //    y_midd[opt.size_x_1 - 1][j] = data.right_bound_func(data.x_2_left + j * opt.h_x_2);
        //}

        if (data.bound[3] == 0) {
            for (size_t i = 0; i < opt.size_x_1; ++i) {
                y_midd[i][opt.size_x_2 - 1] = data.up_bound_func(data.x_1_left + i * opt.h_x_1);
            }
        } else {
            D_1[0] = -(1 - data.bound[0]) * data.left_bound_func(data.x_2_right)
                + data.bound[0] * (2 / opt.tau * y_prev[0][opt.size_x_2 - 1] + /*DD_2(y_prev[0], opt.size_x_2 - 2, opt)*/2 * (y_prev[0][opt.size_x_2 - 2] - y_prev[0][opt.size_x_2 - 1]) / (opt.h_x_2 * opt.h_x_2) + phi[0][opt.size_x_2 - 1]
                    + 2 * data.left_bound_func(data.x_2_right) / opt.h_x_1 + 2 * data.up_bound_func(data.x_1_left) / opt.h_x_2);
            for (size_t i = 1; i < opt.size_x_1 - 1; ++i) {
                D_1[i] = 2 / opt.tau * y_prev[i][opt.size_x_2 - 1] + /*DD_2(y_prev[i], opt.size_x_2 - 2, opt)*/2 * (y_prev[i][opt.size_x_2 - 2] - y_prev[i][opt.size_x_2 - 1]) / (opt.h_x_2 * opt.h_x_2) + phi[i][opt.size_x_2 - 1]
                    + 2 * data.up_bound_func(data.x_1_left + i * opt.h_x_1) / opt.h_x_2;
            }
            D_1[opt.size_x_1 - 1] = -(1 - data.bound[1]) * data.right_bound_func(data.x_2_right)
                + data.bound[1] * (2 / opt.tau * y_prev[opt.size_x_1 - 1][opt.size_x_2 - 1]
                    + /*DD_2(y_prev[opt.size_x_1 - 1], opt.size_x_2 - 2, opt)*/2 * (y_prev[opt.size_x_1 - 1][opt.size_x_2 - 2] - y_prev[opt.size_x_1 - 1][opt.size_x_2 - 1]) / (opt.h_x_2 * opt.h_x_2) + phi[opt.size_x_1 - 1][opt.size_x_2 - 1]
                    + 2 * data.right_bound_func(data.x_2_right) / opt.h_x_1 + 2 * data.up_bound_func(data.x_1_right) / opt.h_x_2);

            right_tridiag_run(vec_1, A_1, B_1, C_1, D_1);
            for (size_t i = 0; i < opt.size_x_1; ++i) {
                y_midd[i][opt.size_x_2 - 1] = vec_1[i];
            }
        }
        //if (k == 1) {
        //    std::cout << " qqq = " << k << "\n";
        //    for (size_t i = 0; i < opt.size_x_1; ++i) {
        //        std::cout << " i = " << i << "\n" << y_midd[i];
        //    }
        //}

        // следующий слой
        if (data.bound[0] == 0) {
            for (size_t j = 0; j < opt.size_x_2; ++j) {
                y_current[0][j] = data.left_bound_func(data.x_2_left + j * opt.h_x_2);
            }
        } else {
            //for (size_t n = 0; n < opt.size_x_1; ++n) {
            //    vec_2[n] = y_midd[n][0];
            //}
            D_2[0] = -(1 - data.bound[2]) * data.down_bound_func(data.x_1_left)
                + data.bound[2] * (2 / opt.tau * y_midd[0][0] + /*DD_1(vec_2, 1, opt)*/ 2 * (y_midd[1][0] - y_midd[0][0]) / (opt.h_x_1 * opt.h_x_1)
                    + phi[0][0] +
                    2 * data.down_bound_func(data.x_1_left) / opt.h_x_2 + 2 * data.left_bound_func(data.x_2_left) / opt.h_x_1);
            for (size_t j = 1; j < opt.size_x_2 - 1; ++j) {
                //for (size_t n = 0; n < opt.size_x_1; ++n) {
                //    vec_2[n] = y_midd[n][j];
                //}
                D_2[j] = 2 / opt.tau * y_midd[0][j] + /*DD_1(vec_2, 1, opt)*/ 2 * (y_midd[1][j] - y_midd[0][j]) / (opt.h_x_1 * opt.h_x_1)
                    + phi[0][j] + 2 * data.left_bound_func(data.x_2_left + j * opt.h_x_2) / opt.h_x_1;
            }
            //for (size_t n = 0; n < opt.size_x_1; ++n) {
            //    vec_2[n] = y_midd[n][opt.size_x_2 - 1];
            //}
            D_2[opt.size_x_2 - 1] = -(1 - data.bound[3]) * data.up_bound_func(data.x_1_left)
                + data.bound[3] * (2 / opt.tau * y_midd[0][opt.size_x_2 - 1] + /*DD_1(vec_2, 1, opt)*/
                    2 * (y_midd[1][opt.size_x_2 - 1] - y_midd[0][opt.size_x_2 - 1]) / (opt.h_x_1 * opt.h_x_1) + phi[0][opt.size_x_2 - 1]
                    + 2 * data.up_bound_func(data.x_1_left) / opt.h_x_2 + 2 * data.left_bound_func(data.x_2_right) / opt.h_x_1);

            right_tridiag_run(y_current[0], A_2, B_2, C_2, D_2);
        }

        for (size_t i = 1; i < opt.size_x_1 - 1; ++i) {
            //for (size_t n = 0; n < opt.size_x_1; ++n) {
            //    vec_2[n] = y_midd[n][0];
            //}
            D_2[0] = -(1 - data.bound[2]) * data.down_bound_func(data.x_1_left + i * opt.h_x_1)
                + data.bound[2] * (2 / opt.tau * y_midd[i][0] + /*DD_1(vec_2, i, opt)*/
                (y_midd[i - 1][0] - 2 * y_midd[i][0] + y_midd[i + 1][0]) / (opt.h_x_1 * opt.h_x_1) + phi[i][0] +
                    2 * data.down_bound_func(data.x_1_left + i * opt.h_x_1) / opt.h_x_2);
            for (size_t j = 1; j < opt.size_x_2 - 1; ++j) {
                //for (size_t n = 0; n < opt.size_x_1; ++n) {
                //    vec_2[n] = y_midd[n][j];
                //}
                D_2[j] = 2 / opt.tau * y_midd[i][j] + /*DD_1(vec_2, i, opt)*/ (y_midd[i - 1][j] - 2 * y_midd[i][j] + y_midd[i + 1][j])
                    / (opt.h_x_1 * opt.h_x_1) + phi[i][j];
            }
            //for (size_t n = 0; n < opt.size_x_1; ++n) {
            //    vec_2[n] = y_midd[n][opt.size_x_2 - 1];
            //}
            D_2[opt.size_x_2 - 1] = -(1 - data.bound[3]) * data.up_bound_func(data.x_1_left + i * opt.h_x_1)
                + data.bound[3] * (2 / opt.tau * y_midd[i][opt.size_x_2 - 1] + /*DD_1(vec_2, i, opt)*/
                (y_midd[i - 1][opt.size_x_2 - 1] - 2 * y_midd[i][opt.size_x_2 - 1] + y_midd[i + 1][opt.size_x_2 - 1]) / (opt.h_x_1 * opt.h_x_1) + phi[i][opt.size_x_2 - 1]
                    + 2 * data.up_bound_func(data.x_1_left + i * opt.h_x_1) / opt.h_x_2);

            right_tridiag_run(y_current[i], A_2, B_2, C_2, D_2);
        }

        if (data.bound[1] == 0) {
            for (size_t j = 0; j < opt.size_x_2; ++j) {
                y_current[opt.size_x_1 - 1][j] = data.right_bound_func(data.x_2_left + j * opt.h_x_2);
            }
        } else {
            //for (size_t n = 0; n < opt.size_x_1; ++n) {
            //    vec_2[n] = y_midd[n][0];
            //}
            D_2[0] = -(1 - data.bound[2]) * data.down_bound_func(data.x_1_right)
                + data.bound[2] * (2 / opt.tau * y_midd[opt.size_x_1 - 1][0] + /*DD_1(vec_2, opt.size_x_1 - 2, opt)*/
                    2 * (y_midd[opt.size_x_1 - 2][0] - y_midd[opt.size_x_1 - 1][0]) / (opt.h_x_1 * opt.h_x_1) + phi[opt.size_x_1 - 1][0]
                    + 2 * data.down_bound_func(data.x_1_right) / opt.h_x_2 + 2 * data.right_bound_func(data.x_2_left) / opt.h_x_1);
            for (size_t j = 1; j < opt.size_x_2 - 1; ++j) {
                //for (size_t n = 0; n < opt.size_x_1; ++n) {
                //    vec_2[n] = y_midd[n][j];
                //}
                D_2[j] = 2 / opt.tau * y_midd[opt.size_x_1 - 1][j] + /*DD_1(vec_2, opt.size_x_1 - 2, opt)*/
                    2 * (y_midd[opt.size_x_1 - 2][j] - y_midd[opt.size_x_1 - 1][j]) / (opt.h_x_1 * opt.h_x_1) + phi[opt.size_x_1 - 1][j]
                    + 2 * data.right_bound_func(data.x_2_left + j * opt.h_x_2) / opt.h_x_1;
            }
            //for (size_t n = 0; n < opt.size_x_1; ++n) {
            //    vec_2[n] = y_midd[n][opt.size_x_2 - 1];
            //}
            D_2[opt.size_x_2 - 1] = -(1 - data.bound[3]) * data.up_bound_func(data.x_1_right)
                + data.bound[3] * (2 / opt.tau * y_midd[opt.size_x_1 - 1][opt.size_x_2 - 1]
                    + /*DD_1(vec_2, opt.size_x_1 - 2, opt)*/ 2 * (y_midd[opt.size_x_1 - 2][opt.size_x_2 - 1] - y_midd[opt.size_x_1 - 1][opt.size_x_2 - 1])
                    / (opt.h_x_1 * opt.h_x_1) + phi[opt.size_x_1 - 1][opt.size_x_2 - 1]
                    + 2 * data.up_bound_func(data.x_1_right) / opt.h_x_2 + 2 * data.right_bound_func(data.x_2_right) / opt.h_x_1);
            //if (k == 1) {
            //    for (auto i = 0; i < opt.size_x_2; ++i) {
            //        std::cout << D_2[i] << "\n";
            //    }
            //}
            right_tridiag_run(y_current[opt.size_x_1 - 1], A_2, B_2, C_2, D_2);
            //if (k == 1) {
            //    std::cout << y_current[opt.size_x_2 - 1];
            //}
        }
        //if (k == 1) {
        //    std::cout << " qqq = " << k << "\n";
        //    for (size_t i = 0; i < opt.size_x_1; ++i) {
        //        std::cout << " i = " << i << "\n" << y_current[i];
        //    }
        //}

        //for (size_t i = 0; i < opt.size_x_1; ++i) {
        //    y_current[i][0] = data.down_bound_func(data.x_1_left + i * opt.h_x_1);
        //    y_current[i][opt.size_x_2 - 1] = data.up_bound_func(data.x_1_left + i * opt.h_x_1);
        //}

//        err = 0.0;
//        for (auto i = 0; i < opt.size_x_1; ++i) {
//            for (auto j = 0; j < opt.size_x_2; ++j) {
//                double diff = std::abs(y_current[i][j] - y_prev[i][j]);
//                if (diff > err) {
//                    err = diff;
//                }
//            }
//        }

        err = 0.0;
        for (auto i = 0; i < opt.size_x_1; ++i) {
            for (auto j = 0; j < opt.size_x_2; ++j) {
                double diff = std::abs(y_current[i][j] - func_anal(data.x_1_left + i * opt.h_x_1,data.x_2_left + j * opt.h_x_2));
                if (diff > err) {
                    err = diff;
                }
            }
        }
//        err_mean += err;

//        clear_file(file_name + "_sol_" + std::to_string(k) + ".dat");
//        export_point(file_name + "_sol_" + std::to_string(k) + ".dat", y_current);

//        clear_file(file_name + "_sol_" + std::to_string(iter * opt.tau) + ".dat");
//        export_point(file_name + "_sol_" + std::to_string(iter * opt.tau) + ".dat", y_current);

        for (size_t i = 0; i < opt.size_x_1; ++i) {
            y_prev[i] = y_current[i];
        }
    } //while (err > eps);
    std::cout << opt.size_t - 1 << "    " << err << std::endl;
//    clear_file(file_name + "_tau_" + std::to_string(opt.tau) + "_sol_" + std::to_string(iter * opt.tau) + ".dat");
//    export_point(file_name + "_tau_" + std::to_string(opt.tau) + "_sol_" + std::to_string(iter * opt.tau) + ".dat", y_current);
    //std::cout << "eps = " << eps * opt.tau << "\niter = " << iter << "\n err = " << std::fixed << std::setprecision(20) << err << "\n";
    //std::cout << "eps = " << eps << "\niter = " << iter << "\n tau = " << opt.tau << "\n time = " << iter * opt.tau << "\n";
    //std::cout << std::fixed << std::setprecision(5) << iter * opt.tau << ", ";

    std::ofstream file("/Users/artem/Desktop/study/вычи/4/export/test_3/test_3_co.dat", std::ios::app);
    file << std::fixed << std::setprecision(10) << err << iter << std::endl;
    file.close();
    file.clear();
}
