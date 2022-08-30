#include <iostream>
#include <string>
#include "methods.h"

void function_copy_pasta(const std::string &file_name, Problem_Data data, Options options, TVector<double>& y_prev,
            Bound_cond_func &bound_cond_func_left, Bound_cond_func &bound_cond_func_right) {
    clear_file(file_name + "_mesh.dat");
    export_mesh(file_name + "_mesh.dat", options);

    //начальные условия - решение на нулевом временном слое
    for (auto i = 0; i < options.size_x; ++i) {
        y_prev[i] = data.init_cond_func(i * options.h, data);
    }
    clear_file(file_name + "_res.dat");
    export_point(file_name + "_res.dat", y_prev);

    //определение функции гран. условий
    bound_cond_func_left = bound_cond_left_1;
    if (data.left_bound_cond_type == BOUND_COND_OF_THE_SECOND_TYPE) {
        bound_cond_func_left = bound_cond_left_2;
    }
    bound_cond_func_right = bound_cond_right_1;
    if (data.right_bound_cond_type == BOUND_COND_OF_THE_SECOND_TYPE) {
        bound_cond_func_right = bound_cond_right_2;
    }
}

void linear_heat_eq(const std::string &file_name, Problem_Data data, Options options, TVector<double>& y_prev,
        TVector<double>& y_next, TVector<double>& A, TVector<double>& B,
        TVector<double>& C, TVector<double>& D, TVector<double>& a) {

    Bound_cond_func bound_cond_func_left;
    Bound_cond_func bound_cond_func_right;

    function_copy_pasta(file_name + "_linear", data, options, y_prev, bound_cond_func_left, bound_cond_func_right);

    std::string f_name = file_name + "_linear_res.dat";
    //a[0] в формулах нигде не фигурирует, но для порядка (чтобы не сбить нумерацию) он есть
    for (auto i = 1; i < options.size_x; ++i) {
        a[i] = 0.5 * (data.K(i * options.h) + data.K((i - 1) * options.h)); // options.K((i - 0.5) * h)
    }
    //заполнение трех диагоналей матрицы
    for (auto i = 1; i < options.size_x - 1; ++i) {
        A[i] = options.tau * options.sigma * a[i];
        B[i] = data.c * data.rho * options.h * options.h + options.tau * options.sigma * (a[i] + a[i + 1]);
        C[i] = options.tau * options.sigma * a[i + 1];
    }

    for (auto j = 1; j < options.size_t; ++j) {
        bound_cond_func_left(A, B, C, D, a, data, options, j * options.tau, y_prev);
        bound_cond_func_right(A, B, C, D, a, data, options, j * options.tau, y_prev);
        for (auto i = 1; i < options.size_x - 1; ++i) {
            D[i] = data.c * data.rho * options.h * options.h * y_prev[i] + (1 - options.sigma) *
                options.tau * (a[i + 1] * (y_prev[i + 1] - y_prev[i]) - a[i] * (y_prev[i] - y_prev[i - 1]));
        }
        right_tridiag_run(y_next, A, B, C, D);
        export_point(f_name, y_next);
        y_prev = y_next;
    }
}

void quasilinear_heat_eq(const std::string &file_name, Problem_Data data, Options options,
        TVector<double>& y_prev, TVector<double>& y_next,
        TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D, TVector<double>& a) {

    Bound_cond_func bound_cond_func_left;
    Bound_cond_func bound_cond_func_right;

    function_copy_pasta(file_name + "_quasilinear", data, options, y_prev, bound_cond_func_left, bound_cond_func_right);

    std::string f_name = file_name + "_quasilinear_res.dat";


    for (auto j = 1; j < options.size_t; ++j) {
        for (auto i = 1; i < options.size_x; ++i) {
            a[i] = 0.5 * (data.K(y_prev[i]) + data.K(y_prev[i - 1]));
        }
        bound_cond_func_left(A, B, C, D, a, data, options, j * options.tau, y_prev);
        bound_cond_func_right(A, B, C, D, a, data, options, j * options.tau, y_prev);
        for (auto i = 1; i < options.size_x - 1; ++i) {
            A[i] = options.tau * options.sigma * a[i];
            B[i] = data.c * data.rho * options.h * options.h + options.tau * options.sigma * (a[i] + a[i + 1]);
            C[i] = options.tau * options.sigma * a[i + 1];
            D[i] = data.c * data.rho * options.h * options.h * y_prev[i] + (1 - options.sigma) * options.tau *
                (a[i + 1] * (y_prev[i + 1] - y_prev[i]) - a[i] * (y_prev[i] - y_prev[i - 1]));
        }
        right_tridiag_run(y_next, A, B, C, D);
        export_point(f_name, y_next);
        y_prev = y_next;
    }
}

void quasilinear_heat_eq_non_linear(const std::string &file_name, Problem_Data data, Options options,
        TVector<double>& y_prev, TVector<double>& y_next,
        TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D, TVector<double>& a) {

    Bound_cond_func bound_cond_func_left;
    Bound_cond_func bound_cond_func_right;

    function_copy_pasta(file_name + "_quasilinear_non_linear", data, options, y_prev, bound_cond_func_left, bound_cond_func_right);

    std::string f_name = file_name + "_quasilinear_non_linear_res.dat";


    for (auto j = 1; j < options.size_t; ++j) {
        y_next = y_prev;
        for (auto s = 0; s < 3; ++s) {
            for (auto i = 1; i < options.size_x; ++i) {
                a[i] = 0.5 * (data.K(y_next[i]) + data.K(y_next[i - 1]));
            }
            bound_cond_func_left(A, B, C, D, a, data, options, j * options.tau, y_prev);
            bound_cond_func_right(A, B, C, D, a, data, options, j * options.tau, y_prev);
            for (auto i = 1; i < options.size_x - 1; ++i) {
                A[i] = options.tau * options.sigma * a[i];
                B[i] = data.c * data.rho * options.h * options.h + options.tau * options.sigma * (a[i] + a[i + 1]);
                C[i] = options.tau * options.sigma * a[i + 1];
                D[i] = data.c * data.rho * options.h * options.h * y_prev[i];
            }
            right_tridiag_run(y_next, A, B, C, D);
        }
        export_point(f_name, y_next);
        y_prev = y_next;
    }
}

void quasilinear_heat_eq_non_linear_iter(const std::string &file_name, Problem_Data data, Options options,
    TVector<double>& y_prev, TVector<double>& y_next,
    TVector<double>& A, TVector<double>& B, TVector<double>& C, TVector<double>& D, TVector<double>& a, double tol) {

    Bound_cond_func bound_cond_func_left;
    Bound_cond_func bound_cond_func_right;

    //начальные условия - решение на нулевом временном слое
    for (auto i = 0; i < options.size_x; ++i) {
        y_prev[i] = data.init_cond_func(i * options.h, data);
    }

    //определение функции гран. условий
    bound_cond_func_left = bound_cond_left_1;
    if (data.left_bound_cond_type == BOUND_COND_OF_THE_SECOND_TYPE) {
        bound_cond_func_left = bound_cond_left_2;
    }
    bound_cond_func_right = bound_cond_right_1;
    if (data.right_bound_cond_type == BOUND_COND_OF_THE_SECOND_TYPE) {
        bound_cond_func_right = bound_cond_right_2;
    }

    std::string f_name = file_name + "_" + std::to_string(tol);
    clear_file(f_name + "_iter.dat");

    clear_file(f_name + "_mesh.dat");
    export_mesh(f_name + "_mesh.dat", options);

    std::ofstream file1("/Users/artem/Desktop/study/вычи/2/export/" + f_name + "_iter.dat", std::ios::app);


    TVector<double> q(y_prev.get_size());
    for (auto j = 1; j < options.size_t; ++j) {
        y_next = y_prev;
        size_t num_iter = 0;
        do {
            q = y_next;
            for (auto i = 1; i < options.size_x; ++i) {
                a[i] = 0.5 * (data.K(y_next[i]) + data.K(y_next[i - 1]));
            }
            bound_cond_func_left(A, B, C, D, a, data, options, j * options.tau, y_prev);
            bound_cond_func_right(A, B, C, D, a, data, options, j * options.tau, y_prev);
            for (auto i = 1; i < options.size_x - 1; ++i) {
                A[i] = options.tau * options.sigma * a[i];
                B[i] = data.c * data.rho * options.h * options.h + options.tau * options.sigma * (a[i] + a[i + 1]);
                C[i] = options.tau * options.sigma * a[i + 1];
                D[i] = data.c * data.rho * options.h * options.h * y_prev[i];
            }
            right_tridiag_run(y_next, A, B, C, D);
            num_iter++;
        } while ((y_next - q).norm(SPHERICAL_NORM) > tol);
        file1 << std::fixed << std::setprecision(10) << num_iter << "\n";
        y_prev = y_next;
    }
    file1.close();
    file1.clear();
}
