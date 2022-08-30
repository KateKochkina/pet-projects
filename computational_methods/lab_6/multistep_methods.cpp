#include <iostream>

#include "multistep_methods.h"

TVector<TVector<double>> runge_kutta_4_step_3(TVector<double> func(double, TVector<double>), double tau,
                                              TVector<double> starting_point, Interval time_interval) {

    size_t size = starting_point.get_size();

    TVector<TVector<double>> desired_point(3);

    TVector<double> k_1(size);
    TVector<double> k_2(size);
    TVector<double> k_3(size);
    TVector<double> k_4(size);

    double sigma_1 = 1.0 / 6.0;
    double sigma_2 = 1.0 / 3.0;
    double sigma_3 = 1.0 / 3.0;
    double sigma_4 = 1.0 / 6.0;
    double a_2 = 0.5;
    double a_3 = 0.5;
    double a_4 = 1.0;
    double b_21 = 0.5;
    double b_31 = 0.0;
    double b_32 = 0.5;
    double b_41 = 0.0;
    double b_42 = 0.0;
    double b_43 = 1.0;
    double step = time_interval.left_border;

    for (size_t i = 0; i < 3; ++i) {
        k_1 = func(step, starting_point);
        k_2 = func(step + a_2 * tau, starting_point + tau * b_21 * k_1);
        k_3 = func(step + a_3 * tau, starting_point + tau * (b_31 * k_1 + b_32 * k_2));
        k_4 = func(step + a_4 * tau, starting_point + tau * (b_41 * k_1 + b_42 * k_2 + b_43 * k_3));

        desired_point[i] = starting_point + tau * (sigma_1 * k_1 + sigma_2 * k_2 + sigma_3 * k_3 + sigma_4 * k_4);

        step += tau;
        starting_point = desired_point[i];
    }

    return desired_point;
}

void adams_bashforth(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                     TVector<double> starting_point, Interval time_interval) {

    size_t size = starting_point.get_size();

    double b_1 = 55.0;
    double b_2 = -59.0;
    double b_3 = 37.0;
    double b_4 = -9.0;

    TVector<TVector<double>> points = runge_kutta_4_step_3(func, tau, starting_point, time_interval);

    TVector<double> point_n_3(starting_point);
    TVector<double> point_n_2(points[0]);
    TVector<double> point_n_1(points[1]);
    TVector<double> point_n(points[2]);
    TVector<double> desired_point(size);

    double step = time_interval.left_border + 3 * tau;

    TVector<double> f_n_3(func(step - 3 * tau, point_n_3));
    TVector<double> f_n_2(func(step - 2 * tau, point_n_2));
    TVector<double> f_n_1(func(step - tau, point_n_1));
    TVector<double> f_n(size);

    std::string filename = file_name + "_ad_bash.dat";
    clear_file(filename);
    export_point(filename, step, starting_point);
    export_point(filename, step, point_n_2);
    export_point(filename, step, point_n_1);
    export_point(filename, step, point_n);
    while (step < time_interval.right_border) {
        f_n = func(step, point_n);

        desired_point = point_n + tau / 24 * (b_1 * f_n + b_2 * f_n_1 + b_3 * f_n_2 + b_4 * f_n_3);
        step += tau;
        export_point(filename, step, desired_point);

        point_n_3 = point_n_2;
        point_n_2 = point_n_1;
        point_n_1 = point_n;
        point_n = desired_point;

        f_n_3 = f_n_2;
        f_n_2 = f_n_1;
        f_n_1 = f_n;
    }
}

void predictor_corrector(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                         TVector<double> starting_point, Interval time_interval) {

    size_t size = starting_point.get_size();

    double b_1_pred = 55.0;
    double b_2_pred = -59.0;
    double b_3_pred = 37.0;
    double b_4_pred = -9.0;

    double b_0_corr = 9.0;
    double b_1_corr = 19.0;
    double b_2_corr = -5.0;
    double b_3_corr = 1.0;

    TVector<TVector<double>> points = runge_kutta_4_step_3(func, tau, starting_point, time_interval);

    TVector<double> point_n_3(starting_point);
    TVector<double> point_n_2(points[0]);
    TVector<double> point_n_1(points[1]);
    TVector<double> point_n(points[2]);

    TVector<double> desired_point_pred(size);
    TVector<double> desired_point(size);

    double step = time_interval.left_border + 3 * tau;

    TVector<double> f_n_3(func(step - 3 * tau, point_n_3));
    TVector<double> f_n_2(func(step - 2 * tau, point_n_2));
    TVector<double> f_n_1(func(step - tau, point_n_1));
    TVector<double> f_n(size);

    std::string filename = file_name + "_pred_corr.dat";
    clear_file(filename);
    export_point(filename, step, starting_point);
    export_point(filename, step, point_n_2);
    export_point(filename, step, point_n_1);
    export_point(filename, step, point_n);
    while (step < time_interval.right_border) {
        f_n = func(step, point_n);

        desired_point_pred = point_n + tau / 24 * (b_1_pred * f_n + b_2_pred * f_n_1 + b_3_pred * f_n_2 + b_4_pred * f_n_3);
        desired_point = point_n + tau / 24 * (b_0_corr * func(step + tau, desired_point_pred) + b_1_corr * f_n + b_2_corr * f_n_1 + b_3_corr * f_n_2);
        step += tau;
        export_point(filename, step, desired_point);

        point_n_3 = point_n_2;
        point_n_2 = point_n_1;
        point_n_1 = point_n;
        point_n = desired_point;

        f_n_3 = f_n_2;
        f_n_2 = f_n_1;
        f_n_1 = f_n;
    }
}