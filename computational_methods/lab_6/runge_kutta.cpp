#include <cmath>
#include <iostream>

#include "runge_kutta.h"

void runge_kutta_2_const_step(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                              TVector<double> starting_point, Interval time_interval) {

    size_t size = starting_point.get_size();

    TVector<double> desired_point(size);
    TVector<double> k_1(size);
    TVector<double> k_2(size);

    double sigma_1 = 0.0;  //0.5
    double sigma_2 = 1.0;  //0.5
    double a_2 = 0.5;      //1.0
    double b_21 = 0.5;     //1.0
    double step = time_interval.left_border;

    std::string filename = file_name + "_r_k_2_const_step.dat";
    clear_file(filename);
    export_point(filename, step, starting_point);
    while (step < time_interval.right_border) {
        k_1 = func(step, starting_point);
        k_2 = func(step + a_2 * tau, starting_point + tau * b_21 * k_1);
        desired_point = starting_point + tau * (sigma_1 * k_1 + sigma_2 * k_2);
        step += tau;
        export_point(filename, step, desired_point);
        starting_point = desired_point;
    }
}

void runge_kutta_2_variable_step(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                                 TVector<double> starting_point, Interval time_interval, Options options) {

    size_t p = 2;
    double tol = 1e-8;

    size_t size = starting_point.get_size();

    TVector<double> intermediate_point(size);
    TVector<double> desired_point(size);
    TVector<double> service_point(size);

    TVector<double> k_1(size);
    TVector<double> k_1_service(size);
    TVector<double> k_2(size);

    double sigma_1 = 0.0;  //0.5
    double sigma_2 = 1.0;  //0.5
    double a_2 = 0.5;      //1.0
    double b_21 = 0.5;     //1.0
    double step = time_interval.left_border;

    double tau_service;
    double error;

    double tau_new;
    bool was_tau_dem = true;

    std::string filename = file_name + "_r_k_2_var_step.dat";
    clear_file(filename);
    export_point(filename, step, starting_point);

    std::string filename_tau = file_name + "_r_k_2_var_step_tau.dat";
    clear_file(filename_tau);
    while (step < time_interval.right_border) {
        if (options.step_recount == false) {
            k_1 = func(step, starting_point);
            k_1_service = k_1;
        }

        k_2 = func(step + a_2 * tau, starting_point + tau * b_21 * k_1_service);
        intermediate_point = starting_point + tau * (sigma_1 * k_1_service + sigma_2 * k_2);

        k_1 = func(step, intermediate_point);
        k_2 = func(step + a_2 * tau, intermediate_point + tau * b_21 * k_1);
        desired_point = intermediate_point + tau * (sigma_1 * k_1 + sigma_2 * k_2);

        tau_service = 2 * tau;
        k_2 = func(step + a_2 * tau_service, starting_point + tau_service * b_21 * k_1_service);
        service_point = starting_point + tau_service * (sigma_1 * k_1_service + sigma_2 * k_2);

        error = ((desired_point - service_point) / (1 - pow(2, p))).norm(SPHERICAL_NORM);
        tau_new = tau * fmin(options.fac_max, fmax(options.fac_min, options.fac * pow(tol / error, 1.0 / (p + 1))));
        ;
        if (error > tol || error < (tol / 100)) {
            tau = tau_new;
            options.step_recount = true;
            was_tau_dem = true;
            continue;
        }

        options.step_recount = false;

        step += tau;
        export_tau(filename_tau, step, tau, error);
        export_point(filename, step, intermediate_point);
        step += tau;
        export_tau(filename_tau, step, tau, error);
        export_point(filename, step, desired_point);
        starting_point = desired_point;

        if (tau_new < tau) {
            tau = tau_new;
            was_tau_dem = true;
            continue;
        } else if (!was_tau_dem) {
            tau = tau_new;
        }
        was_tau_dem = false;
    }
}

void runge_kutta_4_const_step(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                              TVector<double> starting_point, Interval time_interval) {

    size_t size = starting_point.get_size();

    TVector<double> desired_point(size);

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
    double b_32 = 0.5;
    double b_43 = 1.0;
    double step = time_interval.left_border;

    std::string filename = file_name + "_r_k_4_const_step.dat";
    clear_file(filename);
    export_point(filename, step, starting_point);
    while (step < time_interval.right_border) {
        k_1 = func(step, starting_point);
        k_2 = func(step + a_2 * tau, starting_point + tau * b_21 * k_1);
        k_3 = func(step + a_3 * tau, starting_point + tau * b_32 * k_2);
        k_4 = func(step + a_4 * tau, starting_point + tau * b_43 * k_3);

        desired_point = starting_point + tau * (sigma_1 * k_1 + sigma_2 * k_2 + sigma_3 * k_3 + sigma_4 * k_4);
        step += tau;
        export_point(filename, step, desired_point);
        starting_point = desired_point;
    }
}

void runge_kutta_4_variable_step(const std::string& file_name, TVector<double> func(double, TVector<double>), double tau,
                                 TVector<double> starting_point, Interval time_interval, Options options) {

    size_t p = 4;
    double tol = 1e-12;

    size_t size = starting_point.get_size();

    TVector<double> intermediate_point(size);
    TVector<double> desired_point(size);
    TVector<double> service_point(size);

    TVector<double> k_1(size);
    TVector<double> k_1_service(size);
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
    double b_32 = 0.5;
    double b_43 = 1.0;
    double step = time_interval.left_border;

    double tau_service;
    double error;

    double tau_new;
    bool was_tau_dem = true;

    std::string filename = file_name + "_r_k_4_var_step.dat";
    clear_file(filename);
    export_point(filename, step, starting_point);

    std::string filename_tau = file_name + "_r_k_4_var_step_tau.dat";
    clear_file(filename_tau);
    while (step < time_interval.right_border) {
        if (!options.step_recount) {
            k_1 = func(step, starting_point);
            k_1_service = k_1;
        }

        k_2 = func(step + a_2 * tau, starting_point + tau * b_21 * k_1_service);
        k_3 = func(step + a_3 * tau, starting_point + tau * b_32 * k_2);
        k_4 = func(step + a_4 * tau, starting_point + tau * b_43 * k_3);
        intermediate_point = starting_point + tau * (sigma_1 * k_1_service + sigma_2 * k_2 + sigma_3 * k_3 + sigma_4 * k_4);

        k_1 = func(step, intermediate_point);
        k_2 = func(step + a_2 * tau, intermediate_point + tau * b_21 * k_1);
        k_3 = func(step + a_3 * tau, intermediate_point + tau * b_32 * k_2);
        k_4 = func(step + a_4 * tau, intermediate_point + tau * b_43 * k_3);
        desired_point = intermediate_point + tau * (sigma_1 * k_1 + sigma_2 * k_2 + sigma_3 * k_3 + sigma_4 * k_4);

        tau_service = 2 * tau;
        k_2 = func(step + a_2 * tau_service, starting_point + tau_service * b_21 * k_1_service);
        k_3 = func(step + a_3 * tau_service, starting_point + tau_service * b_32 * k_2);
        k_4 = func(step + a_4 * tau_service, starting_point + tau_service * b_43 * k_3);
        service_point = starting_point + tau_service * (sigma_1 * k_1_service + sigma_2 * k_2 + sigma_3 * k_3 + sigma_4 * k_4);

        error = fabs((service_point[0] - desired_point[0]) / (pow(2.0, p) - 1)); //.norm(SPHERICAL_NORM);
        tau_new = tau * fmin(options.fac_max, fmax(options.fac_min, options.fac * pow(tol / error, 1.0 / (p + 1.0))));
        if (error > tol || error < tol / 100) {
            tau = tau_new;
            was_tau_dem = true;
            options.step_recount = true;
            continue;
        }
        options.step_recount = false;

        step += tau;
        export_tau(filename_tau, step, tau, error);
        export_point(filename, step, intermediate_point);
        step += tau;
        export_tau(filename_tau, step, tau, error);
        export_point(filename, step, desired_point);
        starting_point = desired_point;

        if (tau_new < tau) {
            tau = tau_new;
            was_tau_dem = true;
        } else if (was_tau_dem) {
            was_tau_dem = false;
        } else {
            tau = tau_new;
        }
    }
}