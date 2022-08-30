#include <iostream>
#include <string>

#include "euler.h"
#include "multistep_methods.h"
#include "runge_kutta.h"
#include "tests.hpp"
#include "utils.h"

int main() {

    std::string filename =
//            "test";
//            "f_analysis";
            "f_pendulum";
//            "f_nonlinear_oscillatory_circuit";

    /* Const step */
//    Interval t_int = {0.0, 100.0};
//    double tau = 0.01;
//    TVector<double> starting_point = {1.0, 0.0};
//    for (size_t i = 0; i < 1; ++i) {
//        runge_kutta_2_const_step(filename + "_" + std::to_string(tau / pow(2.0, i)), function_pendulum,
//                       tau / pow(2.0, i), starting_point, t_int);
//        std::cout << "step\n";
//    }
//    std::cout << "\n";

    /* Variable step */
    Interval t_int = {0.0, 2.0};
    double tau = 0.0001;
    TVector<double> starting_point = {1.0, 0.0};
    Options options;
    for (size_t i = 0; i < 1; ++i) {
        runge_kutta_4_variable_step(filename + "_" + std::to_string(tau / pow(2.0, i)), function_pendulum, tau / pow(2.0, i), starting_point, t_int, options);
        std::cout << "step\n";
    }
    std::cout << "\n";

//    /* Nonlinear oscillatory circuit */
//    Interval t_int = {0.0, 1.0};
//    double tau = 0.001;
//    size_t points_num = 1;  // кол-во точек на окружности с центром в особой точке
//    double R = 0.5;         // радиус окружности
//    TVector<double> circle_point(2);
//    TVector<double> special_point = {0.0, 0.0};
//    //runge_kutta_4_const_step(filename + "_" + std::to_string(0), function_nonlinear_oscillatory_circuit, tau, special_point, t_int);
//    for (auto i = 0; i < points_num; ++i) {
//        circle_point = {special_point[0] + R * cos(i * 2.0 * PI / points_num), special_point[1] + R * sin(i * 2.0 * PI / points_num)};
////        symmetric_scheme(filename + "_" + std::to_string(i), function_nonlinear_oscillatory_circuit, tau, circle_point, t_int);
//        runge_kutta_4_variable_step(filename + "_" + std::to_string(i), function_nonlinear_oscillatory_circuit, tau, circle_point, t_int, options);
//        std::cout << "step\n";
//    }

    std::cout << "\ndone!\n";
    return 0;
}
