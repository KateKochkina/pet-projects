#include <fstream>

#include "utils.h"

Interval::Interval(double l_b, double r_b) : left_border(l_b), right_border(r_b) {}

Options::Options() : fac_min(0.1), fac_max(4.0), fac(0.9), step_recount(false) {}
Options::Options(double f_min, double f_max, double f, bool step_rec)
    : fac_min(f_min), fac_max(f_max), fac(f), step_recount(step_rec) {}

void export_point(const std::string& file_name, const double step, const TVector<double> point) {
    std::ofstream file("export/" + file_name, std::ios::app);
    file << std::fixed << std::setprecision(6) << step << "    ";
    for (size_t i = 0; i < point.get_size(); ++i) {
        file << std::fixed << std::setprecision(15) << point[i] << "    ";
    }
    file << std::endl;
    file.close();
    file.clear();
}

void export_tau(const std::string& file_name, const double step, double tau, double error) {
    std::ofstream file("export/" + file_name, std::ios::app);
    file << std::fixed << std::setprecision(10) << step << "    " << tau << "    " << std::setprecision(16) << error << "\n";
    file << std::endl;
    file.close();
    file.clear();
}

void clear_file(const std::string& file_name) {
    std::ofstream file("export/" + file_name, std::ios::trunc);
    file.close();
    file.clear();
}

TVector<double> function_euler_implict(double t, TVector<double> x, TVector<double> starting_point,
                                       TVector<double> func(double, TVector<double>), double tau) {
    return x - starting_point - tau * func(t + tau, x);
}

TVector<double> function_symmetric_scheme(double t, TVector<double> x, TVector<double> starting_point,
                                          TVector<double> func(double, TVector<double>), double tau) {
    return x - starting_point - tau / 2 * (func(t + tau, x) + func(t, starting_point));
}

TVector<TVector<double>> get_jacobian(TVector<double> func(double, TVector<double>),
                                      TVector<double> func_meth(double, TVector<double>, TVector<double>, TVector<double> (*)(double, TVector<double>), double),
                                      TVector<double>& x, TVector<double>& starting_point, double t, double tau, double precision) {

    double size = starting_point.get_size();
    TVector<TVector<double>> jacobian(size);
    for (size_t i = 0; i < size; ++i) {
        jacobian[i].resize(size);
    }

    TVector<double> delta_x = x;
    for (size_t i = 0; i < size; ++i) {
        delta_x[i] += precision;
        for (size_t j = 0; j < size; ++j) {
            jacobian[j][i] = (func_meth(t, delta_x, starting_point, func, tau)[j] - func_meth(t, x, starting_point, func, tau)[j]) / precision;
        }
        delta_x[i] -= precision;
    }

    //TVector<double> delta_x_0 = { x[0] + precision, x[1] };
    //TVector<double> delta_x_1 = { x[0], x[1] + precision };
    //
    //jacobian[0] = { (func_meth(t, delta_x_0, starting_point, func, tau)[0] - func_meth(t, x, starting_point, func, tau)[0]) / precision,
    //               (func_meth(t, delta_x_1, starting_point, func, tau)[0] - func_meth(t, x, starting_point, func, tau)[0]) / precision };
    //jacobian[1] = { (func_meth(t, delta_x_0, starting_point, func, tau)[1] - func_meth(t, x, starting_point, func, tau)[1]) / precision,
    //               (func_meth(t, delta_x_1, starting_point, func, tau)[1] - func_meth(t, x, starting_point, func, tau)[1]) / precision };

    return jacobian;
}

TVector<TVector<double>> inverse_jacobian(TVector<TVector<double>>& jacobian) {
    TVector<TVector<double>> matrix(2);

    double denom = -jacobian[0][1] * jacobian[1][0] + jacobian[0][0] * jacobian[1][1];
    matrix[0] = {jacobian[1][1] / denom,
                 -(jacobian[0][1] / denom)};
    matrix[1] = {-(jacobian[1][0] / denom),
                 jacobian[0][0] / denom};

    return matrix;
}

TVector<double> newton_method(TVector<double> func(double, TVector<double>),
                              TVector<double> func_meth(double, TVector<double>, TVector<double>, TVector<double> (*)(double, TVector<double>), double),
                              TVector<double>& starting_point, double t, double tau, double precision) {

    TVector<double> point_next(starting_point);
    TVector<double> point_prev(starting_point.get_size());
    do {
        point_prev = point_next;
        TVector<TVector<double>> jacobian = get_jacobian(func, func_meth, point_prev, starting_point, t, tau);
        point_next = point_prev - inverse_jacobian(jacobian) * func_meth(t, point_prev, starting_point, func, tau);
    } while ((point_next - point_prev).norm(SPHERICAL_NORM) >= precision);
    return point_next;
}