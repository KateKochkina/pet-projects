#include <iostream>
#include <cmath>

#include "utils.h"

std::vector<double> get_function_values(const std::function<double (double)> &f, std::vector<double> mesh) {
    std::vector<double> function_values(mesh.size());

    for (size_t i = 0; i < function_values.size(); ++i) {
        function_values[i] = f(mesh[i]);
    }

    return function_values;
}

std::vector<double> uniform_mesh(double left_border, double right_border, size_t number_of_mesh_points) {
    std::vector<double> mesh(number_of_mesh_points + 1);
    double increment = (right_border - left_border) / number_of_mesh_points;

    for (size_t i = 0; i <= number_of_mesh_points; ++i) {
        mesh[i] = left_border + increment * i;
    }

    return mesh;
}

std::vector<Interval> get_intervals_with_roots(const std::function<double (double)> &f, const Interval &interval,
        size_t number_of_mesh_points) {
    std::vector<Interval> intervals_with_roots;

    auto mesh = uniform_mesh(interval.left_border.x, interval.right_border.x, number_of_mesh_points);
    auto function_values = get_function_values(f, mesh);

    for (size_t i = 0; i < mesh.size() - 1; ++i) {
        if (function_values[i] * function_values[i + 1] <= 0) {
            intervals_with_roots.emplace_back(Point(mesh[i], function_values[i]),
                                              Point(mesh[i + 1], function_values[i + 1]));
        }
    }

    return intervals_with_roots;
}

double get_derivative(const std::function<double (double)> &f, double x, double precision) {
    return (f(x + precision) - f(x)) / precision;
}

std::vector<std::vector<double>> get_jacobian(const std::function<std::vector<double>(std::vector<double>)> &f,
    std::vector<double> x, double precision) {
    std::vector<std::vector<double>> jacobian(2);

    jacobian[0] = {(f({x[0] + precision, x[1]})[0] - f({x[0], x[1]})[0]) / precision,
                   (f({x[0], x[1] + precision})[0] - f({x[0], x[1]})[0]) / precision};
    jacobian[1] = {(f({x[0] + precision, x[1]})[1] - f({x[0], x[1]})[1]) / precision,
                   (f({x[0], x[1] + precision})[1] - f({x[0], x[1]})[1]) / precision};

    return jacobian;
}

std::vector<std::vector<double>> inverse_jacobian(const std::vector<std::vector<double>> &jacobian) {
    std::vector<std::vector<double>> matrix(2);
    
    double denom = -jacobian[0][1] * jacobian[1][0] + jacobian[0][0] * jacobian[1][1];
    matrix[0] = {jacobian[1][1] / denom,
                 -(jacobian[0][1] /denom)};
    matrix[1] = {-(jacobian[1][0] / denom),
                 jacobian[0][0] / denom};

    return matrix;
}

Point::Point(double x) : x(x), y(0) {}

Point::Point(double x, double y) : x(x), y(y) {}

Interval::Interval(Point lb, Point rb) : left_border(lb), right_border(rb) {}

Interval::Interval(double lb, double rb) : left_border(lb), right_border(rb) {}

Meta::Meta(double value, size_t iteration_counter)
        : iteration_counter(iteration_counter) {
    x = { value };
}

Meta::Meta(std::vector<double> &x, size_t iteration_counter) : x(x), iteration_counter(iteration_counter) {}