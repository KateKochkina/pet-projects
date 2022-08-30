#pragma once

#include <functional>
#include <vector>

struct Point {
    explicit Point (double x);
    Point (double x, double y);

    double x;
    double y;
};

struct Interval {
    explicit Interval(Point left_border, Point right_border);
    Interval(double left_border, double right_border);

    Point left_border;
    Point right_border;
};

struct Meta {
    Meta(double value, size_t iteration_counter);
    Meta(std::vector<double> &x, size_t iteration_counter);

    std::vector<double> x;
    size_t iteration_counter;
};

std::vector<Interval> get_intervals_with_roots(const std::function<double (double)> &f, const Interval &interval,
        size_t number_of_mesh_points = 11);

std::vector<double> uniform_mesh(double left_border, double right_border, size_t number_of_mesh_points);

std::vector<double> get_function_values(const std::function<double (double)> &f, std::vector<double> mesh);

double get_derivative(const std::function<double (double)> &f, double x, double precision = 1e-6);

std::vector<std::vector<double>> get_jacobian(const std::function<std::vector<double>(std::vector<double>)> &f,
    std::vector<double> x, double precision = 1e-6);

std::vector<std::vector<double>> inverse_jacobian(const std::vector<std::vector<double>> &jacobian);