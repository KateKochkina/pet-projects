#include "find_root.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

const size_t MAX_ITERATION_NEWTON = 40;

Meta bisection_method(const std::function<double (double)> &f, const Interval &interval, double precision) {
    double double_precision = precision * 2;

    Point answer((interval.left_border.x + interval.right_border.x) / 2);
    size_t iteration_counter = 1;

    Point left_border = interval.left_border;
    Point right_border = interval.right_border;

    while ((right_border.x - left_border.x) >= double_precision) {
        answer.y = f(answer.x);

        if (std::abs(answer.y) < precision) {
            return Meta(answer.x, iteration_counter);
        }

        if (left_border.y * answer.y < 0) {
            right_border = answer;
        } else {
            left_border = answer;
        }

        answer.x = (left_border.x + right_border.x) / 2;
        ++iteration_counter;
    }

    return Meta(answer.x, iteration_counter);
}

Meta newton_method(const std::function<double (double)> &f, const std::function<double (double)> &f_d,
        double prev, const Interval &interval, double precision) {
    double answer = prev - f(prev) / f_d(prev);
    size_t iteration_counter = 1;

    if (answer < interval.left_border.x || answer > interval.right_border.x) {
        std::cout << "*out of range* ";
        return Meta(answer, iteration_counter);
    }

    while (std::abs(answer - prev) >= precision) {
        if (iteration_counter > MAX_ITERATION_NEWTON) {
            std::cout << "*iteration_counter > " << MAX_ITERATION_NEWTON << "* ";
            return Meta(answer, iteration_counter);
        }

        prev = answer;
        answer = prev - f(prev) / f_d(prev);
        ++iteration_counter;

        if (answer < interval.left_border.x || answer > interval.right_border.x) {
            std::cout << "*out of range* ";
            return Meta(answer, iteration_counter);
        }
    }

    return Meta(answer, iteration_counter);
}

Meta newton_method_numerical(const std::function<double (double)> &f,
        double prev, const Interval &interval, double precision) {
    double answer = prev - f(prev) / get_derivative(f, prev, precision);
    size_t iteration_counter = 1;

    if (answer < interval.left_border.x || answer > interval.right_border.x) {
        std::cout << "*out of range* ";
        return Meta(answer, iteration_counter);
    }

    while (std::abs(answer - prev) >= precision) {
        if (iteration_counter > MAX_ITERATION_NEWTON) {
            std::cout << "*iteration_counter > " << MAX_ITERATION_NEWTON << "* ";
            return Meta(answer, iteration_counter);
        }

        prev = answer;
        answer = prev - f(prev) / get_derivative(f, prev, precision);
        ++iteration_counter;

        if (answer < interval.left_border.x || answer > interval.right_border.x) {
            std::cout << "*out of range* ";
            return Meta(answer, iteration_counter);
        }
    }

    return Meta(answer, iteration_counter);
}

//! requires function value at interval points
Meta newton_method_modification(const std::function<double (double)> &f, const std::function<double (double)> &f_d,
        const std::function<double (double)> &f_dd, double prev, const Interval &interval, double precision) {
    double answer = prev - f(prev) / f_d(prev);
    size_t iteration_counter = 1;

    if (answer < interval.left_border.x || answer > interval.right_border.x) {
        answer = get_chord_method_approx(interval);
        ++iteration_counter;
    }

    double chord_border;
    if (f(interval.left_border.x) * f_dd(interval.left_border.x) > 0) {
        chord_border = interval.left_border.x;
    } else {
        chord_border = interval.right_border.x;
    }

    size_t bad_init_approx_counter = 0;

    while (std::abs(answer - prev) >= precision) {
        if (iteration_counter > MAX_ITERATION_NEWTON) {
            ++bad_init_approx_counter;

            if (bad_init_approx_counter == 1) {
                std::cout << "*iteration_counter > " << MAX_ITERATION_NEWTON
                          << ", getting init_approx using the chord method* ";
                answer = get_chord_method_approx(interval);
            } else {
                if (bad_init_approx_counter % 2 == 0) {
                    answer = interval.left_border.x
                             + (interval.right_border.x - interval.left_border.x)
                             / (bad_init_approx_counter);
                } else {
                    answer = interval.right_border.x
                             - (interval.right_border.x - interval.left_border.x)
                             / (bad_init_approx_counter);
                }
                std::cout << "*iteration_counter > " << MAX_ITERATION_NEWTON
                          << ", getting init_approx = " << answer << "* ";

                iteration_counter = 0;
            }
        }

        prev = answer;
        answer = prev - f(prev) / f_d(prev);
        ++iteration_counter;

        if (answer < interval.left_border.x || answer > interval.right_border.x) {
            answer = prev - f(prev) * (prev - chord_border) / (f(prev) - f(chord_border));
            ++iteration_counter;
        }
    }

    return Meta(answer, iteration_counter);
}

//! requires function value at interval points
Meta newton_method_modification_numerical(const std::function<double (double)> &f, double prev,
        const Interval &interval, double precision) {
    double answer = prev - f(prev) / get_derivative(f, prev, precision);
    size_t iteration_counter = 1;

    if (answer < interval.left_border.x || answer > interval.right_border.x) {
        answer = get_chord_method_approx(interval);
        ++iteration_counter;
    }

    double second_derivative_left_border = (get_derivative(f, interval.left_border.x + precision, precision)
        - get_derivative(f, interval.left_border.x, precision)) / precision;

    double chord_border;
    if (f(interval.left_border.x) * second_derivative_left_border > 0) {
        chord_border = interval.left_border.x;
    } else {
        chord_border = interval.right_border.x;
    }

    while (std::abs(answer - prev) >= precision) {
        if (iteration_counter > MAX_ITERATION_NEWTON) {
            std::cout << "*iteration_counter > " << MAX_ITERATION_NEWTON << "* ";
            return Meta(answer, iteration_counter);
        }

        prev = answer;
        answer = prev - f(prev) / get_derivative(f, prev, precision);
        ++iteration_counter;

        if (answer < interval.left_border.x || answer > interval.right_border.x) {
            answer = prev - f(prev) * (prev - chord_border) / (f(prev) - f(chord_border));
            ++iteration_counter;
        }
    }

    return Meta(answer, iteration_counter);
}

//! requires function value at interval points
double get_chord_method_approx(const Interval &interval) {
    return (interval.left_border.y * interval.right_border.x - interval.right_border.y * interval.left_border.x)
           / (interval.left_border.y - interval.right_border.y);
}

Meta newton_method(const std::function<std::vector<double>(std::vector<double>)> &inv_jacobian,
    std::vector<double> prev, const std::vector<Interval> &intervals, double precision) {

    std::vector<double> answer(2);
    answer[0] = prev[0] - inv_jacobian(prev)[0];
    answer[1] = prev[1] - inv_jacobian(prev)[1];
    size_t iteration_counter = 1;

    if (answer[0] < intervals[0].left_border.x || answer[0] > intervals[0].right_border.x ||
        answer[1] < intervals[1].left_border.x || answer[1] > intervals[1].right_border.x) {
        std::cout << "*out of range* ";
        return Meta(answer, iteration_counter);
    }

    while (std::sqrt((answer[0] - prev[0]) * (answer[0] - prev[0]) + (answer[1] - prev[1]) * (answer[1] - prev[1])) >= precision) {
        if (iteration_counter > MAX_ITERATION_NEWTON) {
            std::cout << "*iteration_counter > " << MAX_ITERATION_NEWTON << "* ";
            return Meta(answer, iteration_counter);
        }

        prev = answer;
        answer[0] = prev[0] - inv_jacobian(prev)[0];
        answer[1] = prev[1] - inv_jacobian(prev)[1];
        ++iteration_counter;
        
        if (answer[0] < intervals[0].left_border.x || answer[0] > intervals[0].right_border.x ||
            answer[1] < intervals[1].left_border.x || answer[1] > intervals[1].right_border.x) {
            std::cout << "*out of range* ";
            return Meta(answer, iteration_counter);
        }        
    }

    return Meta(answer, iteration_counter);
}

Meta newton_method_numerical(const std::function<std::vector<double>(std::vector<double>)> &f,
    std::vector<double> prev, const std::vector<Interval> &intervals, double precision) {

    auto inv_jacobian = inverse_jacobian(get_jacobian(f, prev));

    std::vector<double> answer(2);
    answer[0] = prev[0] - (inv_jacobian[0][0] * f(prev)[0] + inv_jacobian[0][1] * f(prev)[1]);
    answer[1] = prev[1] - (inv_jacobian[1][0] * f(prev)[0] + inv_jacobian[1][1] * f(prev)[1]);
    size_t iteration_counter = 1;

    if (answer[0] < intervals[0].left_border.x || answer[0] > intervals[0].right_border.x ||
        answer[1] < intervals[1].left_border.x || answer[1] > intervals[1].right_border.x) {
//        std::cout << "*out of range* ";
        return Meta(answer, MAX_ITERATION_NEWTON);
    }

    while (std::sqrt((answer[0] - prev[0]) * (answer[0] - prev[0]) + (answer[1] - prev[1]) * (answer[1] - prev[1])) >= precision) {
        if (iteration_counter > MAX_ITERATION_NEWTON) {
            std::cout << "*iteration_counter > " << MAX_ITERATION_NEWTON << "* ";
            return Meta(answer, MAX_ITERATION_NEWTON);
        }

        prev = answer;

        inv_jacobian = inverse_jacobian(get_jacobian(f, prev));
        answer[0] = prev[0] - (inv_jacobian[0][0] * f(prev)[0] + inv_jacobian[0][1] * f(prev)[1]);
        answer[1] = prev[1] - (inv_jacobian[1][0] * f(prev)[0] + inv_jacobian[1][1] * f(prev)[1]);
        ++iteration_counter;

        if (answer[0] < intervals[0].left_border.x || answer[0] > intervals[0].right_border.x ||
            answer[1] < intervals[1].left_border.x || answer[1] > intervals[1].right_border.x) {
//            std::cout << "*out of range* ";
            return Meta(answer, MAX_ITERATION_NEWTON);
        }
    }

    return Meta(answer, iteration_counter);
}

void export_newton(const std::string &file_name,
    const std::function<std::vector<double>(std::vector<double>)> &f,
    const std::vector<Interval> &intervals) {
#if (defined (_WIN32) || defined (_WIN64))
    std::ofstream file(file_name);
#else
    std::ofstream file(file_name);
#endif

    size_t number_of_mesh_points = 20;
    auto mesh_x1 = uniform_mesh(intervals[0].left_border.x, intervals[0].right_border.x, number_of_mesh_points);
    auto mesh_x2 = uniform_mesh(intervals[1].left_border.x, intervals[1].right_border.x, number_of_mesh_points);
    
    for (size_t i = 0; i <= number_of_mesh_points; ++i) {
        for (size_t j = 0; j <= number_of_mesh_points; ++j) {
            file << std::setw(10) << mesh_x1[i] << std::setw(10) << mesh_x2[j] << std::setw(8)
                << newton_method_numerical(f, { mesh_x1[i], mesh_x2[j] }, intervals).iteration_counter << std::endl;
        }
    }
}

Meta bisection_method_test(const std::function<double (double)> &f, const Interval &interval, double precision) {
    double double_precision = precision * 2;

    Point answer((interval.left_border.x + interval.right_border.x) / 2);
    size_t iteration_counter = 1;

    Point left_border = interval.left_border;
    Point right_border = interval.right_border;

    double exact = 0.2;
    std::vector<double> err;
    std::vector<double> p;
    err.push_back(std::abs(exact - answer.x));
    p.push_back(0);

    while ((right_border.x - left_border.x) >= double_precision) {
        answer.y = f(answer.x);

        if (std::abs(answer.y) < precision) {
            return Meta(answer.x, iteration_counter);
        }

        if (left_border.y * answer.y < 0) {
            right_border = answer;
        } else {
            left_border = answer;
        }

        answer.x = (left_border.x + right_border.x) / 2;
        ++iteration_counter;

//        err.push_back(std::abs(exact - answer));
//        p.push_back(std::log(*--err.end()) / std::log(*--(--err.end())));
        err.push_back(std::abs(exact - answer.x));
        if (iteration_counter > 2) {
            double x_k_p1 = *--err.end();
            double x_k = *--(--err.end());
            double x_k_m1 = *--(--(--err.end()));
            p.push_back(std::log(x_k_p1 / x_k) / std::log(x_k / x_k_m1));
        } else {
            p.push_back(0);
        }
    }

    for (int i = 0; i < err.size(); ++i) {
        std::cout << i + 1 << " & " << err[i] << " & " << p[i] << '\n';
    }

    return Meta(answer.x, iteration_counter);
}

Meta newton_method_test(const std::function<double (double)> &f,
                             double prev, const Interval &interval, double precision) {
    double answer = prev - f(prev) / get_derivative(f, prev, precision);
    size_t iteration_counter = 1;

    if (answer < interval.left_border.x || answer > interval.right_border.x) {
        std::cout << "*out of range* ";
        return Meta(answer, iteration_counter);
    }

    double exact = 0.2;
    std::vector<double> err;
    std::vector<double> p;
    err.push_back(std::abs(exact - answer));
    p.push_back(0);

    while (std::abs(answer - prev) >= precision) {
        if (iteration_counter > MAX_ITERATION_NEWTON) {
            std::cout << "*iteration_counter > " << MAX_ITERATION_NEWTON << "* ";
            return Meta(answer, iteration_counter);
        }

        prev = answer;
        answer = prev - f(prev) / get_derivative(f, prev, precision);
        ++iteration_counter;

        if (answer < interval.left_border.x || answer > interval.right_border.x) {
            std::cout << "*out of range* ";
            return Meta(answer, iteration_counter);
        }

        err.push_back(std::abs(exact - answer));
        if (iteration_counter > 2) {
            double x_k_p1 = *--err.end();
            double x_k = *--(--err.end());
            double x_k_m1 = *--(--(--err.end()));
            p.push_back(std::log(x_k_p1 / x_k) / std::log(x_k / x_k_m1));
        } else {
            p.push_back(0);
        }
//        err.push_back(std::abs(exact - answer));
//        p.push_back(std::log(*--err.end()) / std::log(*--(--err.end())));
    }

    for (int i = 0; i < err.size(); ++i) {
        std::cout << err[i]  << ' ' << p[i] << '\n';
    }

    return Meta(answer, iteration_counter);
}