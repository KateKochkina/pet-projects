#include <iostream>

#include "tests.hpp"
#include "utils.h"
#include "find_root.h"

int main() {
    {
        std::cout << "test 1\n";

        Interval root_interval = {0, 1};

        auto intervals_with_roots = get_intervals_with_roots(function_test1, root_interval, 11);

        std::cout << "Bisection method:\n";
        for (const auto &interval : intervals_with_roots) {
            auto root = bisection_method(function_test1, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), root = " << root.x[0] << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method (analytic derivative):\n";
        for (const auto &interval : intervals_with_roots) {
            //            double init_approx = get_chord_method_approx(interval);
            double init_approx = (interval.left_border.x + interval.right_border.x) / 2;

            auto root = newton_method(function_test1, function_test1_d, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method (numerical derivative):\n";
        for (const auto &interval : intervals_with_roots) {
            //            double init_approx = get_chord_method_approx(interval);
            double init_approx = (interval.left_border.x + interval.right_border.x) / 2;

            auto root = newton_method_numerical(function_test1, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method modification (analytic derivative):\n";
        for (const auto &interval : intervals_with_roots) {
            double init_approx = (interval.left_border.x + interval.right_border.x) / 2;

            auto root = newton_method_modification(function_test1, function_test1_d, function_test1_dd, init_approx,
                                                   interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method modification (numerical derivative):\n";
        for (const auto &interval : intervals_with_roots) {
            double init_approx = (interval.left_border.x + interval.right_border.x) / 2;

            auto root = newton_method_modification_numerical(function_test1, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << std::endl;
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
        std::cout << "test 2\n";

        Interval root_interval = {-1, 10};

        //! standard
        auto intervals_with_roots = get_intervals_with_roots(function_test2, root_interval, 17);

        //! arbitrary interval (don't forget to initialize function value fields)
        root_interval.left_border.y = function_test2(root_interval.left_border.x);
        root_interval.right_border.y = function_test2(root_interval.right_border.x);
        intervals_with_roots[0] = root_interval;

        std::cout << "Bisection method:\n";
        for (const auto &interval : intervals_with_roots) {
            auto root = bisection_method(function_test2, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), root = " << root.x[0] << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method (analytic derivative):\n";
        for (const auto &interval : intervals_with_roots) {
            double init_approx = 8;

            auto root = newton_method(function_test2, function_test2_d, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method (numerical derivative):\n";
        for (const auto &interval : intervals_with_roots) {
            double init_approx = 8;

            auto root = newton_method_numerical(function_test2, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method modification (analytic derivative):\n";
        for (const auto &interval : intervals_with_roots) {
            double init_approx = 8;

            auto root = newton_method_modification(function_test2, function_test2_d, function_test2_dd, init_approx,
                                                   interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method modification (numerical derivative):\n";
        for (const auto &interval : intervals_with_roots) {
            double init_approx = 8;

            auto root = newton_method_modification_numerical(function_test2, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << std::endl;
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
        std::cout << "test 3\n";

        Interval interval = {0, 1};
        interval.left_border.y = function_test3(interval.left_border.x);
        interval.right_border.y = function_test3(interval.right_border.x);

        // auto intervals_with_roots = get_intervals_with_roots(function_test3, root_interval, 20);

        std::cout << "Bisection method:\n";
        {
            auto root = bisection_method(function_test3, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), root = " << root.x[0] << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method (analytic derivative):\n";
        {
            double init_approx = 0;

            auto root = newton_method(function_test3, function_test3_d, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method (numerical derivative):\n";
        {
            double init_approx = 0;

            auto root = newton_method_numerical(function_test3, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method modification (analytic derivative):\n";
        {
            double init_approx = 0;

            auto root = newton_method_modification(function_test3, function_test3_d, function_test3_dd, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method modification (numerical derivative):\n";
        {
            double init_approx = 0;

            auto root = newton_method_modification_numerical(function_test3, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << std::endl;
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
        std::cout << "test 4";

        std::vector<Interval> intervals = {{-10, 10}, {-10, 10}};
        intervals[0].left_border.y = function_test4({intervals[0].left_border.x, intervals[1].left_border.x})[0];
        intervals[0].right_border.y = function_test4({intervals[0].right_border.x, intervals[1].right_border.x})[0];
        intervals[1].left_border.y = function_test4({intervals[0].left_border.x, intervals[1].left_border.x})[1];
        intervals[1].right_border.y = function_test4({intervals[0].right_border.x, intervals[1].right_border.x})[1];

        // auto intervals_with_roots = get_intervals_with_roots(function_test3, root_interval, 20);

        std::vector<std::vector<double>> init_approximate = {{5, 1}, {-5, 1}};

        for (auto init_approx : init_approximate) {
            std::cout << "\nNewton method (analytic derivative):\n";
            {
                auto root = newton_method(function_test4_newton_an, init_approx, intervals);

                std::cout << "interval: (" << intervals[0].left_border.x << ", " << intervals[0].right_border.x
                          << ") (" << intervals[1].left_border.x << ", " << intervals[1].right_border.x
                          << "), init approx = (" << init_approx[0] << ", " << init_approx[1] << "), root = (" << root.x[0]
                          << ", " << root.x[1] << "), number of iterations: " << root.iteration_counter << "\n";
            }

            std::cout << "\nNewton method (numerical derivative):\n";
            {
                auto root = newton_method_numerical(function_test4, init_approx, intervals);

                std::cout << "interval: (" << intervals[0].left_border.x << ", " << intervals[0].right_border.x
                          << ") (" << intervals[1].left_border.x << ", " << intervals[1].right_border.x
                          << "), init approx = (" << init_approx[0] << ", " << init_approx[1] << "), root = (" << root.x[0]
                          << ", " << root.x[1] << "), number of iterations: " << root.iteration_counter << "\n";
            }
        }

        std::cout << std::endl;
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
        std::cout << "test 5";

        std::vector<Interval> intervals = {{-10, 10}, {-10, 10}};
        intervals[0].left_border.y = function_test5({intervals[0].left_border.x, intervals[1].left_border.x})[0];
        intervals[0].right_border.y = function_test5({intervals[0].right_border.x, intervals[1].right_border.x})[0];
        intervals[1].left_border.y = function_test5({intervals[0].left_border.x, intervals[1].left_border.x})[1];
        intervals[1].right_border.y = function_test5({intervals[0].right_border.x, intervals[1].right_border.x})[1];

        // auto intervals_with_roots = get_intervals_with_roots(function_test3, root_interval, 20);

        std::vector<std::vector<double>> init_approximate = {{-4, 0}, {0, -4}, {1, 3}, {3, 1}};

        for (auto init_approx : init_approximate) {
            std::cout << "\nNewton method (analytic derivative):\n";
            {
                auto root = newton_method(function_test5_newton_an, init_approx, intervals);

                std::cout << "interval: (" << intervals[0].left_border.x << ", " << intervals[0].right_border.x
                          << ") (" << intervals[1].left_border.x << ", " << intervals[1].right_border.x
                          << "), init approx = (" << init_approx[0] << ", " << init_approx[1] << "), root = (" << root.x[0]
                          << ", " << root.x[1] << "), number of iterations: " << root.iteration_counter << "\n";
            }

            std::cout << "\nNewton method (numerical derivative):\n";
            {
                auto root = newton_method_numerical(function_test5, init_approx, intervals);

                std::cout << "interval: (" << intervals[0].left_border.x << ", " << intervals[0].right_border.x
                          << ") (" << intervals[1].left_border.x << ", " << intervals[1].right_border.x
                          << "), init approx = (" << init_approx[0] << ", " << init_approx[1] << "), root = (" << root.x[0]
                          << ", " << root.x[1] << "), number of iterations: " << root.iteration_counter << "\n";
            }
        }

        std::cout << std::endl;
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
        std::cout << "equation 20\n";

        Interval interval = {-1, 0};
        interval.left_border.y = function_variant20(interval.left_border.x);
        interval.right_border.y = function_variant20(interval.right_border.x);

        // auto intervals_with_roots = get_intervals_with_roots(function_test3, root_interval, 20);

        std::cout << "Bisection method:\n";
        {
            auto root = bisection_method(function_variant20, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), root = " << root.x[0] << ", number of iterations: " << root.iteration_counter << "\n";
        }

        //        std::cout << "\nNewton method (analytic derivative):\n";
        //        {
        //            double init_approx = 0;
        //
        //            auto root = newton_method(function_variant20, function_variant20_d, init_approx, interval);
        //
        //            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
        //                      << "), init approx = " << init_approx << ", root = " << root.x[0]
        //                      << ", number of iterations: " << root.iteration_counter << "\n";
        //        }

        std::cout << "\nNewton method (numerical derivative):\n";
        {
            double init_approx = -0.9;

            auto root = newton_method_numerical(function_variant20, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        //        std::cout << "\nNewton method modification (analytic derivative):\n";
        //        {
        //            double init_approx = 0;
        //
        //            auto root = newton_method_modification(function_variant20, function_variant20_d, function_variant20_dd, init_approx, interval);
        //
        //            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
        //                      << "), init approx = " << init_approx << ", root = " << root.x[0]
        //                      << ", number of iterations: " << root.iteration_counter << "\n";
        //        }

        std::cout << "\nNewton method modification (numerical derivative):\n";
        {
            double init_approx = -0.9;

            auto root = newton_method_modification_numerical(function_variant20, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << std::endl;
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
        std::cout << "system 20";

        std::vector<Interval> intervals = {{-10, 10}, {-10, 10}};
        intervals[0].left_border.y = system_variant20({intervals[0].left_border.x, intervals[1].left_border.x})[0];
        intervals[0].right_border.y = system_variant20({intervals[0].right_border.x, intervals[1].right_border.x})[0];
        intervals[1].left_border.y = system_variant20({intervals[0].left_border.x, intervals[1].left_border.x})[1];
        intervals[1].right_border.y = system_variant20({intervals[0].right_border.x, intervals[1].right_border.x})[1];

        //         auto intervals_with_roots = get_intervals_with_roots(function_test3, root_interval, 20);

        std::vector<std::vector<double>> init_approximate = {{-4, 1}, {3, 3}};

        for (auto init_approx : init_approximate) {
            std::cout << "\nNewton method (analytic derivative):\n";
            {
                auto root = newton_method(system_variant20_newton_an, init_approx, intervals);

                std::cout << "interval: (" << intervals[0].left_border.x << ", " << intervals[0].right_border.x
                          << ") (" << intervals[1].left_border.x << ", " << intervals[1].right_border.x
                          << "), init approx = (" << init_approx[0] << ", " << init_approx[1] << "), root = (" << root.x[0]
                          << ", " << root.x[1] << "), number of iterations: " << root.iteration_counter << "\n";
            }

            std::cout << "\nNewton method (numerical derivative):\n";
            {
                auto root = newton_method_numerical(system_variant20, init_approx, intervals);

                std::cout << "interval: (" << intervals[0].left_border.x << ", " << intervals[0].right_border.x
                          << ") (" << intervals[1].left_border.x << ", " << intervals[1].right_border.x
                          << "), init approx = (" << init_approx[0] << ", " << init_approx[1] << "), root = (" << root.x[0]
                          << ", " << root.x[1] << "), number of iterations: " << root.iteration_counter << "\n";
            }
        }

        std::cout << std::endl;
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //    {
    //        std::vector<Interval> intervals = { {-10, 10}, {-10, 10} };
    //        export_newton("test4.dat", function_test4, intervals);
    //    }
    //    {
    //        std::vector<Interval> intervals = { {-10, 10}, {-10, 10} };
    //        export_newton("test5.dat", function_test5, intervals);
    //    }
    //    {
    //        std::vector<Interval> intervals = { {-10, 10}, {-10, 10} };
    //        export_newton("system20.dat", system_variant20, intervals);
    //    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::cout << "Оценка скорости сходимости:\n";
    {
        Interval interval = {0, 1};
        interval.left_border.y = function_test3(interval.left_border.x);
        interval.right_border.y = function_test3(interval.right_border.x);

        std::cout << "Bisection method:\n";
        {
            auto root = bisection_method_test(function_test3, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), root = " << root.x[0] << ", number of iterations: " << root.iteration_counter << "\n";
        }

        std::cout << "\nNewton method (numerical derivative):\n";
        {
            double init_approx = 0.5;

            auto root = newton_method_test(function_test3, init_approx, interval);

            std::cout << "interval (" << interval.left_border.x << ", " << interval.right_border.x
                      << "), init approx = " << init_approx << ", root = " << root.x[0]
                      << ", number of iterations: " << root.iteration_counter << "\n";
        }
    }

    #if (defined (_WIN32) || defined (_WIN64))
        system("PAUSE");
    #endif

    return 0;
}