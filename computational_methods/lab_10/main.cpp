#include<iostream>
#include<cmath>
#include"Integral_equation_solver.h"
#include"Test_functions.h"
#include<string>

int main() {
    //------------------QUADRATURE METHOD AND SIMPLE ITERATION METHOD------------------

    std::string export_path = "/Users/artem/Desktop/study/вычи/5/export/";
    double x_left = 0.0;
    double x_right = 1.0;
    double lambda = 1.0;
    double h = 0.01;

    std::string problem_name = export_path + "test1";
    Problem_data data1(x_left, x_right, f1, K1, lambda);
    Options opt1(h, data1);
    TMatrix<double> system_matr(opt1.DIM_x, opt1.DIM_x);
    TVector<double> right_side_vec(opt1.DIM_x);
    TVector<double> quadrature_coeff(opt1.DIM_x);
    TVector<double> solution(opt1.DIM_x);
    TVector<double> u_0(opt1.DIM_x);
    double TOL = 0.1;

//    Quadrature_method(problem_name + "_h_" + std::to_string(opt1.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f1",
//                      data1, opt1, system_matr, right_side_vec, quadrature_coeff, solution, Trapezoid_quadrature_coefficient);
//    for (auto i = 1; i < 7; ++i) {
//        Simple_Iteration_method(problem_name + "_h_" + std::to_string(opt1.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f1_TOL_" + std::to_string(i + 1),
//                                data1, opt1, u_0, quadrature_coeff, solution, Trapezoid_quadrature_coefficient, TOL / pow(10, i));
//    }
//    data1.right_side_func = g1;
//    Quadrature_method(problem_name + "_h_" + std::to_string(opt1.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f2",
//                      data1, opt1, system_matr, right_side_vec, quadrature_coeff, solution, Trapezoid_quadrature_coefficient);
//    for (auto i = 1; i < 7; ++i) {
//        Simple_Iteration_method(problem_name + "_h_" + std::to_string(opt1.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f2_TOL_" + std::to_string(i + 1),
//                                data1, opt1, u_0, quadrature_coeff, solution, Trapezoid_quadrature_coefficient, TOL / pow(10, i));
//    }
//

    problem_name = export_path + "var17";
    Problem_data data17(x_left, x_right, f17, K17, lambda);
    Options opt17(h, data17);

//    Quadrature_method(problem_name + "_h_" + std::to_string(opt17.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f1",
//    	data17, opt17, system_matr, right_side_vec, quadrature_coeff, solution, Trapezoid_quadrature_coefficient);
//    for (auto i = 1; i < 7; ++i) {
//    	Simple_Iteration_method(problem_name + "_h_" + std::to_string(opt17.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f1_TOL_" + std::to_string(i + 1),
//    		data17, opt17, u_0, quadrature_coeff, solution, Trapezoid_quadrature_coefficient, TOL / pow(10, i));
//    }
//    data17.right_side_func = g17;
//    Quadrature_method(problem_name + "_h_" + std::to_string(opt17.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f2",
//    	data17, opt17, system_matr, right_side_vec, quadrature_coeff, solution, Trapezoid_quadrature_coefficient);
//    for (auto i = 1; i < 7; ++i) {
//    	Simple_Iteration_method(problem_name + "_h_" + std::to_string(opt17.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f2_TOL_" + std::to_string(i + 1),
//    		data17, opt17, u_0, quadrature_coeff, solution, Trapezoid_quadrature_coefficient, TOL / pow(10, i));
//    }


    //------------------DEGENERATE KERNEL EQUATION------------------
    
//    problem_name = export_path + "test1";
//    TVector<func> phi = { phi0, phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8 };
//    TVector<func> psi = { psi0, psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8 };
//    TMatrix<double> system_matr1(phi.get_dim(), phi.get_dim());
//    TVector<double> right_side_vec1(phi.get_dim());
//    TVector<double> coeffs_C(phi.get_dim());
//    for (auto i = 0; i <= 7; ++i) {
//    	data1.right_side_func = f1;
//    	Degenerate_kernel_equation(problem_name + "_h_" + std::to_string(opt1.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f1_dim_" + std::to_string(phi.get_dim()),
//    		data1, opt1, system_matr1, right_side_vec1, quadrature_coeff, solution, Trapezoid_quadrature_coefficient, psi, phi, coeffs_C);
//    	//data1.right_side_func = g1;
//    	//Degenerate_kernel_equation(problem_name + "_h_" + std::to_string(opt1.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f2_dim_" + std::to_string(phi.get_dim()),
//    	//	data1, opt1, system_matr1, right_side_vec1, quadrature_coeff, solution, Trapezoid_quadrature_coefficient, psi, phi, coeffs_C);
//    	phi.remove_last_elems(1);
//    	psi.remove_last_elems(1);
//    	right_side_vec1.remove_last_elems(1);
//    	coeffs_C.remove_last_elems(1);
//    	system_matr1.remove_last_cols(1);
//    	system_matr1.remove_last_rows(1);
//    }
//
//    problem_name = export_path + "var17";
//    TVector<func> phi_var17 = { phi0_var17, phi1_var17, phi2_var17, phi3_var17, phi4_var17, phi5_var17, phi6_var17, phi7_var17, phi8_var17 };
//    TVector<func> psi_var17 = { psi0_var17, psi1_var17, psi2_var17, psi3_var17, psi4_var17, psi5_var17, psi6_var17, psi7_var17, psi8_var17 };
//    TMatrix<double> system_matr17(phi_var17.get_dim(), phi_var17.get_dim());
//    TVector<double> right_side_vec17(phi_var17.get_dim());
//    TVector<double> coeffs_C_var17(phi_var17.get_dim());
//    for (auto i = 0; i <= 7; ++i) {
//    	Degenerate_kernel_equation(problem_name + "_h_" + std::to_string(opt17.h) + "_a_" + std::to_string(x_left) + "_b_" + std::to_string(x_right) + "_f1_dim_" + std::to_string(phi_var17.get_dim()),
//    		data17, opt17, system_matr17, right_side_vec17, quadrature_coeff, solution, Trapezoid_quadrature_coefficient, psi_var17, phi_var17, coeffs_C_var17);
//    	phi_var17.remove_last_elems(1);
//    	psi_var17.remove_last_elems(1);
//    	right_side_vec17.remove_last_elems(1);
//    	coeffs_C_var17.remove_last_elems(1);
//    	system_matr17.remove_last_cols(1);
//    	system_matr17.remove_last_rows(1);
//    }


    //------------------SINGULAR KERNEL EQUATION------------------
    
    TVector<int> pts = { 2, 3, 4, 5, 9, 10, 20, 200 };
    int var = 17;
    for (auto i = 0; i < pts.get_dim(); ++i) {
    	problem_name = export_path + "var17_pts_" + std::to_string(pts[i]);
    	Singular_kernel_equation(problem_name, Q_sing, f_odd, pts[i], var);
    }

    return 0;
}