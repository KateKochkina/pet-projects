#pragma once
#include"SLAE_solve.h"

typedef double(*func) (double);
typedef double(*kernel) (double, double);

struct Problem_data {
	double x_left;
	double x_right;
	func right_side_func;
	kernel kernel_func;
	double lambda;

	Problem_data(double x_left_, double x_right_, func right_side_func_, kernel kernel_func_, double lambda_) {
		x_left = x_left_;
		x_right = x_right_;
		right_side_func = right_side_func_;
		kernel_func = kernel_func_;
		lambda = lambda_;
	}
};

struct Options {
	double h;
	int DIM_x;

	Options(double h_, Problem_data data) {
		h = h_;
		DIM_x = (int)trunc((data.x_right - data.x_left) / h) + 1;
	}
};

void Trapezoid_quadrature_coefficient(TVector<double>& quadrature_coeff, double h);

void Quadrature_method(std::string problem_name, Problem_data data, Options opt, TMatrix<double>& system_matr,
	TVector<double>& right_side_vec, TVector<double>& quadrature_coeff, TVector<double>& solution, void (*Quadratute_formula) (TVector<double>&, double));

void Simple_Iteration_method(std::string problem_name, Problem_data data, Options opt, TVector<double>& u_prev,
	TVector<double>& quadrature_coeff, TVector<double>& solution, void (*Quadratute_formula) (TVector<double>&, double), double TOL);

void Degenerate_kernel_equation(std::string problem_name, Problem_data data, Options opt, TMatrix<double>& system_matr,
	TVector<double>& right_side_vec, TVector<double>& quadrature_coeff, TVector<double>& solution, void (*Quadratute_formula) (TVector<double>&, double),
	TVector<func>& psi, TVector<func>& phi, TVector<double>& coeffs_C);//dim(system_matr) = dim(right_side_vec) = dim(coeffs_C) = dim(psi) = dim(phi) = m; dim(quadrature_coeff) = dim(solution) = N

void Singular_kernel_equation(std::string problem_name, TVector<double>(*Q_sing) (TVector<double>, TVector<double>), double (*right_side_func) (int, double), int N, int var);