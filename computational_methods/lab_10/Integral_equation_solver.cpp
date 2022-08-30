#include"Integral_equation_solver.h"

void filing_mesh(std::string filename, Options opt, Problem_data data) {
	std::ofstream fout;
	fout.open(filename);
	for (auto j = 0; j < opt.DIM_x; ++j) {
		fout << std::setprecision(16) << std::fixed << data.x_left + j * opt.h << '\t';
	}
	fout << std::endl;
	fout.close();
	fout.clear();
}

void filing_matr(std::string filename, TMatrix<double>& matr) {
	std::ofstream fout;
	fout.open(filename);
	for (auto i = 0; i < matr.get_nrows(); ++i) {
		for (auto j = 0; j < matr.get_ncols(); ++j) {
			fout << std::setprecision(16) << std::fixed << matr(i, j) << '\t';
		}
		fout << std::endl;
	}
	fout.close();
	fout.clear();
}

void filing_circle_mesh(std::string filename, TVector<double>* points, int N) {
	std::ofstream fout;
	fout.open(filename);
	for (auto j = 0; j < N; ++j) {
		fout << std::setprecision(16) << std::fixed << points[j][0] << '\t' << points[j][1] << '\n';
	}
	fout.close();
	fout.clear();
}

void Trapezoid_quadrature_coefficient(TVector<double>& quadrature_coeff, double h) {
	quadrature_coeff[0] = 0.5 * h;
	quadrature_coeff[quadrature_coeff.get_dim() - 1] = 0.5 * h;
	for (auto i = 1; i < quadrature_coeff.get_dim() - 1; ++i) {
		quadrature_coeff[i] = h;
	}
}

void Quadrature_method(std::string problem_name, Problem_data data, Options opt, TMatrix<double>& system_matr, 
	TVector<double>& right_side_vec, TVector<double>& quadrature_coeff, TVector<double>& solution, void (*Quadratute_formula) (TVector<double>&, double)) {
	solution.zero_vector();
	clear_file(problem_name + "_QM_mesh.txt");
	filing_mesh(problem_name + "_QM_mesh.txt", opt, data);
	Quadratute_formula(quadrature_coeff, opt.h);
	double x_i;
	for (auto i = 0; i < system_matr.get_nrows(); ++i) {
		x_i = data.x_left + i * opt.h;
		for (auto k = 0; k < system_matr.get_ncols(); ++k) {
			system_matr(i, k) = -data.lambda * quadrature_coeff[k] * data.kernel_func(x_i, data.x_left + k * opt.h) ;
		}
		system_matr(i, i) += 1;
		right_side_vec[i] = data.right_side_func(x_i);
	}
	GaussMeth(problem_name + "_QM_solution.txt", system_matr, right_side_vec, solution, 'y');
	std::cout << "Quadrature method\n";
}

void Simple_Iteration_method(std::string problem_name, Problem_data data, Options opt, TVector<double>& u_prev,
	TVector<double>& quadrature_coeff, TVector<double>& solution, void (*Quadratute_formula) (TVector<double>&, double), double TOL) {
	solution.zero_vector();
	clear_file(problem_name + "_SIM_mesh.txt");
	filing_mesh(problem_name + "_SIM_mesh.txt", opt, data);
	Quadratute_formula(quadrature_coeff, opt.h);
	double quadrature_sum = 0.0;
	int iters = 0;
	double x_i;
	solution = u_prev;
	do {
		swap(solution, u_prev);
		for (auto i = 0; i < u_prev.get_dim(); ++i) {
			x_i = data.x_left + i * opt.h;
			for (auto k = 0; k < u_prev.get_dim(); ++k) {
				quadrature_sum += quadrature_coeff[k] * data.kernel_func(x_i, data.x_left + k * opt.h) * u_prev[k];
			}
			solution[i] = data.right_side_func(x_i) + data.lambda * quadrature_sum;
			quadrature_sum = 0.0;
		}
		iters++;
	} while ((solution - u_prev).norm('3') > TOL);
	filing(problem_name + "_SIM_solution.txt", solution);
	u_prev.zero_vector();
	std::cout << "Simple Iteration method: TOL = " << TOL << ", iters = " << iters << "\n";
}

void Degenerate_kernel_equation(std::string problem_name, Problem_data data, Options opt, TMatrix<double>& system_matr,
	TVector<double>& right_side_vec, TVector<double>& quadrature_coeff, TVector<double>& solution, void (*Quadratute_formula) (TVector<double>&, double), 
	TVector<func>& psi, TVector<func>& phi, TVector<double>& coeffs_C) {
	solution.zero_vector();
	clear_file(problem_name + "_DKE_mesh.txt");
	filing_mesh(problem_name + "_DKE_mesh.txt", opt, data);
	Quadratute_formula(quadrature_coeff, opt.h);
	//matrix mxm construction
	int m = psi.get_dim();
	double alpha_ij;
	double beta_i;
	double s_k;
	for (auto i = 0; i < m; ++i) {
		for (auto j = 0; j < m; ++j) {
			alpha_ij = 0.0;
			for (auto k = 0; k < quadrature_coeff.get_dim(); ++k) {
				s_k = data.x_left + k * opt.h;
				alpha_ij += quadrature_coeff[k] * psi[i](s_k) * phi[j](s_k);
			}
			system_matr(i, j) = -data.lambda * alpha_ij;
		}
		system_matr(i, i) += 1;
		beta_i = 0.0;
		for (auto k = 0; k < quadrature_coeff.get_dim(); ++k) {
			s_k = data.x_left + k * opt.h;
			beta_i += quadrature_coeff[k] * psi[i](s_k) * data.right_side_func(s_k);
		}
		right_side_vec[i] = beta_i;
	}
	clear_file(problem_name + "_DKE_matr.txt");
	filing_matr(problem_name + "_DKE_matr.txt", system_matr);
	coeffs_C = GaussMeth(problem_name + "_DKE_solution.txt", system_matr, right_side_vec, 'n');
	//solution construction
	double x_k;
	for (auto k = 0; k < solution.get_dim(); ++k) {
		x_k = data.x_left + k * opt.h;
		for (auto i = 0; i < m; ++i) {
			solution[k] += coeffs_C[i] * phi[i](x_k);
		}
		solution[k] *= data.lambda;
		solution[k] += data.right_side_func(x_k);
	}
	filing(problem_name + "_DKE_solution.txt", solution);
	std::cout << "Degenerate kernel equation\n";
}

void Singular_kernel_equation(std::string problem_name, TVector<double> (*Q_sing) (TVector<double>, TVector<double>), double (*right_side_func) (int, double), int N, int var) {
	int sol_dim = N + 1;
	TMatrix<double> system_matr(sol_dim, sol_dim);
	TVector<double> right_side_vec(sol_dim);
	TVector<double>* k_points = new TVector<double>[N];
	TVector<double>* c_points = new TVector<double>[N];
	const double PI = 3.14159265358979323846;
	double dl = 2 * PI / N;
	for (auto i = 0; i < N; ++i) {
		k_points[i][0] = cos(2 * PI * (i + 0.5) / N);
		k_points[i][1] = sin(2 * PI * (i + 0.5) / N);
		c_points[i][0] = cos(2 * PI * i / N);
		c_points[i][1] = sin(2 * PI * i / N);
	}
	clear_file(problem_name + "_SKE_k_points.txt");
	clear_file(problem_name + "_SKE_c_points.txt");
	filing_circle_mesh(problem_name + "_SKE_k_points.txt", k_points, N);
	filing_circle_mesh(problem_name + "_SKE_c_points.txt", c_points, N);
	TVector<double> n(2);
	for (auto i = 0; i < N; ++i) {
		n[0] = cos(2 * PI * (i + 0.5) / N);
		n[1] = sin(2 * PI * (i + 0.5) / N);
		for (auto j = 0; j < N; ++j) {
			system_matr(i, j) = scalar_product(n, Q_sing(k_points[i], c_points[j])) * dl;
		}
		right_side_vec[i] = right_side_func(var, dl * (i + 0.5));
		system_matr(sol_dim - 1, i) = dl; //dl_i, but dl_i = dl, so doesn't matter ( right_side_vec[sol_dim - 1] = 0 )
		system_matr(i, sol_dim - 1) = 1;
	}
	TVector<double> solution = GaussMeth(problem_name + "_SKE_solution.txt", system_matr, right_side_vec, 'y');
	std::cout << "Singular kernel equation\n";
	delete[] k_points;
	delete[] c_points;
}