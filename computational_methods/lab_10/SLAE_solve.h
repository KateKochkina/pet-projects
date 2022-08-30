#pragma once
#include <iostream>
#include <iomanip>
#include "Gauss_Method.h"
#include <string>

const int MAX_ITER = 2000;
const double EPS_ZERO = 1e-11;

void filing(std::string filename, double norm);

template <typename MYTYPE> void filing(std::string filename, TVector<MYTYPE> vec) {
	std::ofstream fout;
	fout.open(filename);
	fout << vec.get_dim() << '\n';
	for (auto j = 0; j < vec.get_dim(); ++j) {
		fout << std::setprecision(20) << std::fixed << vec[j] << '\n';
	}
	fout << std::endl;
	fout.close();
	fout.clear();
}

void clear_file(std::string filename);

template <typename MYTYPE> void Jacobi(std::string filename, TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<MYTYPE>& x0, TVector<MYTYPE>& x_real, int power) {
	double TOL = pow(10, (double)power);
	TMatrix<MYTYPE> C(A);
	for (auto i = 0; i < C.get_nrows(); ++i) {
		for (auto j = 0; j < C.get_ncols(); ++j) {
			C(i, j) /= -A(i, i);
		}
		C(i, i) = MYTYPE(0);
	}
	TVector<MYTYPE> y(b);
	for (auto i = 0; i < y.get_dim(); ++i) {
		y[i] /= A(i, i);
	}
	TVector<MYTYPE> x_next(x0);
	TVector<MYTYPE> x_prev(x0);
	int iter = 0;
	do {
		iter++;
		x_prev = x_next;
		x_next = C * x_prev + y;
		filing("norm_Jacobi_" + filename + ".txt", iter);
		filing("norm_Jacobi_" + filename + ".txt", (x_real - x_next).norm('2'));
	} while ((x_next - x_prev).norm('2') > TOL);
	std::cout << "Jacobi " << std::to_string(power) << ": iter = " << iter << "\nr_next.norm('2') = " << (x_real - x_next).norm('2') << "\n";

	clear_file("sol_Jacobi_" + filename + ".txt");
	filing("sol_Jacobi_" + filename + ".txt", x_next);
}

template <typename MYTYPE> void BICGSTAB(std::string filename, TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<MYTYPE>& x0, TVector<MYTYPE>& x_real, int power) {
	double TOL = pow(10, (double)power);
	TVector<MYTYPE> r_prev = b - A * x0;
	TVector<MYTYPE> r_next(r_prev);
	TVector<MYTYPE> r0_ast(r_prev); //arbitrary: (r0_ast,r0) != 0
	r0_ast.single_vector();
	TVector<MYTYPE> p(r_prev);
	TVector<MYTYPE> x(x0);
	int iter = 0;
	MYTYPE rho_prev = scalar_product(r_prev, r0_ast);
	clear_file("norm_BiCGStab_" + filename + "_" + std::to_string(power) + ".txt");
	while (iter < MAX_ITER && r_next.norm('3') >= TOL) {
		iter++;
		TVector<MYTYPE> A_mult_p = A * p;
		MYTYPE alpha = rho_prev / scalar_product(A_mult_p, r0_ast);
		TVector<MYTYPE> s = r_prev - alpha * A_mult_p;
		TVector<MYTYPE> A_mult_s = A * s;
		MYTYPE omega = scalar_product(A_mult_s, s) / scalar_product(A_mult_s, A_mult_s);
		x += alpha * p + omega * s;
		r_next = s - omega * A_mult_s;
		filing("norm_BiCGStab_" + filename + "_" + std::to_string(power) + ".txt", (x_real - x).norm('3'));
		MYTYPE rho_next = scalar_product(r_next, r0_ast);
		if (std::abs(rho_next) < EPS_ZERO) {
			std::cout << "\nBREAK! ";
			std::cout << iter << "\n";
			break;
		}
		MYTYPE beta = (alpha * rho_next) / (omega * rho_prev);
		p = r_next + beta * (p - omega * A_mult_p);
		r_prev = r_next;
		rho_prev = rho_next;
	}
	std::cout << "BISGSTAB " << std::to_string(power) << ": iter = " << iter << "\nr_next.norm('3') = " << r_next.norm('3') << "\n";

	clear_file("sol_BiCGStab_" + filename + ".txt");
	filing("sol_BiCGStab_" + filename + ".txt", x);
}

template <typename MYTYPE> TVector<MYTYPE> small_syst_sol1(TMatrix<MYTYPE>& A, TVector<MYTYPE>& vec) {
	int DIM = vec.get_dim();
	int eq_num;
	TMatrix<MYTYPE> R(3, 3);
	TVector<MYTYPE> v = { - vec[DIM - 3], - vec[DIM - 2], - vec[DIM - 1] };
	TVector<MYTYPE> b(3);
	TVector<MYTYPE> c(3);
	MYTYPE a_i, v_i;
	for (auto i = 0; i < DIM - 3; ++i) {
		if (i % 2 == 0) {
			eq_num = i + 1;
		}
		else {
			eq_num = i - 1;
		}
		a_i = A(eq_num, i);
		b[0] = A(eq_num, DIM - 3); b[1] = A(eq_num, DIM - 2);  b[2] = A(eq_num, DIM - 1);
		c[0] = A(DIM - 3, i); c[1] = A(DIM - 2, i); c[2] = A(DIM - 1, i);
		v_i = vec[eq_num];
		for (auto j = 0; j < 3; ++j) {
			for (auto k = 0; k < 3; ++k) {
				R(j, k) += c[j] * b[k] / a_i;
			}
			v[j] += c[j] * v_i / a_i;
		}
	}
	TVector<MYTYPE> R_sol(3);
	GaussMeth(" ", R, v, R_sol);
	TVector<MYTYPE> sol(DIM);
	for (auto i = 0; i < DIM - 3; ++i) {
		if (i % 2 == 0) {
			eq_num = i + 1;
		}
		else {
			eq_num = i - 1;
		}
		a_i = A(eq_num, i);
		b[0] = A(eq_num, DIM - 3); b[1] = A(eq_num, DIM - 2); b[2] = A(eq_num, DIM - 1);
		v_i = vec[eq_num];
		sol[i] = (v_i - scalar_product(R_sol, b)) / a_i;
	}
	sol[DIM - 3] = R_sol[0];  sol[DIM - 2] = R_sol[1]; sol[DIM - 1] = R_sol[2];
	return sol;
}

template <typename MYTYPE> TVector<MYTYPE> small_syst_sol2(TMatrix<MYTYPE>& A, TVector<MYTYPE>& vec) {
	int DIM = vec.get_dim();
	TMatrix<MYTYPE> R(3, 3);
	TVector<MYTYPE> v = { -vec[DIM - 3], -vec[DIM - 2], -vec[DIM - 1] };
	TVector<MYTYPE> b(3);
	TVector<MYTYPE> c(3);
	MYTYPE a_i, v_i;
	for (auto i = 0; i < DIM - 3; ++i) {
		a_i = A(i, i);
		b[0] = A(i, DIM - 3); b[1] = A(i, DIM - 2);  b[2] = A(i, DIM - 1);
		c[0] = A(DIM - 3, i); c[1] = A(DIM - 2, i); c[2] = A(DIM - 1, i);
		v_i = vec[i];
		for (auto j = 0; j < 3; ++j) {
			for (auto k = 0; k < 3; ++k) {
				R(j, k) += c[j] * b[k] / a_i;
			}
			v[j] += c[j] * v_i / a_i;
		}
	}
	TVector<MYTYPE> R_sol(3);
	GaussMeth(" ", R, v, R_sol);
	TVector<MYTYPE> sol(DIM);
	for (auto i = 0; i < DIM - 3; ++i) {
		a_i = A(i, i);
		b[0] = A(i, DIM - 3); b[1] = A(i, DIM - 2); b[2] = A(i, DIM - 1);
		v_i = vec[i];
		sol[i] = (v_i - scalar_product(R_sol, b)) / a_i;
	}
	sol[DIM - 3] = R_sol[0];  sol[DIM - 2] = R_sol[1]; sol[DIM - 1] = R_sol[2];
	return sol;
}

template <typename MYTYPE> void BICGSTAB_PRECOND(std::string filename, TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<MYTYPE>& x0, TVector<MYTYPE>& x_real, int power, int precond_type) {
	double TOL = pow(10, (double)power);
	TVector<MYTYPE>(*small_syst_sol)(TMatrix<MYTYPE>&, TVector<MYTYPE>&) = small_syst_sol1;
	if (precond_type == 2)//если система с новым порядком - small_syst_sol2(A, p)
		small_syst_sol = small_syst_sol2;
	TVector<MYTYPE> x0_prec = x0;
	TVector<MYTYPE> r_prev = b - A * x0_prec;
	TVector<MYTYPE> r_next(r_prev);
	TVector<MYTYPE> r0_ast(r_prev); //arbitrary: (r0_ast,r0) != 0
	r0_ast.single_vector();
	TVector<MYTYPE> p(r_prev);
	TVector<MYTYPE> x(x0_prec);
	int iter = 0;
	MYTYPE rho_prev = scalar_product(r_prev, r0_ast);
	clear_file("norm_BiCGStab_precond_" + filename + "_" + std::to_string(power) + ".txt");
	while (iter < MAX_ITER && r_next.norm('3') >= TOL) {
		iter++;
		TVector<MYTYPE> y = small_syst_sol(A, p);
		TVector<MYTYPE> A_mult_y = A * y;
		MYTYPE alpha = rho_prev / scalar_product(A_mult_y, r0_ast);
		TVector<MYTYPE> s = r_prev - alpha * A_mult_y;
		TVector<MYTYPE> z = small_syst_sol(A, s);
		TVector<MYTYPE> A_mult_z = A * z;
		MYTYPE omega = scalar_product(A_mult_z, s) / scalar_product(A_mult_z, A_mult_z);
		x += alpha * y + omega * z;
		r_next = s - omega * A_mult_z;
		filing("norm_BiCGStab_precond_" + filename + "_" + std::to_string(power) + ".txt", (x_real - x).norm('3'));
		MYTYPE rho_next = scalar_product(r_next, r0_ast);
		if (std::abs(rho_next) < EPS_ZERO) {
			std::cout << "\nBREAK! ";
			std::cout << iter << "\n";
			break;
		}
		MYTYPE beta = (alpha * rho_next) / (omega * rho_prev);
		p = r_next + beta * (p - omega * A_mult_y);
		r_prev = r_next;
		rho_prev = rho_next;
	}
	std::cout << "BISGSTAB_precond" << std::to_string(power) << ": iter = " << iter << "\nr_next.norm('3') = " << r_next.norm('3') << "\n";

	clear_file("sol_BiCGStab_precond_" + filename + ".txt");
	filing("sol_BiCGStab_precond_" + filename + ".txt", x);
}

template <typename MYTYPE> void GMRES(std::string filename, TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<MYTYPE>& x0, TVector<MYTYPE>& x, int SUBSP_DIM, TVector<MYTYPE>& x_real) {
	int NEW_SUBSP_DIM = SUBSP_DIM; //if exist j: H(j + 1, j) = 0 => NEW_SUBSP_DIM = j
	TVector<MYTYPE> r0 = b - A * x0;
	MYTYPE beta = r0.norm('2');
	TVector<MYTYPE> gamma_prev = r0;
	TVector<MYTYPE> gamma_next = gamma_prev;
	TVector<MYTYPE> v = 1 / beta * r0;
	TMatrix<MYTYPE> V(x0.get_dim(), SUBSP_DIM);
	V.set_as_column(v, 0);
	TMatrix<MYTYPE> H(SUBSP_DIM + 1, SUBSP_DIM);
	for (auto j = 0; j < SUBSP_DIM; ++j) {
		TVector<MYTYPE> w = A * V.get_column(j);
		for (auto i = 0; i < j + 1; ++i) {
			H(i, j) = scalar_product(w, V.get_column(i));
			w -= H(i, j) * V.get_column(i);
		}
		H(j + 1, j) = w.norm('2');
		if (std::abs(H(j + 1, j)) < EPS_ZERO) {
			NEW_SUBSP_DIM = j + 1;
			break;
		}
		v = 1 / H(j + 1, j) * w;
		if (j + 1 < V.get_ncols()) {
			V.set_as_column(v, j + 1);
		}
		//std::cout << "\nH(j + 1, j) = " << H(j + 1, j) << "\n";
	}
	if (NEW_SUBSP_DIM != SUBSP_DIM) {
		V.remove_last_cols(SUBSP_DIM - NEW_SUBSP_DIM);
		H.remove_last_cols(SUBSP_DIM - NEW_SUBSP_DIM);
		H.remove_last_rows(SUBSP_DIM - NEW_SUBSP_DIM);
		SUBSP_DIM = NEW_SUBSP_DIM;
	}
	//std::cout << "\nV = \n" << V << "\nH = \n" << H << "norms = ";
	//for (auto i = 0; i < V.get_ncols(); ++i) {
	//	std::cout << V.get_column(i).norm('2') <<"; ";
	//}
	//std::cout << "\nscalar products = ";
	//for (auto i = 0; i < V.get_ncols(); ++i) {
	//	for (auto j = i + 1; j < V.get_ncols(); ++j) {
	//		std::cout << scalar_product(V.get_column(i), V.get_column(j)) << "; ";
	//	}
	//}

	TVector<MYTYPE> y = Residual_minimization(beta, H, SUBSP_DIM);
	x = x0 + V * y;
	filing("norm_GMRES_" + filename + ".txt", SUBSP_DIM);
	filing("norm_GMRES_" + filename + ".txt", (x_real - x).norm('3'));
	std::cout << "GMRES " << SUBSP_DIM <<": iter = " << NEW_SUBSP_DIM << "\n(x_real - x).norm('3') = " << (x_real - x).norm('3') << "\n";

	clear_file("sol_GMRES_" + filename + ".txt");
	filing("sol_GMRES_" + filename + ".txt", x);
}

template <typename MYTYPE> void GMRES_PRECOND(std::string filename, TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<MYTYPE>& x0, TVector<MYTYPE>& x, int SUBSP_DIM, TVector<MYTYPE>& x_real, int precond_type) {
	int NEW_SUBSP_DIM = SUBSP_DIM; //if exist j: H(j + 1, j) = 0 => NEW_SUBSP_DIM = j
	TVector<MYTYPE>(*small_syst_sol)(TMatrix<MYTYPE>&, TVector<MYTYPE>&) = small_syst_sol1;
	if (precond_type == 2)//если система с новым порядком - small_syst_sol2(A, p)
		small_syst_sol = small_syst_sol2;
	TVector<MYTYPE> x0_prec = small_syst_sol(A, x0);
	TVector<MYTYPE> r0 = b - A * x0_prec;
	r0 = small_syst_sol(A, r0);
	MYTYPE beta = r0.norm('2');
	TVector<MYTYPE> gamma_prev = r0;
	TVector<MYTYPE> gamma_next = gamma_prev;
	TVector<MYTYPE> v = 1 / beta * r0;
	TMatrix<MYTYPE> V(x0_prec.get_dim(), SUBSP_DIM);
	V.set_as_column(v, 0);
	TMatrix<MYTYPE> H(SUBSP_DIM + 1, SUBSP_DIM);
	for (auto j = 0; j < SUBSP_DIM; ++j) {
		TVector<MYTYPE> w = A * V.get_column(j);
		w = small_syst_sol(A, w);
		for (auto i = 0; i < j + 1; ++i) {
			H(i, j) = scalar_product(w, V.get_column(i));
			w -= H(i, j) * V.get_column(i);
		}
		H(j + 1, j) = w.norm('2');
		if (std::abs(H(j + 1, j)) < EPS_ZERO) {
			NEW_SUBSP_DIM = j + 1;
			break;
		}
		v = 1 / H(j + 1, j) * w;
		if (j + 1 < V.get_ncols()) {
			V.set_as_column(v, j + 1);
		}
		//std::cout << "\nH(j + 1, j) = " << H(j + 1, j) << "\n";
	}
	if (NEW_SUBSP_DIM != SUBSP_DIM) {
		V.remove_last_cols(SUBSP_DIM - NEW_SUBSP_DIM);
		H.remove_last_cols(SUBSP_DIM - NEW_SUBSP_DIM);
		H.remove_last_rows(SUBSP_DIM - NEW_SUBSP_DIM);
		SUBSP_DIM = NEW_SUBSP_DIM;
	}
	//std::cout << "\nV = \n" << V << "\nH = \n" << H << "norms = ";
	//for (auto i = 0; i < V.get_ncols(); ++i) {
	//	std::cout << V.get_column(i).norm('2') <<"; ";
	//}
	//std::cout << "\nscalar products = ";
	//for (auto i = 0; i < V.get_ncols(); ++i) {
	//	for (auto j = i + 1; j < V.get_ncols(); ++j) {
	//		std::cout << scalar_product(V.get_column(i), V.get_column(j)) << "; ";
	//	}
	//}

	TVector<MYTYPE> y = Residual_minimization(beta, H, SUBSP_DIM);
	x = x0_prec + V * y;
	filing("norm_GMRES_precond_" + filename + ".txt", SUBSP_DIM);
	filing("norm_GMRES_precond_" + filename + ".txt", (x_real - x).norm('3'));
	std::cout << "GMRES_precond " << SUBSP_DIM << ": iter = " << NEW_SUBSP_DIM << "\n(x_real - x).norm('3') = " << (x_real - x).norm('3') << "\n";

	clear_file("sol_GMRES_precond_" + filename + ".txt");
	filing("sol_GMRES_precond_" + filename + ".txt", x);
}

template <typename MYTYPE> void Apply_Givens_rotation(TMatrix<MYTYPE>& H, TVector<MYTYPE>& g, int SUBSP_DIM) {
	double temp_i, temp_i_p_1;
	for (auto i = 0; i < SUBSP_DIM; ++i) {
		double znam = sqrt(H(i, i) * H(i, i) + H(i + 1, i) * H(i + 1, i));
		double c = H(i, i) / znam;
		double s = H(i + 1, i) / znam;
		H(i, i) = c * H(i, i) + s * H(i + 1, i);
		H(i + 1, i) = 0.0;
		//Changing two elements of the matrix of the j-th column in i and (i + 1) rows
		for (auto j = i + 1; j < SUBSP_DIM; ++j) {
			temp_i = c * H(i, j) + s * H(i + 1, j);
			temp_i_p_1 = -s * H(i, j) + c * H(i + 1, j);
			H(i, j) = temp_i;
			H(i + 1, j) = temp_i_p_1;
		}
		//Change of the i-th and (i + 1)-th elements of the vector
		temp_i = c * g[i];
		temp_i_p_1 = -s * g[i];
		g[i] = temp_i;
		g[i + 1] = temp_i_p_1;
	}
}

//Minimize ||beta * e1 - H * y||_2
template <typename MYTYPE> TVector<MYTYPE> Residual_minimization(MYTYPE beta, TMatrix<MYTYPE>& H, int SUBSP_DIM) {
	TVector<MYTYPE> g(H.get_nrows());
	g[0] = beta;
	//Givens_rotation
	Apply_Givens_rotation(H, g, SUBSP_DIM);
	//Here H = R_m and g = g_m

	//std::cout << "\nH = \n" << H << "\ng = \n" << g;
	H.remove_last_rows(1);
	g.remove_last_elems(1);

	TVector<MYTYPE> y_m(SUBSP_DIM);
	//Reverse Gauss

	for (auto i = SUBSP_DIM - 1; i >= 0; --i) {
		double sum = 0.0;
		for (auto j = SUBSP_DIM - 1; j >= i + 1; --j) {
			sum += H(i, j) * y_m[j];
		}
		y_m[i] = (g[i] - sum) / H(i, i);
	}

	return y_m;
}

template <typename MYTYPE> TVector<MYTYPE> GMRES_m(TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<MYTYPE>& x0, int SUBSP_DIM) {
	const double EPS = 1e-5;
	int rest_num = 0;
	TVector<MYTYPE> x = GMRES(A, b, x0, SUBSP_DIM);
	TVector<MYTYPE> r = b - A * x;
	while (r.norm('2') >= EPS && rest_num <= MAX_ITER) {
		rest_num++;
		//std::cout << "RESTART!\n";
		GMRES(A, b, x, SUBSP_DIM);
		r = b - A * x;
	}
	//std::cout << "rest_num = " << rest_num << "\n";
	return x;
}

//BICGSTAB
	//TVector<MYTYPE> r_prev = b - A * x0;
	//TVector<MYTYPE> r_next(r_prev);
	//TVector<MYTYPE> r0_ast(r_prev); //arbitrary: (r0_ast,r0) != 0
	//MYTYPE rho_prev = MYTYPE(1);
	//MYTYPE alpha = MYTYPE(1);
	//MYTYPE omega_prev = MYTYPE(1);
	//TVector<MYTYPE> p_prev(x0.get_dim());
	//TVector<MYTYPE> v_prev(x0.get_dim());
	//TVector<MYTYPE> x_prev(x0);
	//TVector<MYTYPE> x_next(x0);
	//int iter = 0;
	//clear_file("norm1.txt");
	//while (iter < MAX_ITER && r_next.norm('2') >= TOL) {
	//	iter++;
	//	MYTYPE rho_next = scalar_product(r_prev, r0_ast);
	//	MYTYPE beta = (alpha * rho_next) / (omega_prev * rho_prev);
	//	TVector<MYTYPE> p_next = r_prev + beta * (p_prev - omega_prev * v_prev);
	//	TVector<MYTYPE> v_next = A * p_next;
	//	alpha = rho_next / scalar_product(r0_ast, v_next);
	//	TVector<MYTYPE> s = r_prev - alpha * v_next;
	//	TVector<MYTYPE> t = A * s;
	//	MYTYPE omega_next = scalar_product(t, s) / scalar_product(t, t);
	//	x_next = x_prev + alpha * p_next + omega_next * s;
	//	r_next = s - omega_next * t;
	//	filing("norm1.txt", r_next.norm('2'));
	//	rho_prev = rho_next;
	//	p_prev = p_next;
	//	v_prev = v_next;
	//	x_prev = x_next;
	//	omega_prev = omega_next;
	//	r_prev = r_next;
	//}
	//std::cout << "iter = " << iter << "\nr_next.norm('2') = " << r_next.norm('2') << "\n";
	//filing("newSOL.txt", x_next);
	//return x_next;

//template <typename MYTYPE> void DQGMRES(TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<MYTYPE>& x0, TVector<MYTYPE>& x, int SUBSP_DIM, TVector<MYTYPE>& x_real) {
//		int NEW_SUBSP_DIM = SUBSP_DIM; //if exist j: H(j + 1, j) = 0 => NEW_SUBSP_DIM = j
//		TVector<MYTYPE> r0 = b - A * x0;
//		MYTYPE beta = r0.norm('2');
//		MYTYPE gamma_prev = r0.norm('2');
//		MYTYPE gamma_next = gamma_prev;
//		TVector<MYTYPE> v = 1 / gamma_prev * r0;
//		TMatrix<MYTYPE> V(x0.get_dim(), SUBSP_DIM);
//		V.set_as_column(v, 0);
//		TMatrix<MYTYPE> H(SUBSP_DIM + 1, SUBSP_DIM);
//		for (auto j = 0; j < SUBSP_DIM; ++j) {
//			TVector<MYTYPE> w = A * V.get_column(j);
//			for (auto i = 0; i < j + 1; ++i) {
//				H(i, j) = scalar_product(w, V.get_column(i));
//				w -= H(i, j) * V.get_column(i);
//			}
//			//H(j + 1, j) = w.norm('2');
//			//v = 1 / H(j + 1, j) * w;
//			//if (j + 1 < V.get_ncols()) {
//			//	V.set_as_column(v, j + 1);
//			//}
//
//			Apply_Givens_rotation_m(H, j +1);
//
//			double znam = sqrt(H(j, j) * H(j, j) + H(j + 1, j) * H(j + 1, j));
//			double c = H(j, j) / znam;
//			double s = H(j + 1, j) / znam;
//			gamma_next = -s * gamma_prev;
//			gamma_prev = c * gamma_prev;
//			H(j, j) = c * H(j, j) + s * H(j + 1, j);
//
//			if (std::abs(gamma_next) < EPS_ZERO) {
//				NEW_SUBSP_DIM = j + 1;
//				break;
//			}
//		}
//		if (NEW_SUBSP_DIM != SUBSP_DIM) {
//			V.remove_last_cols(SUBSP_DIM - NEW_SUBSP_DIM);
//			H.remove_last_cols(SUBSP_DIM - NEW_SUBSP_DIM);
//			H.remove_last_rows(SUBSP_DIM - NEW_SUBSP_DIM);
//			SUBSP_DIM = NEW_SUBSP_DIM;
//		}
//
//		TVector<MYTYPE> y = Residual_minimization(gamma_next, H, SUBSP_DIM);
//		x = x0 + V * y;
//		std::cout << "iter = " << NEW_SUBSP_DIM << "\n(x_real - x).norm('2') = " << (x_real - x).norm('2') << "\n";
//
//		filing("GMRES_sol.txt", x);
//	}
