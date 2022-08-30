#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
#include "omp.h"

int N = 2048;

void separation_LU(double* matr, double* L, double* U);
double max_elem(double* matr);
void sub_matr(double* result_matr, double* source_matr);

void init_trig_matr(double*& matr1, double*& matr2);
void init_system(std::string filename, double*& matr, double*& vec);
void init_matr(std::string filename, double*& matr);
void init_vec(std::string filename, double*& vec);
void write_rand_matr(std::string filename);
void read_matr(std::string filename, double*& matr);
void read_vec(std::string filename, double*& vec);
void init_zero_matr(double*& matr);
void clear_mem(double* obj);
void print_matr(double* matr);
void print_vec(double* vec);

void LU_decomposition(double* matr, const int N, const int proc_num = 1);
void LU_block_decomposition(double* matr, const int N, const int bl_size, const int proc_num = 1);
void LU_decomposition_WIKI(double* matr, double* L, double* U);

void matr_product(double* A, double* B, double* C, int nth = 1);
void matr_product_transp(double* A, double* B_transp, double* C);
void matr_product_transpOMP(double* A, double* B_transp, double* C, int nth = omp_get_max_threads());

void reverse_Gauss_L(double* matr, double* vec);
void reverse_Gauss_U(double* matr, double* vec);
void SLAE_solve(double* A, double* b);

int main() {

	//------------------РЕШЕНИЕ СИСТЕМ------------------
	//double * A = NULL, * b = NULL;
	//init_system("test5.txt", A, b);
	//std::cout << "A:\n"; print_matr(A);
	//LU_decomposition(A, N);
	//std::cout << "LU:\n"; print_matr(A);
	//read_matr("test5_matr.txt", A);
	//SLAE_solve(A, b);
	//std::cout << "LU:\n"; print_matr(A);
	//std::cout << "Sol:\n"; print_vec(b);

	//clear_mem(A);
	//clear_mem(b);

	//std::cout << "N = " << N << '\n';

	//------------------LU-РАЗЛОЖЕНИЕ (ПРЯМОЕ И БЛОЧНОЕ)------------------
	// write_rand_matr("rand_matr" + std::to_string(N) + ".txt");
	std::cout << "\nN = " << N << "\n";
	double t1, t2;
	double* A = NULL;
	std::string filename = "rand_matr2048.txt";//"rand_matr" + std::to_string(N) + ".txt"
	init_matr(filename, A);
	t1 = omp_get_wtime();
	LU_decomposition(A, N);
	t2 = omp_get_wtime();
	//std::cout << "A:\n"; print_matr(A);
	double time_seq = t2 - t1;
	std::cout << "time_seq = " << time_seq << " sec\n";
	//std::cout << "A[N/2 + 1][N/2] = " << A[(N / 2 + 1) * N + N / 2] << "\n\n";

	int nth = 4;

	read_matr(filename, A);
	t1 = omp_get_wtime();
	LU_decomposition(A, N, nth);
	t2 = omp_get_wtime();
	//std::cout << "A:\n"; print_matr(A);
	double time_parall = t2 - t1;
	std::cout << "time_parall = " << time_parall << " sec\n";
	//std::cout << "A[N/2 + 1][N/2] = " << A[(N / 2 + 1) * N + N / 2] << "\n";
	std::cout << "time_seq / time_parall = " << time_seq / time_parall << "\n\n";

	double best_time_seq = 1000.0, best_time_parall = 1000.0;
	int best_bs_seq = 0, best_bs_parall = 0;
	for (int _bs = 4; _bs <= 512; _bs *= 2) {
		read_matr(filename, A);
		t1 = omp_get_wtime();
		LU_block_decomposition(A, N, _bs);
		t2 = omp_get_wtime();
		//std::cout << "A:\n"; print_matr(A);
		double time_block = t2 - t1;
		std::cout << "bs = " << _bs << "\ntime_block = " << time_block << " sec\n";
		//std::cout << "A[N/2 + 1][N/2] = " << A[(N / 2 + 1) * N + N / 2] << "\n";
		std::cout << "time_seq / time_block = " << time_seq / time_block << "\n";

		read_matr(filename, A);
		t1 = omp_get_wtime();
		LU_block_decomposition(A, N, _bs, nth);
		t2 = omp_get_wtime();
		//std::cout << "A:\n"; print_matr(A);
		double time_block_parall = t2 - t1;
		std::cout << "time_block_parall = " << time_block_parall << " sec\n";
		//std::cout << "A[N/2 + 1][N/2] = " << A[(N / 2 + 1) * N + N / 2] << "\n";
		std::cout << "time_block / time_block_parall = " << time_block / time_block_parall << "\n\n";
		if (time_block < best_time_seq) {
			best_time_seq = time_block;
			best_bs_seq = _bs;
		}
		if (time_block_parall < best_time_parall) {
			best_time_parall = time_block_parall;
			best_bs_parall = _bs;
		}
	}
	printf("best_time_seq = %f, best_bs_seq = %d\nbest_time_parall = %f, best_bs_parall = %d\naccelerateon = %f\n\n",
		   best_time_seq, best_bs_seq, best_time_parall, best_bs_parall, best_time_seq / best_time_parall);

	//double* L = NULL;
	//double* U = NULL;
	//double* res = NULL;
	//init_zero_matr(L);
	//init_zero_matr(U);
	//init_zero_matr(res);
	//separation_LU(A, L, U);
	//matr_product(L, U, res);
	//read_matr(filename, A);
	//sub_matr(res, A);
	//std::cout << "max elem = " << max_elem(res);

	//int _bs = 64;
	//read_matr(filename, A);
	//t1 = omp_get_wtime();
	//LU_block_decomposition(A, N, _bs);
	//t2 = omp_get_wtime();
	////std::cout << "A:\n"; print_matr(A);
	//double time_block = t2 - t1;
	//std::cout << "bs = " << _bs << ", time (block) = " << time_block << " sec\n";
	//std::cout << "A[N/2 + 1][N/2] = " << A[(N / 2 + 1) * N + N / 2] << "\n";
	//std::cout << "time_seq / time_block = " << time_seq / time_block << "\n\n";

	//_bs = 128;
	//read_matr(filename, A);
	//t1 = omp_get_wtime();
	//LU_block_decomposition(A, N, _bs, nth);
	//t2 = omp_get_wtime();
	////std::cout << "A:\n"; print_matr(A);
	//double time_block_parall = t2 - t1;
	//std::cout << "bs = " << _bs << ", time (block_parall) = " << time_block_parall << " sec\n";
	//std::cout << "A[N/2 + 1][N/2] = " << A[(N / 2 + 1) * N + N / 2] << "\n";
	//std::cout << "time_block / time_block_parall = " << time_block / time_block_parall << "\n\n";


	clear_mem(A);
	//clear_mem(L);
	//clear_mem(U);
	//clear_mem(res);

	////------------------МАТРИЧНОЕ УМНОЖЕНИЕ (ПРЯМОЕ И OMP)------------------
	//N = 2048;
	//double * A = NULL, * B = NULL, * C = NULL;
	//double t1, t2;
	//
	//init_trig_matr(A, B);
	////init_matr("rand_matr" + std::to_string(N) + "A.txt", A);
	////init_matr("rand_matr" + std::to_string(N) + "B.txt", B);
	//init_zero_matr(C);

	//t1 = omp_get_wtime();
	//matr_product_transp(A, B, C);
	//t2 = omp_get_wtime();
	//double time_seq = t2 - t1;
	//std::cout << "Time (seq) = " << time_seq << " sec\n";
	////std::cout << "C:\n"; print_matr(C);

	//for (int p = 1; p <= 4; ++p) {
	//	t1 = omp_get_wtime();
	//	matr_product_transpOMP(A, B, C, p);
	//	t2 = omp_get_wtime();
	//	double time_omp = t2 - t1;
	//	std::cout << "Time (OMP) = " << time_omp << " sec\n";
	//	//std::cout << "C:\n"; print_matr(C);
	//	std::cout << "time_seq / time_omp = " << time_seq / time_omp << "\n\n";
	//}

	//clear_mem(A);
	//clear_mem(B);
	//clear_mem(C);

	//------------------ЗАПИСЬ РАНДОМНЫХ МАТРИЦ------------------
	//write_rand_matr("rand_matr" + std::to_string(N) + "A.txt");
	//write_rand_matr("rand_matr" + std::to_string(N) + "B.txt");

	return 0;
}

void separation_LU(double* matr, double* L, double* U) {
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j <= i; ++j) {
			L[i * N + j] = matr[i * N + j];
			U[j * N + i] = matr[j * N + i];
		}
		L[i * N + i] = 1.0;
	}
}

double max_elem(double* matr) {
	double max = 0.0;
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			double temp = fabs(matr[i * N + j]);
			if (max < temp) {
#pragma omp critical
				{
					if (max < temp) {
						max = temp;
					}
				}
			}
		}
	}
	return max;
}

void sub_matr(double* result_matr, double* source_matr) {
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			result_matr[i * N + j] -= source_matr[i * N + j];
		}
	}
}

void init_system(std::string filename, double*& matr, double*& vec)
{
	std::ifstream fin;
	fin.open(filename.c_str(), std::ios::in);
	if (!fin.is_open()) {
		std::cout << "Error file!\n";
	}
	else {
		fin >> N;
		matr = new double[N * N];
		vec = new double[N];
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				fin >> matr[i * N + j];
			}
		}
		for (int i = 0; i < N; ++i) {
			fin >> vec[i];
		}
	}
	fin.close();
}

void init_matr(std::string filename, double*& matr)
{
	std::ifstream fin;
	fin.open(filename.c_str(), std::ios::in);
	if (!fin.is_open()) {
		std::cout << "Error file!\n";
	}
	else {
		fin >> N;
		matr = new double[N * N];
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				fin >> matr[i * N + j];
			}
		}
	}
	fin.close();
}

void init_vec(std::string filename, double*& vec)
{
	std::ifstream fin;
	fin.open(filename.c_str(), std::ios::in);
	if (!fin.is_open()) {
		std::cout << "Error file!\n";
	}
	else {
		fin >> N;
		vec = new double[N];
		for (int i = 0; i < N; ++i) {
			fin >> vec[i];
		}
	}
	fin.close();
}

//без выделения памяти (повторное считывание)
void read_matr(std::string filename, double*& matr) {
	std::ifstream fin;
	fin.open(filename.c_str(), std::ios::in);
	if (!fin.is_open()) {
		std::cout << "Error file!\n";
	}
	else {
		fin >> N;
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				fin >> matr[i * N + j];
			}
		}
	}
	fin.close();
}

void write_rand_matr(std::string filename)
{
	std::ofstream fout;
	fout.open(filename.c_str());
	if (!fout.is_open()) {
		std::cout << "Error file!\n";
	}
	else {
		fout << N << '\n';
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				fout << rand() % N + 1.0 << '\t';
			}
			fout << '\n';
		}
	}
	fout.close();
}

void read_vec(std::string filename, double*& vec) {
	std::ifstream fin;
	fin.open(filename.c_str(), std::ios::in);
	if (!fin.is_open()) {
		std::cout << "Error file!\n";
	}
	else {
		fin >> N;
		for (int i = 0; i < N; ++i) {
			fin >> vec[i];
		}
	}
	fin.close();
}

void init_trig_matr(double*& matr1, double*& matr2) {
	matr1 = new double[N * N];
	matr2 = new double[N * N];
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			matr1[i * N + j] = cos(i - j);
			matr2[i * N + j] = sin(i + j);
		}
	}
}

void init_zero_matr(double*& matr) {
	matr = new double[N * N];
	for (int i = 0; i < N * N; ++i) {
		matr[i] = 0.0;
	}
}

double norm_matr(double*& matr, int N) {
	double norm = 0.0;
	for (int i = 0; i < N; ++i) {
		double sum = 0.0;
		for (int j = 0; j < N; ++j) {
			sum += fabs(matr[i * N + j]);
		}
		if (sum > norm) {
			norm = sum;
		}
	}
	return norm;
}

void clear_mem(double* obj) {
	delete[] obj;
}

void print_matr(double* matr) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			std::cout << std::setprecision(8) << std::fixed << matr[i * N + j] << '\t';
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

void print_vec(double* vec) {
	for (int i = 0; i < N; ++i) {
		std::cout << std::setprecision(8) << std::fixed << vec[i] << '\t';
	}
	std::cout << '\n';
}

void LU_decomposition(double* matr, const int N, const int proc_num) {
	for (int i = 0; i < N - 1; ++i) {
#pragma omp parallel for num_threads(proc_num) if(proc_num > 1)
		for (int j = i + 1; j < N; ++j) {
			double temp = matr[j * N + i] / matr[i * N + i];
			for (int k = i + 1; k < N; ++k) {
				matr[j * N + k] -= temp * matr[i * N + k];
			}
			matr[j * N + i] = temp;
		}
	}
}

void LU_block_decomposition(double* matr, const int N, const int bl_size, const int proc_num) {
	double* A11 = new double[bl_size * bl_size];
	double* A12 = new double[bl_size * (N - bl_size)];//храним по столбцам
	double* A21 = new double[(N - bl_size) * bl_size];//храним по строкам
	for (int bi = 0; bi < N - 1; bi += bl_size) {
		//заполняем блоки//подумать, как сделать это оптимальн
#pragma omp parallel for num_threads(proc_num) if((proc_num > 1) && (bl_size * bl_size > 4000))
		for (int i = 0; i < bl_size; ++i) {
			for (int j = 0; j < bl_size; ++j) {
				A11[i * bl_size + j] = matr[(i + bi) * N + (j + bi)];
			}
		}
#pragma omp parallel for num_threads(proc_num) if((proc_num > 1) && (bl_size * (N - bi - bl_size) > 4000))
		for (int j = 0; j < (N - bi - bl_size); ++j) {
			for (int i = 0; i < bl_size; ++i) {
				A12[j * bl_size + i] = matr[(i + bi) * N + (j + bi + bl_size)];
			}
		}
#pragma omp parallel for num_threads(proc_num) if((proc_num > 1) && (bl_size * (N - bi - bl_size) > 4000))
		for (int i = 0; i < (N - bi - bl_size); ++i) {
			for (int j = 0; j < bl_size; ++j) {
				A21[i * bl_size + j] = matr[(i + bi + bl_size) * N + (j + bi)];
			}
		}
		//обычное LU-разложение для A11
		LU_decomposition(A11, bl_size, proc_num);
		//ищем U12: L11.U12 = A12
		double sum;
		//#pragma omp parallel for num_threads(proc_num) private(sum) if((proc_num > 1) && (bl_size * bl_size * (N - bi - bl_size) > 4000))
#pragma omp parallel for num_threads(proc_num) private(sum) if((proc_num > 1) && ((bl_size + 1) * bl_size / 2 * (N - bi - bl_size) > 4000))
		for (int p = 0; p < (N - bi - bl_size); ++p) {//цикл по количеству решаемых систем (они же - количество столбцов U12)
			for (int j = 1; j < bl_size; ++j) {
				sum = 0.0;
				for (int k = 0; k <= j - 1; ++k) {
					sum += A11[j * bl_size + k] * A12[p * bl_size + k];
				}
				A12[p * bl_size + j] -= sum;
			}
		}
		//ищем L21: L21.U11 = A21
//#pragma omp parallel for num_threads(proc_num) private(sum) if((proc_num > 1) && (bl_size * bl_size * (N - bi - bl_size) > 4000))
#pragma omp parallel for num_threads(proc_num) private(sum) if((proc_num > 1) && ((bl_size + 1) * bl_size / 2 * (N - bi - bl_size) > 4000))
		for (int p = 0; p < (N - bi - bl_size); ++p) {//цикл по количеству решаемых систем (они же - количество строк L21)
			for (int j = 0; j < bl_size; ++j) {
				sum = 0.0;
				for (int k = 0; k <= j - 1; ++k) {
					sum += A21[p * bl_size + k] * A11[k * bl_size + j];
				}
				A21[p * bl_size + j] = (A21[p * bl_size + j] - sum) / A11[j * bl_size + j];
			}
		}
		//ищем \tilde{A22} = A22 - L21.U12
#pragma omp parallel for num_threads(proc_num) private(sum) if((proc_num > 1) && (bl_size * (N - bi - bl_size) * (N - bi - bl_size) > 4000))
		for (int i = 0; i < (N - bi - bl_size); ++i) {
			for (int j = 0; j < (N - bi - bl_size); ++j) {
				sum = 0.0;
				for (int k = 0; k < bl_size; ++k) {
					sum += A21[i * bl_size + k] * A12[j * bl_size + k];
				}
				matr[(i + bi + bl_size) * N + (j + bi + bl_size)] -= sum;
			}
		}
		//запись блоков L11U11, L21 и U12 в исходную матрицу
#pragma omp parallel for num_threads(proc_num) if((proc_num > 1) && (bl_size * bl_size > 4000))
		for (int i = 0; i < bl_size; ++i) {
			for (int j = 0; j < bl_size; ++j) {
				matr[(i + bi) * N + (j + bi)] = A11[i * bl_size + j];
			}
		}
#pragma omp parallel for num_threads(proc_num) if((proc_num > 1) && (bl_size * (N - bi - bl_size) > 4000))
		for (int j = 0; j < (N - bi - bl_size); ++j) {
			for (int i = 0; i < bl_size; ++i) {
				matr[(i + bi) * N + (j + bi + bl_size)] = A12[j * bl_size + i];
			}
		}
#pragma omp parallel for num_threads(proc_num) if((proc_num > 1) && (bl_size * (N - bi - bl_size) > 4000))
		for (int i = 0; i < (N - bi - bl_size); ++i) {
			for (int j = 0; j < bl_size; ++j) {
				matr[(i + bi + bl_size) * N + (j + bi)] = A21[i * bl_size + j];
			}
		}
	}
	clear_mem(A12);
	clear_mem(A11);
	clear_mem(A21);
}

void LU_decomposition_WIKI(double* matr, double* L, double* U) {
	for (int j = 0; j < N; ++j) {
		U[0 * N + j] = matr[0 * N + j];
		L[j * N + 0] = matr[j * N + 0] / U[0 * N + 0];
	}
	for (int i = 1; i < N; ++i) {
		for (int j = i; j < N; ++j) {
			double temp1 = 0.0;
			double temp2 = 0.0;
			for (int k = 0; k < i; ++k) {
				temp1 += L[i * N + k] * U[k * N + j];
				temp2 += L[j * N + k] * U[k * N + i];
			}
			U[i * N + j] = matr[i * N + j] - temp1;
			L[j * N + i] = (matr[j * N + i] - temp2) / U[i * N + i];
		}
	}
}

void reverse_Gauss_L(double* matr, double* vec/*, double* sol*/) {
	double sum;
	for (int i = 0; i < N; ++i) {
		sum = 0.0;
		for (int j = 0; j <= i - 1; ++j) {
			sum += matr[i * N + j] * vec[j]; //sum += matr[i * N + j] * sol[j];
		}
		vec[i] = (vec[i] - sum); //sol[i] = (vec[i] - sum) / matr[i * N + i];
	}
}

void reverse_Gauss_U(double* matr, double* vec/*, double* sol*/) {
	double sum;
	for (int i = N - 1; i >= 0; --i) {
		sum = 0.0;
		for (int j = N - 1; j >= i + 1; --j) {
			sum += matr[i * N + j] * vec[j]; //sum += matr[i * N + j] * sol[j];
		}
		vec[i] = (vec[i] - sum) / matr[i * N + i]; //sol[i] = (vec[i] - sum) / matr[i * N + i];
	}
}

void SLAE_solve(double* A, double* b) {
	LU_decomposition(A, N);
	reverse_Gauss_L(A, b);
	reverse_Gauss_U(A, b);
}

void matr_product(double* A, double* B, double* C, int nth) {
#pragma omp parallel for num_threads(nth) if(nth > 1)
	for (int i = 0; i < N; ++i) {
		for (int k = 0; k < N; ++k) {
			for (int j = 0; j < N; ++j) {
				C[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}
}

void matr_product_transp(double* A, double* B_transp, double* C) {
	double sum = 0.0;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			sum = 0.0;
			for (int k = 0; k < N; ++k) {
				sum += A[i * N + k] * B_transp[j * N + k];
			}
			C[i * N + j] = sum;
		}
	}
	std::cout << "C[N/2 + 1][N/2] = " << C[(N / 2 + 1) * N + N / 2] << '\n';
}

void matr_product_transpOMP(double* A, double* B_transp, double* C, int nth) {
	omp_set_num_threads(nth);
	std::cout << "nth = " << nth << '\n';
#pragma omp parallel
	{
		//std::cout << "nth = " << omp_get_num_threads() << '\n';
		int th = omp_get_thread_num();
		//printf("th = %d\n", th);
		for (int i = th; i < N; i += nth) {
			for (int j = 0; j < N; ++j) {
				double sum = 0.0;
				for (int k = 0; k < N; ++k) {
					sum += A[i * N + k] * B_transp[j * N + k];
				}
				C[i * N + j] = sum;
			}
			//printf("i = %d\n", i);
		}
	}
	std::cout << "C[N/2 + 1][N/2] = " << C[(N / 2 + 1) * N + N / 2] << '\n';
}
