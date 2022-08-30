#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cassert>
#include "omp.h"

int N = 2048;
int N_THREADS = 4;


void read_matr(std::string filename, double*& matr) {
	std::ifstream fin;
	fin.open(filename.c_str(), std::ios::in);
	if (!fin.is_open()) {
		std::cout << "Error file!\n";
		return;
	}
	fin >> N;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			fin >> matr[i * N + j];
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
		return;
	}
	fout << N << '\n';
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			fout << rand() % N + 1.0 << '\t';
		}
		fout << '\n';
	}
	fout.close();
}

void print_matr(double* matr, const int N) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			std::cout << std::setprecision(8) << std::fixed << matr[i * N + j] << '\t';
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

void separation_LU(double* matr, double* L, double* U, const int N) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j <= i; ++j) {
			L[i * N + j] = matr[i * N + j];
			U[j * N + i] = matr[j * N + i];
		}
		L[i * N + i] = 1.0;
	}
}

void LU_decomposition(double* matr, const int N, const int proc_num = 1) {
	for (int i = 0; i < N - 1; ++i) {
#pragma omp parallel for num_threads(proc_num)
		for (int j = i + 1; j < N; ++j) {
			double temp = matr[j * N + i] / matr[i * N + i];
			for (int k = i + 1; k < N; ++k) {
				matr[j * N + k] -= temp * matr[i * N + k];
			}
			matr[j * N + i] = temp;
		}
	}
}

void LU_block_decomposition(double* matr, const int N, const int bl_size, const int proc_num = 1) {
	double* A11 = new double[bl_size * bl_size];
	double* A12 = new double[bl_size * (N - bl_size)];
	double* A21 = new double[(N - bl_size) * bl_size];

	for (int bi = 0; bi < N - 1; bi += bl_size) {
#pragma omp parallel for num_threads(proc_num)
		for (int i = 0; i < bl_size; ++i) {
			for (int j = 0; j < bl_size; ++j) {
				A11[i * bl_size + j] = matr[(i + bi) * N + (j + bi)];
			}
		}

#pragma omp parallel for num_threads(proc_num)
		for (int j = 0; j < (N - bi - bl_size); ++j) {
			for (int i = 0; i < bl_size; ++i) {
				A12[j * bl_size + i] = matr[(i + bi) * N + (j + bi + bl_size)];
			}
		}

#pragma omp parallel for num_threads(proc_num)
		for (int i = 0; i < (N - bi - bl_size); ++i) {
			for (int j = 0; j < bl_size; ++j) {
				A21[i * bl_size + j] = matr[(i + bi + bl_size) * N + (j + bi)];
			}
		}

		LU_decomposition(A11, bl_size, proc_num);
		double sum;

#pragma omp parallel for num_threads(proc_num) private(sum)
		for (int p = 0; p < (N - bi - bl_size); ++p) {
			for (int j = 1; j < bl_size; ++j) {
				sum = 0.0;
				for (int k = 0; k <= j - 1; ++k) {
					sum += A11[j * bl_size + k] * A12[p * bl_size + k];
				}
				A12[p * bl_size + j] -= sum;
			}
		}

#pragma omp parallel for num_threads(proc_num) private(sum)
		for (int p = 0; p < (N - bi - bl_size); ++p) {
			for (int j = 0; j < bl_size; ++j) {
				sum = 0.0;
				for (int k = 0; k <= j - 1; ++k) {
					sum += A21[p * bl_size + k] * A11[k * bl_size + j];
				}
				A21[p * bl_size + j] = (A21[p * bl_size + j] - sum) / A11[j * bl_size + j];
			}
		}

#pragma omp parallel for num_threads(proc_num) private(sum)
		for (int i = 0; i < (N - bi - bl_size); ++i) {
			for (int j = 0; j < (N - bi - bl_size); ++j) {
				sum = 0.0;
				for (int k = 0; k < bl_size; ++k) {
					sum += A21[i * bl_size + k] * A12[j * bl_size + k];
				}
				matr[(i + bi + bl_size) * N + (j + bi + bl_size)] -= sum;
			}
		}

#pragma omp parallel for num_threads(proc_num)
		for (int i = 0; i < bl_size; ++i) {
			for (int j = 0; j < bl_size; ++j) {
				matr[(i + bi) * N + (j + bi)] = A11[i * bl_size + j];
			}
		}

#pragma omp parallel for num_threads(proc_num)
		for (int j = 0; j < (N - bi - bl_size); ++j) {
			for (int i = 0; i < bl_size; ++i) {
				matr[(i + bi) * N + (j + bi + bl_size)] = A12[j * bl_size + i];
			}
		}

#pragma omp parallel for num_threads(proc_num)
		for (int i = 0; i < (N - bi - bl_size); ++i) {
			for (int j = 0; j < bl_size; ++j) {
				matr[(i + bi + bl_size) * N + (j + bi)] = A21[i * bl_size + j];
			}
		}
	}

	delete[] A12;
	delete[] A11;
	delete[] A21;
}

int main(int argc, char const *argv[]) {
	if (argc >= 2) {
		N = atoi(argv[1]);
	}
	N_THREADS = omp_get_max_threads();
	
	std::cout << "N = " << N << "\n";
	std::cout << "num_threads: " << N_THREADS << "\n\n";

	double* A = new double[N * N];
	double t1, t2;
	int check_idx = (N / 2 + 1) * N + N / 2;
	double check_val;
	std::string filename = "rand_matr" + std::to_string(N) + ".txt";
	write_rand_matr(filename);

	read_matr(filename, A);
	t1 = omp_get_wtime();
	LU_decomposition(A, N);
	t2 = omp_get_wtime();
	double time_seq = t2 - t1;
	std::cout << "time_seq = " << time_seq << " sec\n";
	check_val = A[check_idx];

	read_matr(filename, A);
	t1 = omp_get_wtime();
	LU_decomposition(A, N, N_THREADS);
	t2 = omp_get_wtime();
	double time_parall = t2 - t1;
	std::cout << "time_parall = " << time_parall << " sec\n";
	assert(A[check_idx] == check_val);
	std::cout << "time_seq / time_parall = " << time_seq / time_parall << "\n\n";

	double best_time_seq = 1000.0, best_time_parall = 1000.0;
	int best_bs_seq = 0, best_bs_parall = 0;
	for (int bs = 4; bs <= 512; bs *= 2) {
		read_matr(filename, A);
		t1 = omp_get_wtime();
		LU_block_decomposition(A, N, bs);
		t2 = omp_get_wtime();
		double time_block = t2 - t1;
		check_val = A[check_idx];

		read_matr(filename, A);
		t1 = omp_get_wtime();
		LU_block_decomposition(A, N, bs, N_THREADS);
		t2 = omp_get_wtime();
		double time_block_parall = t2 - t1;
		assert(A[check_idx] == check_val);

		if (time_block < best_time_seq) {
			best_time_seq = time_block;
			best_bs_seq = bs;
		}
		if (time_block_parall < best_time_parall) {
			best_time_parall = time_block_parall;
			best_bs_parall = bs;
		}
	}
	printf("best_time_seq = %f, best_bs_seq = %d\nbest_time_parall = %f, best_bs_parall = %d\nacceleration = %f\n\n",
		   best_time_seq, best_bs_seq, best_time_parall, best_bs_parall, best_time_seq / best_time_parall);

	delete[] A;
	return 0;
}
