#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <omp.h>

const int N = 2048;

void init_matr(double*& matr);          //создание
void comp_matr(double* matr);           //рандомные числа
void comp_matr_zero(double* matr);      //нулевые элементы
void clear_mem(double* obj);            //отчистка памяти
void print_matr(double* matr);          //печать матрицы
void print_element(double* matr);       //печать элемента matr[n/2][n/2]
void LU_decomposition(double* matr);											//обычная LU
void LU_decomposition_parallel(double* matr);									//параллельная LU
void block_LU_decomposition(double* matr, const int bs);						//блочная LU
void block_LU_decomposition_parallel(double* matr, const int bs, int num_th);	//блочная параллельная LU
void block_LU_decomposition_un(double* matr, const int bs, int num_th);
void matr_product(double* A, double* B, double* C);								//умножение матриц
void separation_LU(double* matr, double* L, double* U);
void matr_equal(double* result_matr, double* source_matr);
void sub_matr(double* result_matr, double* source_matr);
double max_matr(double* matr);
void read_matr(std::string namefile, double*& matr);  //запись из файла
void write_matr(std::string namefile, double*& matr); //запись в файл


int main() {
	double* matr = nullptr, * matr_copy = nullptr, * L = nullptr, * U = nullptr;

	init_matr(matr);
	init_matr(matr_copy);
	init_matr(L);
	init_matr(U);

	int flag = 1;
	std::string file_name = "rand_matr2048.txt";

	comp_matr(matr);
	// write_matr(file_name, matr); //если  матрица есть, то не нужно

	read_matr(file_name, matr);
	// read_matr(file_name, matr_copy);

	double t1, t2;
	double max;

	t1 = omp_get_wtime();
	LU_decomposition(matr);
	t2 = omp_get_wtime();
	std::cout << "LU, time = " << (double)(t2 - t1) << "\n";
	print_element(matr);
	double t_LU = double(t2 - t1);

	read_matr(file_name, matr);

	t1 = omp_get_wtime();
	LU_decomposition_parallel(matr);
	t2 = omp_get_wtime();
	std::cout << "LU parallel, time = " << (double)(t2 - t1) << "\n";
	print_element(matr);
	double t_LU_par = double(t2 - t1);


	std::cout << "time seq / time par = " << t_LU / t_LU_par << "\n";

	//separation_LU(matr, L, U);
	//comp_matr_zero(matr);
	//matr_product(L, U, matr);
	//sub_matr(matr, matr_copy);
	//max = max_matr(matr);
	//std::cout << "max (matr - L * U) = " << max << "\n";

	int best_bs = 0, best_bs_parallel = 0;
	double best_time = 1000.0, best_time_parallel = 1000.0;

	for (int bs = 16; bs <= 65; bs *= 2) {
		//int bs = 64;
		std::cout << "----------------------------\nbs = " << bs << std::endl;
		read_matr(file_name, matr);
		t1 = omp_get_wtime();
		block_LU_decomposition(matr, bs);
		t2 = omp_get_wtime();
		std::cout << "block LU, time = " << (double)(t2 - t1) << std::endl;
		print_element(matr);
		double t_LU_bl = double(t2 - t1);

		//separation_LU(matr, L, U);
		//comp_matr_zero(matr);
		//matr_product(L, U, matr);
		//sub_matr(matr, matr_copy);
		//max = max_matr(matr);
		//std::cout << "max (matr - L * U) = " << max << "\n";

		std::cout << std::endl;
		read_matr(file_name, matr);
		t1 = omp_get_wtime();
		block_LU_decomposition_parallel(matr, bs, 4);
		t2 = omp_get_wtime();
		std::cout << "\nblock LU parallel, time = " << (double)(t2 - t1) << "\n";
		print_element(matr);
		double t_LU_bl_par = double(t2 - t1);

		std::cout << "time seq / time par = " << t_LU_bl / t_LU_bl_par << "\n";

		//separation_LU(matr, L, U);
		//comp_matr_zero(matr);
		//matr_product(L, U, matr);
		//sub_matr(matr, matr_copy);
		//max = max_matr(matr);
		//std::cout << "max (matr - L * U) = " << max << "\n";

		if (best_time > t_LU_bl) {
			best_time = t_LU_bl;
			best_bs = bs;
		}
		if (best_time_parallel > t_LU_bl_par) {
			best_time_parallel = t_LU_bl_par;
			best_bs_parallel = bs;
		}
	}

	std::cout << "\nbest bs seq = " << best_bs << "; best time seq = " << best_time << "\n";
	std::cout << "best bs par = " << best_bs_parallel << "; best time par = " << best_time_parallel << "\n";
	std::cout << "time seq bl/ time bl par = " << best_time / best_time_parallel << "\n";
	std::cout << "============================\n " << std::endl;
	std::cout << "time seq / time bl = " << t_LU / best_time << "\n";
	std::cout << "time seq par/ time bl par= " << t_LU_par / best_time_parallel << "\n";

	//std::cout << "\n";
	//read_matr(file_name, matr);
	//t1 = omp_get_wtime();
	//block_LU_decomposition(matr, 16);
	//t2 = omp_get_wtime();
	//std::cout << "block LU, time = " << (double)(t2 - t1) << "\n";
	//double tBl = double(t2 - t1);
	//std::cout << "time / time block = " << t / tBl << "\n";
	////std::cout << "matr[n/2 + 1, n/2] = " << matr[(N / 2 + 1) * N + N / 2] << "\n";

	//separation_LU(matr, L, U);
	//comp_matr_zero(matr);
	//matr_product(L, U, matr);
	//sub_matr(matr, matr_copy);
	//max = max_matr(matr);
	//std::cout << "max (matr - L * U) = " << max << "\n";

	//std::cout << "\n";
	//read_matr(file_name, matr);
	//t1 = omp_get_wtime();
	//block_LU_decomposition_parallel(matr, 64, 4);
	//t2 = omp_get_wtime();
	//std::cout << "block LU parallel, time = " << (double)(t2 - t1) << "\n";
	//double tBlPar = double(t2 - t1);
	//std::cout << "\n";
	//std::cout << "time / time block parallel = " << t / tBlPar << "\n";
	//std::cout << "time block / time block parallel = " << tBl / tBlPar << "\n";
	////std::cout << "matr[n/2 + 1, n/2] = " << matr[(N / 2 + 1) * N + N / 2] << "\n";

	//separation_LU(matr, L, U);
	//comp_matr_zero(matr);
	//matr_product(L, U, matr);
	//sub_matr(matr, matr_copy);
	//max = max_matr(matr);
	//std::cout << "max (matr - L * U) = " << max << "\n";

	clear_mem(matr);
	clear_mem(matr_copy);
	clear_mem(L);
	clear_mem(U);

	return 0;
}

void init_matr(double*& matr) {
	matr = new double[N * N];
}

void comp_matr(double* matr) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			matr[i * N + j] = rand() / 10000.0;
		}
	}
}

void comp_matr_zero(double* matr) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			matr[i * N + j] = 0.0;
		}
	}
}

void clear_mem(double* obj) {
	delete[] obj;
}

void print_matr(double* matr) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			std::cout << std::fixed << std::setprecision(10) << matr[i * N + j] << "\t";
		}
		std::cout << "\n";
	}
}

void print_element(double* matr) {
	std::cout << "matr[n/2][n/2] = " << matr[N / 2 * N + N / 2] << std::endl;
}

void LU_decomposition(double* matr) {
	for (int i = 0; i < N - 1; ++i) {
		for (int j = i + 1; j < N; ++j) {
			double temp = matr[j * N + i] / matr[i * N + i];
			for (int k = i + 1; k < N; ++k) {
				matr[j * N + k] -= temp * matr[i * N + k];
			}
			matr[j * N + i] = temp;
		}
	}
}

void LU_decomposition_parallel(double* matr) {
#pragma omp parallel
	for (int i = 0; i < N - 1; ++i) {
#pragma omp for schedule(static, 1)
		for (int j = i + 1; j < N; ++j) {
			double temp = matr[j * N + i] / matr[i * N + i];
			for (int k = i + 1; k < N; ++k) {
				matr[j * N + k] -= temp * matr[i * N + k];
			}
			matr[j * N + i] = temp;
		}
	}
}

void block_LU_decomposition(double* matr, const int bs) {
	// bs  Размер блока
	double* A11 = nullptr, * U12 = nullptr, * L21 = nullptr; // Диаг., U_12, L_21 блоки соответственно
	A11 = new double[bs * bs];
	U12 = new double[bs * (N - bs)];
	L21 = new double[(N - bs) * bs];

	for (int bi = 0; bi < N - 1; bi += bs) {
		double temp; // Временная переменная
		// Копирование диагонального блока
		for (int i = 0; i < bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				A11[i * bs + j] = matr[(i + bi) * N + (j + bi)];
			}
		}
		// Копирование блока U_12
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				U12[i * bs + j] = matr[(j + bi) * N + (i + bi + bs)];
			}
		}
		// Копирование блока L_21
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				L21[i * bs + j] = matr[(i + bi + bs) * N + (j + bi)];
			}
		}
		// LU разложение диагонального блока
		for (int i = 0; i < bs - 1; ++i) {
			for (int j = i + 1; j < bs; ++j) {
				temp = A11[j * bs + i] / A11[i * bs + i];
				for (int k = i + 1; k < bs; ++k) {
					A11[j * bs + k] = A11[j * bs + k] - temp * A11[i * bs + k];
				}
				A11[j * bs + i] = temp;
			}
		}
		// Заполнение блока U_12
		for (int j = 0; j < N - bi - bs; ++j) {
			for (int i = 1; i < bs; ++i) {
				temp = 0.0;
				for (int k = 0; k <= i - 1; ++k) {
					temp += A11[i * bs + k] * U12[j * bs + k];
				}
				U12[j * bs + i] -= temp;
			}
		}
		// Заполнение блока L_21
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				temp = 0.0;
				for (int k = 0; k <= j - 1; ++k) {
					temp += L21[i * bs + k] * A11[k * bs + j];
				}
				L21[i * bs + j] = (L21[i * bs + j] - temp) / A11[j * bs + j];
			}
		}
		// Вычитание перед рекурсией
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < N - bi - bs; ++j) {
				temp = 0.0;
				for (int k = 0; k < bs; ++k) {
					temp += L21[i * bs + k] * U12[j * bs + k];
				}
				matr[(i + bi + bs) * N + (j + bi + bs)] -= temp;
			}
		}
		// Перенос лок. массивов в матрицу
		// Диаг. блок
		for (int i = 0; i < bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				matr[(i + bi) * N + (j + bi)] = A11[i * bs + j];
			}
		}
		// Блок U_12
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				matr[(j + bi) * N + (i + bi + bs)] = U12[i * bs + j];
			}
		}
		// Блок L_21
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				matr[(i + bi + bs) * N + (j + bi)] = L21[i * bs + j];
			}
		}
	}

	clear_mem(A11);
	clear_mem(U12);
	clear_mem(L21);
}

void block_LU_decomposition_parallel(double* matr, const int bs, int num_th) {
	// bs  Размер блока
	double* A11 = nullptr, * U12 = nullptr, * L21 = nullptr; // Диаг., U_12, L_21 блоки соответственно
	A11 = new double[bs * bs];
	U12 = new double[bs * (N - bs)];
	L21 = new double[(N - bs) * bs];

	for (int bi = 0; bi < N - 1; bi += bs) {
		double temp; // Временная переменная
		// Копирование диагонального блока
		for (int i = 0; i < bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				A11[i * bs + j] = matr[(i + bi) * N + (j + bi)];
			}
		}
		// Копирование блока U_12
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				U12[i * bs + j] = matr[(j + bi) * N + (i + bi + bs)];
			}
		}

		// Копирование блока L_21
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				L21[i * bs + j] = matr[(i + bi + bs) * N + (j + bi)];
			}
		}

		//std::cout << "num_th = " << num_th << std::endl;

		// LU разложение диагонального блока
#pragma omp parallel if(((bs * bs + 1) * bs / 2) > 4000)
		for (int i = 0; i < bs - 1; ++i) {
#pragma omp for shedule(static, 1) private(temp)
			for (int j = i + 1; j < bs; ++j) {
				temp = A11[j * bs + i] / A11[i * bs + i];
				for (int k = i + 1; k < bs; ++k) {
					A11[j * bs + k] -= - temp * A11[i * bs + k];
				}
				A11[j * bs + i] = temp;

				//int max_th = omp_get_max_threads();
				//std::cout << "max_th = " << max_th << std::endl;
				//int now_th = omp_get_num_threads();
				//std::cout << "now_th = " << now_th << std::endl;
			}
		}
		// Заполнение блока U_12
#pragma omp parallel if(((N - bi - bs) * (bs + 1) * bs / 2) > 4000)
		for (int j = 0; j < N - bi - bs; ++j) {
#pragma omp for shedule(static, 1) private(temp)
			for (int i = 1; i < bs; ++i) {
				temp = 0.0;
				for (int k = 0; k <= i - 1; ++k) {
					temp += A11[i * bs + k] * U12[j * bs + k];
				}
				U12[j * bs + i] -= temp;
			}
		}

		// Заполнение блока L_21
#pragma omp parallel if(((N - bi - bs) * (bs + 1) * bs / 2) > 4000)
		for (int i = 0; i < N - bi - bs; ++i) {
#pragma omp for shedule(static, 1) private(temp)
			for (int j = 0; j < bs; ++j) {
				temp = 0.0;
				for (int k = 0; k <= j - 1; ++k) {
					temp += L21[i * bs + k] * A11[k * bs + j];
				}
				L21[i * bs + j] = (L21[i * bs + j] - temp) / A11[j * bs + j];
			}
		}

		// Вычисление A22~
#pragma omp parallel if(((N - bi - bs) * (N - bi - bs) * bs) > 4000)
		for (int i = 0; i < N - bi - bs; ++i) {
#pragma omp for shedule(static, 1) private(temp)
			for (int j = 0; j < N - bi - bs; ++j) {
				temp = 0.0;
				for (int k = 0; k < bs; ++k) {
					temp += L21[i * bs + k] * U12[j * bs + k];
				}
				matr[(i + bi + bs) * N + (j + bi + bs)] -= temp;
			}
		}

		// Перенос лок. массивов в матрицу
		// Диаг. блок
		for (int i = 0; i < bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				matr[(i + bi) * N + (j + bi)] = A11[i * bs + j];
			}
		}
		// Блок U_12
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				matr[(j + bi) * N + (i + bi + bs)] = U12[i * bs + j];
			}
		}
		// Блок L_21
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				matr[(i + bi + bs) * N + (j + bi)] = L21[i * bs + j];
			}
		}
	}

	clear_mem(A11);
	clear_mem(U12);
	clear_mem(L21);
}

void block_LU_decomposition_un(double* matr, const int bs, int num_th) {
	// bs  Размер блока
	double* A11 = nullptr, * U12 = nullptr, * L21 = nullptr; // Диаг., U_12, L_21 блоки соответственно
	A11 = new double[bs * bs];
	U12 = new double[bs * (N - bs)];
	L21 = new double[(N - bs) * bs];

	for (int bi = 0; bi < N - 1; bi += bs) {
		double temp; // Временная переменная
		// Копирование диагонального блока
		for (int i = 0; i < bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				A11[i * bs + j] = matr[(i + bi) * N + (j + bi)];
			}
		}
		// Копирование блока U_12
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				U12[i * bs + j] = matr[(j + bi) * N + (i + bi + bs)];
			}
		}
		// Копирование блока L_21
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				L21[i * bs + j] = matr[(i + bi + bs) * N + (j + bi)];
			}
		}
		// LU разложение диагонального блока
		for (int i = 0; i < bs - 1; ++i) {
#pragma omp parallel for private(temp) if(num_th > 1)//((bs * bs + 1) * bs / 2) > 4000)
			for (int j = i + 1; j < bs; ++j) {
				temp = A11[j * bs + i] / A11[i * bs + i];
				for (int k = i + 1; k < bs; ++k) {
					A11[j * bs + k] = A11[j * bs + k] - temp * A11[i * bs + k];
				}
				A11[j * bs + i] = temp;
			}
		}
		// Заполнение блока U_12
#pragma omp parallel for private(temp) if(num_th > 1)//((N - bi - bs) * (bs + 1) * bs / 2) > 4000)
		for (int j = 0; j < N - bi - bs; ++j) {
			for (int i = 1; i < bs; ++i) {
				temp = 0.0;
				for (int k = 0; k <= i - 1; ++k) {
					temp += A11[i * bs + k] * U12[j * bs + k];
				}
				U12[j * bs + i] -= temp;
			}
		}
		// Заполнение блока L_21
#pragma omp parallel for private(temp) if(num_th > 1)//((N - bi - bs) * (bs + 1) * bs / 2) > 4000)
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				temp = 0.0;
				for (int k = 0; k <= j - 1; ++k) {
					temp += L21[i * bs + k] * A11[k * bs + j];
				}
				L21[i * bs + j] = (L21[i * bs + j] - temp) / A11[j * bs + j];
			}
		}
		// Вычисление A22~
#pragma omp parallel for private(temp) if(num_th > 1)//((N - bi - bs) * (N - bi - bs) * bs) > 4000)
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < N - bi - bs; ++j) {
				temp = 0.0;
				for (int k = 0; k < bs; ++k) {
					temp += L21[i * bs + k] * U12[j * bs + k];
				}
				matr[(i + bi + bs) * N + (j + bi + bs)] -= temp;
			}
		}

		// Перенос лок. массивов в матрицу
		// Диаг. блок
		for (int i = 0; i < bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				matr[(i + bi) * N + (j + bi)] = A11[i * bs + j];
			}
		}
		// Блок U_12
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				matr[(j + bi) * N + (i + bi + bs)] = U12[i * bs + j];
			}
		}
		// Блок L_21
		for (int i = 0; i < N - bi - bs; ++i) {
			for (int j = 0; j < bs; ++j) {
				matr[(i + bi + bs) * N + (j + bi)] = L21[i * bs + j];
			}
		}
	}

	clear_mem(A11);
	clear_mem(U12);
	clear_mem(L21);
}

void matr_product(double* A, double* B, double* C) {
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		for (int k = 0; k < N; ++k) {
			for (int j = 0; j < N; ++j) {
				C[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}
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

void matr_equal(double* result_matr, double* source_matr) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			result_matr[i * N + j] = source_matr[i * N + j];
		}
	}
}

void sub_matr(double* result_matr, double* source_matr) {
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			result_matr[i * N + j] -= source_matr[i * N + j];
		}
	}
}

double max_matr(double* matr) {
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

void read_matr(std::string namefile, double*& matr) {
	std::ifstream fin;
	fin.open(namefile.c_str(), std::ios::in);
	if (!fin.is_open()) {
		std::cout << "Error file!\n";
	}
	else {
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				fin >> matr[i * N + j];
			}
		}
	}
	fin.close();
}

void write_matr(std::string namefile, double*& matr) {
	std::ofstream fout;
	fout.open(namefile.c_str());
	if (!fout.is_open()) {
		std::cout << "Error file!\n";
	}
	else {
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				fout << matr[i * N + j] << "\t";
			}
			fout << "\n";
		}
	}
	fout.close();
	std::cout << "write matr complete\n";
}
