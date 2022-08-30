#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "omp.h"

const int dim = 1024;

void init_matr(double *&matr);
void comp_matr(double *matr);
void comp_matr_zero(double *matr);
void clear_mem(double *obj);
void print_matr(double *matr);
void LU_decomposition(double *matr);
void block_LU_decomposition(double *matr, const int bs);
void block_LU_decomposition_parallel(double* matr, const int bs, int num_th);
void block_LU_decomposition_un(double *matr, const int bs, int num_th);
void matr_product(double *A, double *B, double *C);
void separation_LU(double *matr, double *L, double *U);
void matr_equal(double *result_matr, double *source_matr);
void sub_matr(double *result_matr, double *source_matr);
double max_matr(double *matr);
void read_matr(std::string namefile, double*& matr);
void write_matr(std::string namefile, double*& matr);

int main() {
    double *matr = nullptr, *matr_copy = nullptr, *L = nullptr, *U = nullptr;

    init_matr(matr);
    init_matr(matr_copy);
    init_matr(L);
    init_matr(U);

    std::string file_name = "rand_matr4096.txt";

    comp_matr(matr);
    //write_matr(file_name, matr);

    // read_matr(file_name, matr);
    //read_matr(file_name, matr_copy);

    double t1, t2;
    double max;

    t1 = omp_get_wtime();
    LU_decomposition(matr);
    t2 = omp_get_wtime();
    std::cout << "LU, time = " << (double)(t2 - t1) << "\n";
    double t_seq = double(t2 - t1);
    //std::cout << "matr[n/2 + 1, n/2] = " << matr[(dim / 2 + 1) * dim + dim / 2] << "\n";

    //separation_LU(matr, L, U);
    //comp_matr_zero(matr);
    //matr_product(L, U, matr);
    //sub_matr(matr, matr_copy);
    //max = max_matr(matr);
    //std::cout << "max (matr - L * U) = " << max << "\n";

    int best_bs = 0, best_bs_parallel = 0;
    double best_time = 1000.0, best_time_parallel = 1000.0;

    for (int bs = 64; bs <= 65; bs *= 2) {
    //int bs = 64;
        std::cout << "============================\nbs = " << bs << "\n";
        read_matr(file_name, matr);
        t1 = omp_get_wtime();
        block_LU_decomposition(matr, bs);
        t2 = omp_get_wtime();
        std::cout << "block LU, time = " << (double)(t2 - t1) << "\n";
        double tBl = double(t2 - t1);

        //separation_LU(matr, L, U);
        //comp_matr_zero(matr);
        //matr_product(L, U, matr);
        //sub_matr(matr, matr_copy);
        //max = max_matr(matr);
        //std::cout << "max (matr - L * U) = " << max << "\n";

        std::cout << "\n";
        read_matr(file_name, matr);
        t1 = omp_get_wtime();
        block_LU_decomposition_parallel(matr, bs, 2);
        t2 = omp_get_wtime();
        std::cout << "block LU parallel, time = " << (double)(t2 - t1) << "\n";
        double tBlPar = double(t2 - t1);

        //separation_LU(matr, L, U);
        //comp_matr_zero(matr);
        //matr_product(L, U, matr);
        //sub_matr(matr, matr_copy);
        //max = max_matr(matr);
        //std::cout << "max (matr - L * U) = " << max << "\n";

        if (best_time > tBl) {
            best_time = tBl;
            best_bs = bs;
        }
        if (best_time_parallel > tBlPar) {
            best_time_parallel = tBlPar;
            best_bs_parallel = bs;
        }
    }

    std::cout << "\n\ntime seq = " << t_seq << "\n";
    std::cout << "best bs seq = " << best_bs << " best time seq = " << best_time << "\n";
    std::cout << "best bs par = " << best_bs_parallel << " best time par = " << best_time_parallel << "\n";
    std::cout << "time seq / time par = " << best_time / best_time_parallel << "\n";

    //std::cout << "\n";
    //read_matr(file_name, matr);
    //t1 = omp_get_wtime();
    //block_LU_decomposition(matr, 16);
    //t2 = omp_get_wtime();
    //std::cout << "block LU, time = " << (double)(t2 - t1) << "\n";
    //double tBl = double(t2 - t1);
    //std::cout << "time / time block = " << t / tBl << "\n";
    ////std::cout << "matr[n/2 + 1, n/2] = " << matr[(dim / 2 + 1) * dim + dim / 2] << "\n";

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
    ////std::cout << "matr[n/2 + 1, n/2] = " << matr[(dim / 2 + 1) * dim + dim / 2] << "\n";

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

void init_matr(double *&matr) {
    matr = new double[dim * dim];
}

void comp_matr(double *matr) {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            matr[i * dim + j] = rand() % dim + 1.0/*cos(i + j) + 2.0*//*int(cos(i + j) * 10)*/;
        }
    }
}

void comp_matr_zero(double *matr) {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            matr[i * dim + j] = 0.0;
        }
    }
}

void clear_mem(double *obj) {
    delete[] obj;
}

void print_matr(double *matr) {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            std::cout << std::fixed << std::setprecision(10) << matr[i * dim + j] << "\t";
        }
        std::cout << "\n";
    }
}

void LU_decomposition(double *matr) {
    for (int i = 0; i < dim - 1; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            double temp = matr[j * dim + i] / matr[i * dim + i];
            for (int k = i + 1; k < dim; ++k) {
                matr[j * dim + k] -= temp * matr[i * dim + k];
            }
            matr[j * dim + i] = temp;
        }
    }
}

void block_LU_decomposition(double *matr, const int bs) {

    // bs  ������ �����

    double *A11 = nullptr, *U12 = nullptr, *L21 = nullptr; // ����., U_12, L_21 ����� ��������������
    A11 = new double[bs * bs];
    U12 = new double[bs * (dim - bs)];
    L21 = new double[(dim - bs) * bs];

    for (int bi = 0; bi < dim - 1; bi += bs) {
        double temp; // ��������� ����������
        // ����������� ������������� �����
        for (int i = 0; i < bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                A11[i * bs + j] = matr[(i + bi) * dim + (j + bi)];
            }
        }
        // ����������� ����� U_12
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                U12[i * bs + j] = matr[(j + bi) * dim + (i + bi + bs)];
            }
        }
        // ����������� ����� L_21
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                L21[i * bs + j] = matr[(i + bi + bs) * dim + (j + bi)];
            }
        }
        // LU ���������� ������������� �����
        for (int i = 0; i < bs - 1; ++i) {
            for (int j = i + 1; j < bs; ++j) {
                temp = A11[j * bs + i] / A11[i * bs + i];
                for (int k = i + 1; k < bs; ++k) {
                    A11[j * bs + k] = A11[j * bs + k] - temp * A11[i * bs + k];
                }
                A11[j * bs + i] = temp;
            }
        }
        // ���������� ����� U_12
        for (int j = 0; j < dim - bi - bs; ++j) {
            for (int i = 1; i < bs; ++i) {
                temp = 0.0;
                for (int k = 0; k <= i - 1; ++k) {
                    temp += A11[i * bs + k] * U12[j * bs + k];
                }
                U12[j * bs + i] -= temp;
            }
        }
        // ���������� ����� L_21
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                temp = 0.0;
                for (int k = 0; k <= j - 1; ++k) {
                    temp += L21[i * bs + k] * A11[k * bs + j];
                }
                L21[i * bs + j] = (L21[i * bs + j] - temp) / A11[j * bs + j];
            }
        }
        // ��������� ����� ���������
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < dim - bi - bs; ++j) {
                temp = 0.0;
                for (int k = 0; k < bs; ++k) {
                     temp += L21[i * bs + k] * U12[j * bs + k];
                }
                matr[(i + bi + bs) * dim + (j + bi + bs)] -= temp;
            }
        }
        // ������� ���. �������� � �������
        // ����. ����
        for (int i = 0; i < bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(i + bi) * dim + (j + bi)] = A11[i * bs + j];
            }
        }
        // ���� U_12
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(j + bi) * dim + (i + bi + bs)] = U12[i * bs + j];
            }
        }
        // ���� L_21
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(i + bi + bs) * dim + (j + bi)] = L21[i * bs + j];
            }
        }
    }

    clear_mem(A11);
    clear_mem(U12);
    clear_mem(L21);
}

void block_LU_decomposition_parallel(double* matr, const int bs, int num_th) {
    // bs  ������ �����
    double* A11 = nullptr, * U12 = nullptr, * L21 = nullptr; // ����., U_12, L_21 ����� ��������������
    A11 = new double[bs * bs];
    U12 = new double[bs * (dim - bs)];
    L21 = new double[(dim - bs) * bs];

    for (int bi = 0; bi < dim - 1; bi += bs) {
        double temp; // ��������� ����������
        // ����������� ������������� �����
        for (int i = 0; i < bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                A11[i * bs + j] = matr[(i + bi) * dim + (j + bi)];
            }
        }
        // ����������� ����� U_12
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                U12[i * bs + j] = matr[(j + bi) * dim + (i + bi + bs)];
            }
        }
        // ����������� ����� L_21
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                L21[i * bs + j] = matr[(i + bi + bs) * dim + (j + bi)];
            }
        }
        // LU ���������� ������������� �����
        for (int i = 0; i < bs - 1; ++i) {
#pragma omp parallel for private(temp) num_threads(num_th) if(((bs * bs + 1) * bs / 2) > 4000)
            for (int j = i + 1; j < bs; ++j) {
                temp = A11[j * bs + i] / A11[i * bs + i];
                for (int k = i + 1; k < bs; ++k) {
                    A11[j * bs + k] = A11[j * bs + k] - temp * A11[i * bs + k];
                }
                A11[j * bs + i] = temp;
            }
        }
        // ���������� ����� U_12
#pragma omp parallel for private(temp) num_threads(num_th) if(((dim - bi - bs) * (bs + 1) * bs / 2) > 4000)
        for (int j = 0; j < dim - bi - bs; ++j) {
            for (int i = 1; i < bs; ++i) {
                temp = 0.0;
                for (int k = 0; k <= i - 1; ++k) {
                    temp += A11[i * bs + k] * U12[j * bs + k];
                }
                U12[j * bs + i] -= temp;
            }
        }
        // ���������� ����� L_21
#pragma omp parallel for private(temp) num_threads(num_th) if(((dim - bi - bs) * (bs + 1) * bs / 2) > 4000)
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                temp = 0.0;
                for (int k = 0; k <= j - 1; ++k) {
                    temp += L21[i * bs + k] * A11[k * bs + j];
                }
                L21[i * bs + j] = (L21[i * bs + j] - temp) / A11[j * bs + j];
            }
        }
        // ���������� A22~
#pragma omp parallel for private(temp) num_threads(num_th) if(((dim - bi - bs) * (dim - bi - bs) * bs) > 4000)
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < dim - bi - bs; ++j) {
                temp = 0.0;
                for (int k = 0; k < bs; ++k) {
                    temp += L21[i * bs + k] * U12[j * bs + k];
                }
                matr[(i + bi + bs) * dim + (j + bi + bs)] -= temp;
            }
        }

        // ������� ���. �������� � �������
        // ����. ����
        for (int i = 0; i < bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(i + bi) * dim + (j + bi)] = A11[i * bs + j];
            }
        }
        // ���� U_12
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(j + bi) * dim + (i + bi + bs)] = U12[i * bs + j];
            }
        }
        // ���� L_21
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(i + bi + bs) * dim + (j + bi)] = L21[i * bs + j];
            }
        }
    }

    clear_mem(A11);
    clear_mem(U12);
    clear_mem(L21);
}

void block_LU_decomposition_un(double *matr, const int bs, int num_th) {
    // bs  ������ �����
    double *A11 = nullptr, *U12 = nullptr, *L21 = nullptr; // ����., U_12, L_21 ����� ��������������
    A11 = new double[bs * bs];
    U12 = new double[bs * (dim - bs)];
    L21 = new double[(dim - bs) * bs];

    for (int bi = 0; bi < dim - 1; bi += bs) {
        double temp; // ��������� ����������
        // ����������� ������������� �����
        for (int i = 0; i < bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                A11[i * bs + j] = matr[(i + bi) * dim + (j + bi)];
            }
        }
        // ����������� ����� U_12
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                U12[i * bs + j] = matr[(j + bi) * dim + (i + bi + bs)];
            }
        }
        // ����������� ����� L_21
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                L21[i * bs + j] = matr[(i + bi + bs) * dim + (j + bi)];
            }
        }
        // LU ���������� ������������� �����
        for (int i = 0; i < bs - 1; ++i) {
#pragma omp parallel for private(temp) num_threads(num_th) if(num_th > 1)//((bs * bs + 1) * bs / 2) > 4000)
            for (int j = i + 1; j < bs; ++j) {
                temp = A11[j * bs + i] / A11[i * bs + i];
                for (int k = i + 1; k < bs; ++k) {
                    A11[j * bs + k] = A11[j * bs + k] - temp * A11[i * bs + k];
                }
                A11[j * bs + i] = temp;
            }
        }
        // ���������� ����� U_12
#pragma omp parallel for private(temp) num_threads(num_th) if(num_th > 1)//((dim - bi - bs) * (bs + 1) * bs / 2) > 4000)
        for (int j = 0; j < dim - bi - bs; ++j) {
            for (int i = 1; i < bs; ++i) {
                temp = 0.0;
                for (int k = 0; k <= i - 1; ++k) {
                    temp += A11[i * bs + k] * U12[j * bs + k];
                }
                U12[j * bs + i] -= temp;
            }
        }
        // ���������� ����� L_21
#pragma omp parallel for private(temp) num_threads(num_th) if(num_th > 1)//((dim - bi - bs) * (bs + 1) * bs / 2) > 4000)
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                temp = 0.0;
                for (int k = 0; k <= j - 1; ++k) {
                    temp += L21[i * bs + k] * A11[k * bs + j];
                }
                L21[i * bs + j] = (L21[i * bs + j] - temp) / A11[j * bs + j];
            }
        }
        // ���������� A22~
#pragma omp parallel for private(temp) num_threads(num_th) if(num_th > 1)//((dim - bi - bs) * (dim - bi - bs) * bs) > 4000)
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < dim - bi - bs; ++j) {
                temp = 0.0;
                for (int k = 0; k < bs; ++k) {
                    temp += L21[i * bs + k] * U12[j * bs + k];
                }
                matr[(i + bi + bs) * dim + (j + bi + bs)] -= temp;
            }
        }

        // ������� ���. �������� � �������
        // ����. ����
        for (int i = 0; i < bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(i + bi) * dim + (j + bi)] = A11[i * bs + j];
            }
        }
        // ���� U_12
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(j + bi) * dim + (i + bi + bs)] = U12[i * bs + j];
            }
        }
        // ���� L_21
        for (int i = 0; i < dim - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(i + bi + bs) * dim + (j + bi)] = L21[i * bs + j];
            }
        }
    }

    clear_mem(A11);
    clear_mem(U12);
    clear_mem(L21);
}

void matr_product(double *A, double *B, double *C) {
#pragma omp parallel for
    for (int i = 0; i < dim; ++i) {
        for (int k = 0; k < dim; ++k) {
            for (int j = 0; j < dim; ++j) {
                C[i * dim + j] += A[i * dim + k] * B[k * dim + j];
            }
        }
    }
}

void separation_LU(double *matr, double *L, double *U) {
#pragma omp parallel for
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j <= i; ++j) {
            L[i * dim + j] = matr[i * dim + j];
            U[j * dim + i] = matr[j * dim + i];
        }
        L[i * dim + i] = 1.0;
    }
}

void matr_equal(double *result_matr, double *source_matr) {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            result_matr[i * dim + j] = source_matr[i * dim + j];
        }
    }
}

void sub_matr(double *result_matr, double *source_matr) {
#pragma omp parallel for
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            result_matr[i * dim + j] -= source_matr[i * dim + j];
        }
    }
}

double max_matr(double *matr) {
    double max = 0.0;
#pragma omp parallel for
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            double temp = fabs(matr[i * dim + j]);
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
    } else {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                fin >> matr[i * dim + j];
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
    } else {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                fout << matr[i * dim + j] << "\t";
            }
            fout << "\n";
        }
    }
    fout.close();
    std::cout << "write matr complete\n";
}
