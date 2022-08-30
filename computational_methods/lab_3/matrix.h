#pragma once

#include "constants.h"


double computeMatrixDet(const MyType* const* const A, const int n)  // Вычисление определителя
{
    MyType** A_copy = new MyType*[n];
    for (int i = 0; i < n; ++i) {
        A_copy[i] = new MyType[n];
        for (int j = 0; j < n; ++j) {
            A_copy[i][j] = A[i][j];
        }
    }

    double det = 1.0;
    for (int k = 0; k < n; ++k) {
        int iMax = k;
        for (int i = k + 1; i < n; ++i) {
            if (fabs(A_copy[i][k]) > fabs(A_copy[iMax][k])) {
                iMax = i;
            }
        }
        if (fabs(A_copy[iMax][k]) < EPS0) {
            return 0.0;
        }
        if (iMax != k) {
            MyType* temp = A_copy[k];
            A_copy[k] = A_copy[iMax];
            A_copy[iMax] = temp;
            det *= -1;
        }

        det *= A_copy[k][k];
        for (int i = k + 1; i < n; ++i) {
            double t = A_copy[i][k] / A_copy[k][k];
            for (int j = k; j < n; ++j) {
                A_copy[i][j] -= A_copy[k][j] * t;
            }
        }
    }
    return det;
}


void multiplyMatrix(const MyType* const* const A, const MyType* const* const B,
                    MyType* const* const result, const int n)  // Перемножение матриц
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            MyType sum = 0.0;
            for (int k = 0; k < n; ++k) {
                sum += A[i][k] * B[k][j];
            }
            result[i][j] = sum;
        }
    }
}


void transpose(MyType* const* const A, const int n)  // Транспонирование на месте
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            MyType temp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = temp;
        }
    }
}


void transposed(const MyType* const* const A, MyType* const* const result,
                const int n)  // Транспонирование
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = A[j][i];
        }
    }
}


int qrDecomposition(const MyType* const* const A, MyType* const* const Q,
                    MyType* const* const R, const int n)  // QR-разложение
{
    // Инициализируем матрицы
    // T == tranposed(Q)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Q[i][j] = (i == j) ? 1.0 : 0.0;
            R[i][j] = A[i][j];  // R = T * A
        }
    }

    int iter_count = 0;
    // Далее: k - диагональ, i - строка, j - столбец
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            if (fabs(R[i][k]) < EPS0) {
                continue;  // Вращение не требуется
            }

            MyType denominator = sqrt(R[k][k] * R[k][k] + R[i][k] * R[i][k]);

            MyType c = R[k][k] / denominator;
            MyType s = R[i][k] / denominator;
            R[k][k] = denominator;  // Диагональный элемент будет равен знаменателю
            R[i][k] = 0.0;          // А этот должен обнулиться
            iter_count += 5;

            // Домножаем матрицу A слева на матрицу T_ij
            for (int j = k + 1; j < n; ++j) {
                MyType Rkj = R[k][j];
                MyType Rij = R[i][j];
                R[k][j] = c * Rkj + s * Rij;
                R[i][j] = -s * Rkj + c * Rij;
                iter_count += 4;
            }

            // Домножаем матрицу T слева на матрицу T_ij
            for (int j = 0; j < n; ++j) {
                MyType Tkj = Q[k][j];
                MyType Tij = Q[i][j];
                Q[k][j] = c * Tkj + s * Tij;
                Q[i][j] = -s * Tkj + c * Tij;
                iter_count += 4;
            }
        }
    }

    transpose(Q, n);
    return iter_count;
}

int qrDecompositionHessenberg(const MyType* const* const A, MyType* const* const Q,
                               MyType* const* const R, const int n)  // QR-разложение
{
    // Инициализируем матрицы
    // T == tranposed(Q)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Q[i][j] = (i == j) ? 1.0 : 0.0;
            R[i][j] = A[i][j];  // R = T * A
        }
    }

    int iter_count = 0;
    // Далее: k - диагональ, i - строка, j - столбец
    for (int k = 0; k < n - 1; ++k) {
        int i = k + 1;

        if (fabs(R[i][k]) < EPS0) {
            continue;  // Вращение не требуется
        }

        MyType denominator = sqrt(R[k][k] * R[k][k] + R[i][k] * R[i][k]);

        MyType c = R[k][k] / denominator;
        MyType s = R[i][k] / denominator;
        R[k][k] = denominator;  // Диагональный элемент будет равен знаменателю
        iter_count += 5;
        R[i][k] = 0.0;          // А этот должен обнулиться

        // Домножаем матрицу A слева на матрицу T_ij
        for (int j = k + 1; j < fmin(k + 3, n); ++j) {
            MyType Rkj = R[k][j];
            MyType Rij = R[i][j];
            R[k][j] = c * Rkj + s * Rij;
            R[i][j] = -s * Rkj + c * Rij;
            iter_count += 4;
        }

        // Домножаем матрицу T слева на матрицу T_ij
        for (int j = 0; j < n; ++j) {
            MyType Tkj = Q[k][j];
            MyType Tij = Q[i][j];
            Q[k][j] = c * Tkj + s * Tij;
            Q[i][j] = -s * Tkj + c * Tij;
            iter_count += 4;
        }
    }

    transpose(Q, n);
    return iter_count;
}


void matrixMulByVector(const MyType* const* const A, const MyType* const b,
                       MyType* const result, int n)  // Умножение матрицы на вектор
{
    for (int i = 0; i < n; ++i) {
        MyType sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A[i][j] * b[j];
        }
        result[i] = sum;
    }
}