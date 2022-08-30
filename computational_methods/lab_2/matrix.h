#pragma once

#include "constants.h"


void transpose(MyType* const* const A, const int n)  // Транспонирование на месте
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n / 2; ++j) {
            MyType temp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = temp;
        }
    }
}


void transposed(const MyType* const* const A, MyType* const* const B,
                const int n)  // Транспонирование
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            B[i][j] = A[j][i];
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

    // Далее: k - диагональ, i - строка, j - столбец
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            if (R[i][k] == 0) {
                continue;  // Вращение не требуется
            }

            MyType denominator = sqrt(R[k][k] * R[k][k] + R[i][k] * R[i][k]);

            MyType c = R[k][k] / denominator;
            MyType s = R[i][k] / denominator;
            R[k][k] = denominator;  // Диагональный элемент будет равен знаменателю
            R[i][k] = 0.0;          // А этот должен обнулиться

            // Домножаем матрицу A слева на матрицу T_ij
            for (int j = k + 1; j < n; ++j) {
                MyType Rkj = R[k][j];
                MyType Rij = R[i][j];
                R[k][j] = c * Rkj + s * Rij;
                R[i][j] = -s * Rkj + c * Rij;
            }

            // Домножаем матрицу T слева на матрицу T_ij
            for (int j = 0; j < n; ++j) {
                MyType Tkj = Q[k][j];
                MyType Tij = Q[i][j];
                Q[k][j] = c * Tkj + s * Tij;
                Q[i][j] = -s * Tkj + c * Tij;
            }
        }
    }

    transpose(Q, n);

    return 0;
}